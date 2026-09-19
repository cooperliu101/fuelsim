"""Manual cylinder heat-flux accuracy study using unchanged production cards.

The first-order reconstruction is a closed-form discrete analytical solution,
not a second numerical solver. Report raw integration-point flux, never recovered
or smoothed flux in place of the production quantity.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from compare import metric
from study import close, output, run


def inspect(folder, name):
    data = output(folder, name)
    conn = data['connect1']-1
    radius = data['coordx']
    temperature = data['temperature'][-1]
    linear = conn.shape[1] == 4
    exact_temperature = 400-100*np.log(radius/.01)/np.log(3)
    q, exact_q, positions, transverse = [], [], [], []
    for point in range(4 if linear else 9):
        r = data[f'point_x_q{point}'][-1]
        positions.append(r)
        q.append(data[f'heat_flux_x_q{point}'][-1])
        exact_q.append(1000/(r*np.log(3)))
        for direction in 'yz':
            transverse.extend(data[f'heat_flux_{direction}_q{point}'][-1])
    q, exact_q, positions = map(np.asarray, [q, exact_q, positions])
    relative = np.abs(q/exact_q-1)
    index = np.unravel_index(np.argmax(relative), relative.shape)
    heat = np.sum(data['heat_reaction'][-1, np.isclose(radius, .01)])
    exact_heat = 20*np.pi/np.log(3)
    fields = {'temperature': metric(temperature, exact_temperature),
              'temperature_rise': metric(temperature-300, exact_temperature-300),
              'radial_flux': metric(q.flatten(), exact_q.flatten()),
              'transverse_flux': metric(transverse, np.zeros(len(transverse)), zero=True),
              'total_hot_boundary_heat': metric([heat], [exact_heat])}
    result = {'case': name, 'elements': len(conn), 'temperature_unknowns_before_constraints': len(radius),
              'fields': fields, 'all_reported_metrics_below_threshold': all(v['passed'] for v in fields.values()),
              'maximum_flux_error_radius_m': float(positions[index]),
              'maximum_flux_error_element': int(index[1]+1)}
    if linear:
        a, b = np.min(radius[conn], axis=1), np.max(radius[conn], axis=1)
        h, middle = b-a, (a+b)/2
        series = np.sum(h/middle)
        predicted_q = np.broadcast_to(1000/(series*middle), q.shape)
        close(q, predicted_q, 'DCAX4 differs from its discrete analytical flux', atol=1e-6, rtol=1e-9)
        result['discrete_flux_reconstruction_max_relative_error'] = float(np.max(np.abs(q/predicted_q-1)))
        result['discrete_total_heat_relative_error'] = float(abs(heat/(20*np.pi/series)-1))
        result['global_heat_bias'] = float(np.log(3)/series-1)
        result['maximum_element_half_width_over_mean_radius'] = float(np.max(h/(2*middle)))
        # For this source-free cylinder, ∫r*dT/dr is integrated exactly. The
        # sole remaining first-order point-flux error follows from the basis.
        predicted_ratio = np.log(3)/series*positions/middle
        close(q/exact_q, predicted_ratio, 'Flux-error decomposition', atol=1e-9, rtol=0)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fuelsim', type=Path, required=True)
    p.add_argument('--source', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--work', type=Path, required=True)
    args = p.parse_args()
    cases = ['study_dcax4_n16', 'study_dcax4_n64', 'study_dcax4_n128',
             'study_dcax4_geometric_n16', 'study_dcax4_geometric_n64',
             'study_dcax8_n4', 'study_dcax8_n8', 'study_dcax8_n16']
    results = []
    for name in cases:
        folder, log = run(args, name)
        result = inspect(folder, name)
        results.append(result)
        print(name, 'raw flux error (%) =', 100*result['fields']['radial_flux']['maximum_pointwise_relative'])
    accepted = ['study_dcax4_n128', 'study_dcax4_geometric_n64', 'study_dcax8_n8', 'study_dcax8_n16']
    passed = all(r['all_reported_metrics_below_threshold'] for r in results if r['case'] in accepted)
    summary = {'scope': 'Fixed cylindrical geometry, constant conductivity; all nodes and all Gauss points.',
               'relative_threshold': .005, 'zero_flux_absolute_threshold_W_m2': 1e-7,
               'passed': passed, 'accepted_meshes': accepted, 'cases': results}
    (args.work/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if not passed:
        raise AssertionError('At least one proposed mesh misses the unchanged accuracy threshold')


if __name__ == '__main__':
    main()
