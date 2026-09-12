"""Read-only inner-body force/energy audit for the prescribed GPS contact case.

This compares actual radial reactions, axial sectional forces, mean stresses
and elastic strain energy after the mean-hoop formulation change. It uses the
complete ten-step production history and the unchanged native Abaqus reference.
"""
import argparse
import hashlib
import math
from pathlib import Path

import netCDF4

from analyze_contact_stress import COMPONENTS, ROOT, metrics, read_csv


def dot(a, b):
    return sum(x*y*(2 if i == 3 else 1) for i, (x, y) in enumerate(zip(a, b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    nodes = read_csv(ROOT/'gps_two_slice_contact_nodes.csv', ('time', 'node'))
    points = read_csv(ROOT/'gps_two_slice_contact_points.csv', ('time', 'element', 'point'))
    errors = dict.fromkeys(('mean_stress', 'axial_force', 'radial_reaction', 'elastic_energy'), 0.0)
    groups = {name: [] for name in errors}
    samples = 0
    first = {}
    with netCDF4.Dataset(args.results) as data:
        if data['time_whole'][:].tolist() != list(range(11)) or len(nodes) != 160 or len(points) != 160:
            raise RuntimeError('Expected complete ten-step native and production histories')
        names = netCDF4.chartostring(data['name_elem_var'][:]).tolist()
        nodal_names = netCDF4.chartostring(data['name_nod_var'][:]).tolist()

        def field(name, time, layer):
            value = float(data['vals_elem_var%deb1' % (names.index(name)+1)][time, layer])
            if not math.isfinite(value):
                raise RuntimeError('Nonfinite production element field')
            return value

        def reaction(time, layer, end):
            value = float(data['vals_nod_var%d' % (nodal_names.index('reaction_force_r')+1)][time, 4*layer+end])
            if not math.isfinite(value):
                raise RuntimeError('Nonfinite production reaction')
            return value

        for time in range(1, 11):
            for layer in range(2):
                element = 1+2*layer
                n = [nodes[(time, 8*layer+i)] for i in range(1, 5)]
                height = n[3]['z']-n[0]['z']
                u = sum(v['ur'] for v in n)/4
                if max(abs(v['ur']-u) for v in n) > 1e-13 or max(abs(v['uz']-n[0]['uz']) for v in n) > 1e-13:
                    raise RuntimeError('Audit requires constant prescribed displacements within each element')
                volume = 0.0
                gps_integral, native_integral = [0.0]*4, [0.0]*4
                gps_energy, native_energy = 0.0, 0.0
                for q in range(2):
                    suffix = '_q%d' % q
                    radius = field('reference_r'+suffix, time, layer)
                    weight = field('reference_measure'+suffix, time, layer)
                    volume += weight
                    gps_stress = [field('stress_'+c+suffix, time, layer) for c in COMPONENTS]
                    gps_elastic = [field('elastic_'+c+suffix, time, layer) for c in COMPONENTS]
                    gps_energy += 0.5*dot(gps_stress, gps_elastic)*weight
                    for c in range(4):
                        gps_integral[c] += gps_stress[c]*weight
                    native_weight = 0.0
                    for point in (q+1, q+3):
                        p = points[(time, element, point)]
                        if abs(p['r']-radius) > 1e-14:
                            raise RuntimeError('Native and production radial material points differ')
                        native_weight += p['volume']
                        native_stress = [p['stress_'+c] for c in COMPONENTS]
                        native_elastic = [p['elastic_'+c] for c in COMPONENTS]
                        native_energy += 0.5*dot(native_stress, native_elastic)*p['volume']
                        for c in range(4):
                            native_integral[c] += native_stress[c]*p['volume']
                    if abs(native_weight-weight) > 1e-18:
                        raise RuntimeError('Native/production reference measures differ')
                mean_difference = max(abs(gps_integral[c]-native_integral[c])/volume for c in range(4))
                axial_force = field('axial_force', time, layer)
                errors['mean_stress'] = max(errors['mean_stress'], mean_difference)
                errors['axial_force'] = max(errors['axial_force'], abs(axial_force-native_integral[1]/height))
                gps_reaction = [reaction(time, layer, end) for end in range(2)]
                native_reaction = [n[0]['rf_r']+n[3]['rf_r'], n[1]['rf_r']+n[2]['rf_r']]
                delta_reaction = [a-b for a, b in zip(gps_reaction, native_reaction)]
                errors['radial_reaction'] = max(errors['radial_reaction'], max(abs(v) for v in delta_reaction))
                delta_energy = gps_energy-native_energy
                errors['elastic_energy'] = max(errors['elastic_energy'], abs(delta_energy))
                groups['mean_stress'].append(([v/volume for v in gps_integral],
                                               [v/volume for v in native_integral]))
                groups['axial_force'].append(([axial_force], [native_integral[1]/height]))
                for actual, reference in zip(gps_reaction, native_reaction):
                    groups['radial_reaction'].append(([actual], [reference]))
                groups['elastic_energy'].append(([gps_energy], [native_energy]))
                samples += 1
                if time == 1 and layer == 0:
                    first = {'mean_stress': [v/volume for v in gps_integral],
                             'axial_force': axial_force, 'gps_radial_reaction': gps_reaction,
                             'abaqus_radial_reaction': native_reaction, 'radial_reaction_difference': delta_reaction,
                             'gps_elastic_energy': gps_energy, 'abaqus_elastic_energy': native_energy,
                             'elastic_energy_difference': delta_energy,
                             'elastic_energy_relative_difference': delta_energy/native_energy}
    print('scope=current_mean_hoop_production_outputs; native_force_and_energy_comparison')
    print('inner_slice_step_samples=%d' % samples)
    for name, value in first.items():
        values = value if isinstance(value, list) else [value]
        print('first_slice_closed_'+name+'='+','.join('%.17g' % v for v in values))
    print('nonzero_field_metrics=relative_l2,relative_absolute_peak,maximum_pointwise_relative')
    print('metric_units=fraction; multiply by 100 for percent')
    for name, values in groups.items():
        values_metrics = metrics(values)
        print(name+'='+','.join('%.17g' % v for v in values_metrics))
        if max(values_metrics) >= 0.001:
            raise RuntimeError('Native force/energy comparison failed: '+name)
    for name, value in errors.items():
        print(name+'_maximum_absolute_error=%.17g' % value)
        tolerance = 1e-4 if name == 'mean_stress' else (1e-12 if name == 'elastic_energy' else 1e-8)
        if value > tolerance:
            raise RuntimeError('Independent audit failed: '+name)
    if samples != 20:
        raise RuntimeError('Incomplete slice-step coverage')
    if len(groups['radial_reaction']) != 40 or any(len(groups[name]) != 20
            for name in ('mean_stress', 'axial_force', 'elastic_energy')):
        raise RuntimeError('Incomplete native field coverage')
    print('production_results_sha256='+hashlib.sha256(args.results.read_bytes()).hexdigest())
    print('inner_force_and_energy_comparison=passed')


if __name__ == '__main__':
    main()
