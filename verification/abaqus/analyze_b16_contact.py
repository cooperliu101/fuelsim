"""Independent output-only diagnostics of CAX4T contact pressure differences.

NETCDF_LIBRARY must identify the NetCDF C library when it is not on the loader path.
Inputs are literal production cards; this script never generates or modifies them.
"""
import csv
import gzip
from pathlib import Path
import sys

from analyze_b15_time import Exodus, metrics, np, ROOT, reference, production

OUT = ROOT / 'verification/abaqus/b16_contact'


def rows(directory, kind):
    p = directory / ('b13_small_cax4t_' + kind + '.csv')
    stream = p.open() if p.exists() else gzip.open(str(p) + '.gz', 'rt')
    with stream:
        yield from csv.DictReader(stream)


def write(name, values):
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def compare(label, native, result):
    ref = [{k: float(v) for k, v in r.items()} for r in rows(native, 'contact')]
    e = Exodus(result)
    times, names = e.get('time_whole'), e.names('name_nod_var')
    rt = np.array([r['time'] for r in ref])
    ti = np.searchsorted(times, rt)
    assert len(times) == 2561 and len(ref) == 12800
    assert np.max(abs(times[ti] - rt)) < 1e-12
    nodes = np.array([int(r['node']) - 1 for r in ref])

    def nodal(name):
        return e.get('vals_nod_var' + str(names.index(name) + 1))

    comparisons, samples = [], {}
    for field in ['pressure', 'gap', 'normal_force']:
        a = nodal('contact_' + field + '_fuel_cladding')[ti, nodes]
        b = np.array([np.hypot(r['normal_r'], r['normal_z']) if field == 'normal_force' else r[field]
                      for r in ref])
        m = metrics(a[:, None], b[:, None])
        nz = b != 0
        ids = np.flatnonzero(nz)
        worst = ids[np.argmax(abs(a[nz] - b[nz]) / abs(b[nz]))]
        comparisons.append(dict(case=label, field=field, **m, time=rt[worst], node=nodes[worst] + 1,
                                actual=a[worst], reference=b[worst]))
        samples[field] = a, b
    a, b = samples['pressure']
    ag, bg = samples['gap']
    comparisons.append(dict(case=label, field='abaqus_pressure_law',
                            **metrics(b[:, None], (1e14 * np.maximum(-bg, 0))[:, None]),
                            time=None, node=None, actual=None, reference=None))
    comparisons.append(dict(case=label, field='fuelsim_pressure_law',
                            **metrics(a[:, None], (1e14 * np.maximum(-ag, 0))[:, None]),
                            time=None, node=None, actual=None, reference=None))
    timeline = [dict(case=label, time=rt[i], node=nodes[i] + 1, actual_pressure=a[i], reference_pressure=b[i],
                     actual_gap=ag[i], reference_gap=bg[i]) for i in range(len(ref))
                if (8.9 <= rt[i] <= 9.2 and nodes[i] == 27) or (6.5 <= rt[i] <= 6.65 and nodes[i] == 34)]

    # Reconstruct the straight primary segment and secondary nodal area from
    # each program's own output at the original worst-error sample.
    focus = 9.0390625
    native_nodes = {}
    for r in rows(native, 'nodes'):
        if abs(float(r['time']) - focus) < 1e-12:
            native_nodes[int(r['node']) - 1] = [float(r['ur']), float(r['uz'])]
        if float(r['time']) > focus:
            break
    assert len(native_nodes) == 53
    coords = np.stack([e.get('coordx'), e.get('coordy')], axis=-1)
    index = int(np.searchsorted(times, focus))
    fs_u = np.stack([nodal('displacement_r')[index], nodal('displacement_z')[index]], axis=-1)
    ab_u = np.array([native_nodes[i] for i in range(53)])
    geometry = []
    for solver, u, pressure, gap, force in [
            ('fuelsim', fs_u, a, ag, samples['normal_force'][0]),
            ('abaqus', ab_u, b, bg, samples['normal_force'][1])]:
        x = coords + u
        # Original mesh: node 28 projects to primary segment 46--49;
        # secondary adjacent edges are 21--28 and 28--35.
        s, first, last = x[27], x[45], x[48]
        t = last - first
        normal = np.array([t[1], -t[0]]) / np.linalg.norm(t)
        alpha = np.dot(s - first, t) / np.dot(t, t)
        raw_gap = np.dot(first - s, normal)
        area = sum(np.pi * np.linalg.norm(s - x[j]) * (2 * s[0] + x[j, 0]) / 3 for j in [20, 34])
        i = int(np.flatnonzero((rt == focus) & (nodes == 27))[0])
        geometry.append(dict(case=label, solver=solver, time=focus, node=28, primary_alpha=alpha,
                             straight_gap=raw_gap, output_gap=gap[i], straight_pressure=max(-raw_gap, 0) * 1e14,
                             output_pressure=pressure[i], current_area=area,
                             output_force_over_pressure=force[i] / pressure[i] if pressure[i] else None,
                             radial_displacement=u[27, 0], axial_displacement=u[27, 1]))
    e.close()
    return comparisons, timeline, geometry


def integration_at_peak(native, result):
    """Identify the local time rule from independently extracted stress and CEEQ."""
    old_time, time = 9.03125, 9.0390625
    states = {}
    for row in rows(native, 'points'):
        t, element, point = float(row['time']), int(row['element']), int(row['point'])
        if t in [old_time, time] and element >= 25:
            states[t, element, point] = {k: float(v) for k, v in row.items()}
        if t > time:
            break
    e = Exodus(result)
    times, names = e.get('time_whole'), e.names('name_elem_var')
    indices = [int(np.searchsorted(times, t)) for t in [old_time, time]]
    cache = {}
    for q in [0, 1, 3, 2]:
        for key in ['stress_rr', 'stress_zz', 'stress_hoop', 'stress_rz', 'equiv_creep', 'equiv_plastic']:
            name = key + '_q' + str(q)
            cache[name] = e.get('vals_elem_var%deb2' % (names.index(name) + 1))[indices]
    e.close()

    def mises(state):
        x, y, z, w = [state['stress_' + c] for c in ['rr', 'zz', 'hoop', 'rz']]
        return np.sqrt(.5 * ((x-y)**2 + (y-z)**2 + (z-x)**2) + 3*w*w)

    output = []
    for element in range(25, 35):
        for point in range(1, 5):
            q = [0, 1, 3, 2][point - 1]
            for solver in ['abaqus', 'fuelsim']:
                if solver == 'abaqus':
                    old, current = states[old_time, element, point], states[time, element, point]
                else:
                    old, current = [{key: cache[key + '_q' + str(q)][i, element - 25]
                                     for key in ['stress_rr', 'stress_zz', 'stress_hoop', 'stress_rz',
                                                 'equiv_creep', 'equiv_plastic']} for i in [0, 1]]
                increment = current['equiv_creep'] - old['equiv_creep']
                forward = (time-old_time) * 8e-26 * mises(old)**3
                backward = (time-old_time) * 8e-26 * mises(current)**3
                assert increment > 0 and forward > 0 and backward > 0
                ef = abs(increment-forward) / max(increment, forward)
                eb = abs(increment-backward) / max(increment, backward)
                mode = ('indistinguishable' if ef < 1e-7 and eb < 1e-7 else
                        'explicit' if ef < 1e-7 else 'implicit' if eb < 1e-7 else 'other')
                output.append(dict(solver=solver, time=time, element=element, point=point, mode=mode,
                                   creep_increment=increment, old_stress_prediction=forward,
                                   new_stress_prediction=backward, forward_relative_difference=ef,
                                   backward_relative_difference=eb, equiv_plastic=current['equiv_plastic']))
    write('creep_integration_at_peak.csv', output)


def common_fields(label, native, result):
    ref, cn = reference(native)
    actual, fixed = production(result, cn)
    assert np.max(abs(ref['displacement'][fixed])) < 1e-12
    ref['displacement'][fixed] = 0.
    output = []
    for field in actual:
        if np.any(ref[field] != 0):
            value = metrics(actual[field], ref[field])
            defined = True
        else:
            difference = float(np.linalg.norm(actual[field] - ref[field], axis=1).max())
            value = dict(relative_l2_percent=None, relative_peak_percent=None, maximum_pointwise_percent=None,
                         maximum_absolute_difference=difference, zero_reference_count=len(ref[field]),
                         maximum_zero_reference_absolute_difference=difference)
            defined = False
        output.append(dict(case=label, field=field, relative_defined=defined, **value))
    return output


def main():
    native_root = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
    base = ROOT / 'build/b15_time/dt00078125/transient_b15_time_00078125_results.e'
    nocreep = ROOT / 'verification/fuelsim/transient_b16_nocreep_results.e'
    if not nocreep.exists():
        nocreep = ROOT / 'build/b16_contact/transient_b16_nocreep_results.e'
    cases = [('baseline', ROOT / 'verification/abaqus/b15_time/dt00078125', base),
             ('smooth0', native_root / 'smooth0', base),
             ('nocreep', native_root / 'nocreep', nocreep),
             ('k002', native_root / 'k002', base),
             ('nocreep_k002', native_root / 'nocreep_k002', nocreep)]
    all_metrics, timeline, geometry, fields = [], [], [], []
    for label, native, result in cases:
        m, t, g = compare(label, native, result)
        all_metrics.extend(m)
        timeline.extend(t)
        geometry.extend(g)
        fields.extend(common_fields(label, native, result))
    write('contact_metrics.csv', all_metrics)
    write('contact_activation_history.csv', timeline)
    write('same_state_geometry.csv', geometry)
    write('common_time_fields.csv', fields)
    integration_at_peak(ROOT / 'verification/abaqus/b15_time/dt00078125', base)


if __name__ == '__main__':
    main()
