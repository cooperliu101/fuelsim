"""Reconstruct the two prescribed-motion stress fields without solving either model.

Requires netCDF4 for Exodus I/O. All constitutive calculations use Python floats.
Counterfactual fields below are algebraic diagnostics, not new solver results.
"""
import argparse
import csv
import hashlib
import math
from pathlib import Path

import netCDF4

COMPONENTS = ('rr', 'zz', 'hoop', 'rz')
E, NU, ALPHA, T0 = 1e9, 0.3, 1e-5, 600.0
LAME = E * NU / ((1 + NU) * (1 - 2 * NU))
SHEAR = E / (2 * (1 + NU))
ROOT = Path(__file__).resolve().parent


def read_csv(path, keys):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    result = {}
    for row in rows:
        key = tuple(float(row[k]) for k in keys)
        if key in result:
            raise RuntimeError('Duplicate reference key: ' + str(key))
        result[key] = {k: float(v) for k, v in row.items()}
        if not all(math.isfinite(v) for v in result[key].values()):
            raise RuntimeError('Nonfinite native reference value')
    return result


def stress(elastic):
    pressure = LAME * sum(elastic[:3])
    return [2*SHEAR*v + (pressure if i < 3 else 0) for i, v in enumerate(elastic)]


def elastic(hoop, temperature):
    thermal = ALPHA*(temperature-T0)
    return [-thermal, -thermal, hoop-thermal, 0.0]


def norm(values):
    return math.sqrt(sum(v*v*(2 if i == 3 else 1) for i, v in enumerate(values)))


def difference(a, b):
    return [x-y for x, y in zip(a, b)]


def metrics(pairs):
    reference = [norm(b) for a, b in pairs]
    delta = [norm(difference(a, b)) for a, b in pairs]
    if not all(math.isfinite(x) and x > 0 for x in reference):
        raise RuntimeError('This prescribed hot case requires nonzero tensor references')
    return (math.sqrt(sum(d*d for d in delta)/sum(r*r for r in reference)),
            max(delta)/max(reference),
            max(d/r for d, r in zip(delta, reference)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, help='Unmodified production contact-small Exodus result')
    args = parser.parse_args()
    nodes = read_csv(ROOT/'gps_two_slice_contact_nodes.csv', ('time', 'node'))
    points = read_csv(ROOT/'gps_two_slice_contact_points.csv', ('time', 'element', 'point'))
    if len(nodes) != 160 or len(points) != 160:
        raise RuntimeError('Expected all ten frames and sixteen native nodes/material points per frame')
    groups = {name: [] for name in ('stress', 'elastic', 'inner_stress', 'outer_stress',
                                   'diagnostic_old_point_temperature', 'both_averages')}
    peak = None
    error = dict.fromkeys(('gps_stress', 'abaqus_stress', 'gps_elastic', 'abaqus_elastic',
                           'difference_decomposition', 'point_temperature'), 0.0)
    with netCDF4.Dataset(args.results) as data:
        times = data['time_whole'][:].tolist()
        if times != list(range(11)):
            raise RuntimeError('Expected initial frame and all ten accepted increments')
        names = netCDF4.chartostring(data['name_elem_var'][:]).tolist()

        def field(name, step, element):
            # Native elements 1/3 are inner, 2/4 outer; Exodus stores by named block.
            block, local = (element-1) % 2 + 1, (element-1)//2
            value = float(data['vals_elem_var%deb%d' % (names.index(name)+1, block)][step, local])
            if not math.isfinite(value):
                raise RuntimeError('Nonfinite production field: '+name)
            return value

        for key, row in points.items():
            time, element, point = map(int, key)
            if key != (float(time), float(element), float(point)):
                raise RuntimeError('Unexpected native indices or time')
            if element not in (1, 2, 3, 4) or point not in (1, 2, 3, 4):
                raise RuntimeError('Unexpected native material-point mapping')
            n = [nodes[(time, 4*(element-1)+i)] for i in range(1, 5)]
            r0, r1 = n[0]['r'], n[1]['r']
            u = sum(v['ur'] for v in n)/4
            # Every mechanical displacement is uniform within this native element.
            if max(abs(v['ur']-u) for v in n) > 1e-13 or max(abs(v['uz']-n[0]['uz']) for v in n) > 1e-13:
                raise RuntimeError('This diagnostic is restricted to the prescribed-translation case')
            q = (point-1) % 2
            eta = 0.5*(1 + (-1 if q == 0 else 1)/math.sqrt(3))
            radius = (1-eta)*r0 + eta*r1
            if abs(row['r']-radius) > 1e-14:
                raise RuntimeError('Radial Gauss correspondence changed')
            temperature = (1-eta)*n[0]['temperature'] + eta*n[1]['temperature']
            mean_temperature = sum(v['temperature'] for v in n)/4
            suffix = '_q%d' % q
            if abs(field('reference_r'+suffix, time, element)-radius) > 1e-14:
                raise RuntimeError('Production/native Gauss correspondence changed')
            error['point_temperature'] = max(error['point_temperature'],
                abs(field('temperature'+suffix, time, element)-temperature))
            # In this rectangle the r-weighted average of u/r equals u/r_mid.
            hoop, mean_hoop = u/radius, u/((r0+r1)/2)
            gps_elastic = elastic(hoop, mean_temperature)
            abaqus_elastic = elastic(mean_hoop, mean_temperature)
            gps_stress, abaqus_stress = stress(gps_elastic), stress(abaqus_elastic)
            actual_stress = [field('stress_'+c+suffix, time, element) for c in COMPONENTS]
            actual_elastic = [field('elastic_'+c+suffix, time, element) for c in COMPONENTS]
            native_stress = [row['stress_'+c] for c in COMPONENTS]
            native_elastic = [row['elastic_'+c] for c in COMPONENTS]
            for name, a, b in [('gps_stress', actual_stress, gps_stress),
                               ('abaqus_stress', native_stress, abaqus_stress),
                               ('gps_elastic', actual_elastic, gps_elastic),
                               ('abaqus_elastic', native_elastic, abaqus_elastic)]:
                error[name] = max(error[name], max(abs(v) for v in difference(a, b)))
            mechanical = stress([0.0, 0.0, hoop-mean_hoop, 0.0])
            thermal = [0.0]*4
            observed = difference(actual_stress, native_stress)
            error['difference_decomposition'] = max(error['difference_decomposition'],
                max(abs(observed[i]-mechanical[i]-thermal[i]) for i in range(4)))
            groups['stress'].append((actual_stress, native_stress))
            groups['elastic'].append((actual_elastic, native_elastic))
            groups['inner_stress' if element % 2 else 'outer_stress'].append((actual_stress, native_stress))
            # This reconstructs the former rule only as an algebraic diagnostic.
            groups['diagnostic_old_point_temperature'].append((stress(elastic(hoop, temperature)), native_stress))
            groups['both_averages'].append((abaqus_stress, native_stress))
            relative = norm(observed)/norm(native_stress)
            if peak is None or relative > peak[0]:
                peak = (relative, key, temperature, mean_temperature, actual_stress, native_stress,
                        actual_elastic, native_elastic, mechanical, thermal)

    print('scope=postprocessing_of_unmodified_production_and_native_outputs')
    print('sample_count=160')
    print('nonzero_tensor_metrics=relative_l2,relative_absolute_peak,maximum_pointwise_relative')
    print('metric_units=fraction; multiply by 100 for percent')
    for name, values in groups.items():
        print(name+'='+','.join('%.17g' % v for v in metrics(values)))
    for name, value in error.items():
        print(name+'_maximum_absolute_component_error=%.17g' % value)
    print('peak_relative=%.17g' % peak[0])
    print('peak_time_element_point='+','.join(str(int(v)) for v in peak[1]))
    print('peak_gps_expansion_temperature=%.17g' % peak[3])
    print('peak_gps_material_point_temperature=%.17g' % peak[2])
    print('peak_abaqus_expansion_temperature=%.17g' % peak[3])
    for name, values in zip(('gps_stress', 'abaqus_stress', 'gps_elastic', 'abaqus_elastic',
                             'mechanical_stress_difference', 'thermal_stress_difference'), peak[4:]):
        print('peak_'+name+'='+','.join('%.17g' % v for v in values))
    print('former_point_temperature_rule_outer_exact_relative=%.17g' % (1/math.sqrt(3)))
    print('production_results_sha256='+hashlib.sha256(args.results.read_bytes()).hexdigest())
    for name in ('gps_two_slice_contact_nodes.csv', 'gps_two_slice_contact_points.csv'):
        print(name+'_sha256='+hashlib.sha256((ROOT/name).read_bytes()).hexdigest())
    for name, value in error.items():
        gate = 1e-9 if name == 'point_temperature' else (1e-13 if 'elastic' in name else 1e-4)
        if value > gate:
            raise RuntimeError('Independent reconstruction failed: '+name)
    print('stress_difference_reconstruction=passed')


if __name__ == '__main__':
    main()
