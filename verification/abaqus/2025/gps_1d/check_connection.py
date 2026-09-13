"""Compare native prescribed-motion cases differing only in middle-node connectivity."""
import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent


def rows(stem):
    with (root / (stem + '_points.csv')).open(newline='') as stream:
        result = list(csv.DictReader(stream))
    if len(result) != 160:
        raise RuntimeError('Expected four material points per element in all ten increments')
    keys = [(float(r['time']), int(r['element']), int(r['point'])) for r in result]
    if len(set(keys)) != 160 or len({k[0] for k in keys}) != 10:
        raise RuntimeError('Missing or duplicate native samples')
    return dict(zip(keys, result))


independent = rows('gps_two_slice_contact')
connected = rows('gps_two_slice_connected')
if independent.keys() != connected.keys():
    raise RuntimeError('Native frames or material points differ')
for prefix, absolute_tolerance in [('stress_', 1e-4), ('elastic_', 1e-13)]:
    differences, references = [], []
    pointwise = []
    for key, before in independent.items():
        delta, reference = [], []
        for component in ['rr', 'zz', 'hoop', 'rz']:
            weight = math.sqrt(2) if component == 'rz' else 1
            a = float(before[prefix + component]) * weight
            b = float(connected[key][prefix + component]) * weight
            if not math.isfinite(a) or not math.isfinite(b):
                raise RuntimeError('Nonfinite native tensor')
            delta.append(b - a)
            reference.append(a)
        dn = math.sqrt(sum(v*v for v in delta))
        rn = math.sqrt(sum(v*v for v in reference))
        differences.extend(delta)
        references.extend(reference)
        if rn:
            pointwise.append(dn/rn)
        elif dn > absolute_tolerance:
            raise RuntimeError('Zero-reference tensor changed')
    l2 = math.sqrt(sum(v*v for v in differences)/sum(v*v for v in references))
    maximum = max(abs(v) for v in differences)
    print(prefix + 'connection_change_relative_l2=%.17g' % l2)
    print(prefix + 'connection_change_maximum_pointwise_relative=%.17g' % max(pointwise))
    print(prefix + 'connection_change_maximum_absolute_component=%.17g' % maximum)
    if l2 >= 0.001 or max(pointwise) >= 0.001 or maximum > absolute_tolerance:
        raise RuntimeError('Changing middle-node connectivity changed prescribed-motion body tensors')
print('connection_only_native_comparison=passed')
