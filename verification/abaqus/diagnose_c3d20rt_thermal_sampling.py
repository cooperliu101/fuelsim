"""Compare the identified face sampling rule with an independent native unit-flux probe."""
import csv
from pathlib import Path
import re
import numpy as np

source = Path('src/core/cartesian3d_assembly.cpp').read_text() + Path('elements/src/quad8_face.cpp').read_text()
rules = []
for name in ('corner', 'edge'):
    begin = source.index('static const std::array<AbaqusQuad8TransferSample',
                         source.index('abaqus_quad8_' + name + '_transfer_rule()'))
    end = source.index('return value;', begin)
    rule = np.array([float(value) for value in re.findall(
        r'-?\d+\.\d+(?:e[+-]?\d+)?', source[begin:end])]).reshape(-1, 3)
    rule[:, 2] /= rule[:, 2].sum()
    rules.append(rule)
weights = np.zeros(5)
for lower, upper in ((0.0, 1.2), (1.2, 2.0)):
    for constraint in range(8):
        fraction = 1.0 / 24.0 if constraint < 4 else 5.0 / 24.0
        for first, second, weight in rules[constraint // 4]:
            if constraint in (1, 2):
                first = 1.0 - first
            elif constraint == 5:
                first = 1.0 - second
            elif constraint == 7:
                first = second
            y = lower + (upper - lower) * first
            primary = min(int(y / 0.5), 3)
            coordinate = y / 0.5 - primary
            measure = (upper - lower) * fraction * weight
            weights[primary] += measure * (1.0 - coordinate)
            weights[primary + 1] += measure * coordinate
with open('verification/abaqus/c3d20rt_contact_unit_flux_nodes.csv') as stream:
    rows = list(csv.DictReader(stream))
reference = np.array([-2.0 * float(rows[node - 1]['heat_reaction']) for node in (2, 3, 21, 33, 45)])
exact = np.array([0.25, 0.5, 0.5, 0.5, 0.25])
print('identified_primary_row_weights', weights.tolist())
print('native_primary_row_weights', reference.tolist())
print('native_sampling_maximum_absolute_error', np.max(abs(weights - reference)))
print('native_vs_exact_integration_maximum_difference', np.max(abs(reference - exact)))
print('total_area', weights.sum())
