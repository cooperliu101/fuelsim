"""Compare direct primary-gap perturbation with reaction-derived transfer."""
import csv
from pathlib import Path
import numpy as np
from analyze_h20_28_finite_transfer import extract

root = Path(__file__).resolve().parent
reference = extract('h20_28_finite_transfer_32_probe')
columns = np.flatnonzero(np.all(abs(reference[3] - [.125, .140625]) < 1e-12, axis=1))
assert len(columns) == 1
with (root / 'h20_28_primary_signed_gap_32_probe_gap.csv').open() as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 24
labels = [int(row['node']) for row in rows if row['step'] == 'BASE']
assert len(labels) == 8 and labels == sorted(set(labels))
values = {}
for step in ['PRIMARY_PLUS', 'PRIMARY_MINUS']:
    selected = [row for row in rows if row['step'] == step]
    assert [int(row['node']) for row in selected] == labels
    values[step] = np.array([float(row['gap']) for row in selected])
direct = -(values['PRIMARY_PLUS'] - values['PRIMARY_MINUS']) / 2e-7
force = reference[2][:, columns[0]]
print('diagnostic_primary_node', 889)
print('direct_gap_transfer', direct.tolist())
print('reaction_transfer', force.tolist())
print('maximum_absolute_difference', np.max(abs(direct - force)))
assert direct[0] < -0.04 and np.max(abs(direct - force)) < 1e-6
