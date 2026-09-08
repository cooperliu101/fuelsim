"""Recover the contact matrix from prescribed-temperature native output, without solving."""
import csv
import itertools
import sys
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent
stem = sys.argv[1] if len(sys.argv) > 1 else 'c3d20rt_contact_thermal_matrix'
if stem not in ('c3d20rt_contact_thermal_matrix', 'c3d20rt_contact_matching_matrix', 'c3d20rt_finite_matching_matrix'):
    raise RuntimeError('Unknown matrix probe')
matching = stem in ('c3d20rt_contact_matching_matrix', 'c3d20rt_finite_matching_matrix')
nodes = {}
elements = []
section = ''
for line in (root / (stem + '.inp')).read_text().splitlines():
    if line.startswith('*'):
        section = line.split(',')[0].lower()
    elif section == '*node':
        values = line.split(',')
        nodes[int(values[0])] = list(map(float, values[1:4]))
    elif section == '*element':
        elements.append(list(map(int, line.split(',')))[1:9])
labels = [2, 3, 6, 7, 57, 60, 61, 64] if matching else [2, 3, 6, 7, 21, 23, 33, 35, 45, 47, 57, 60, 61, 64, 78, 80]
signs = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                  [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
operators = []
for element in elements:
    coordinates = np.array([nodes[node] for node in element])
    stiffness = np.zeros((8, 8))
    capacity = np.zeros((8, 8))
    for point in itertools.product([-1/np.sqrt(3), 1/np.sqrt(3)], repeat=3):
        factors = 1 + signs * np.array(point)
        shape = np.prod(factors, axis=1) / 8
        derivative = np.array([signs[:,axis] * np.prod(np.delete(factors, axis, axis=1), axis=1) / 8
                               for axis in range(3)]).T
        jacobian = coordinates.T @ derivative
        gradient = derivative @ np.linalg.inv(jacobian)
        weight = np.linalg.det(jacobian)
        stiffness += 10 * weight * gradient @ gradient.T
        capacity += weight * np.outer(shape, shape)
    operators.append((np.array(element)-1, stiffness, capacity))
with (root / 'c3d20rt_contact_unit_flux_nodes.csv').open() as stream:
    baseline = np.array([float(row['heat_reaction']) for row in csv.DictReader(stream)])
with (root / (stem + '_nodes.csv')).open() as stream:
    rows = list(csv.DictReader(stream))
if matching:
    baseline = np.array([float(row['heat_reaction']) for row in rows[:88]])
    rows = rows[88:]
if len(rows) != len(labels) * 88:
    raise RuntimeError('Expected sixteen complete prescribed-temperature increments')
old = np.array([300.0 if node < 57 else 301.0 for node in range(1, 89)])
columns = []
for column in range(len(labels)):
    frame = rows[column*88:(column+1)*88]
    temperature = np.array([float(row['temperature']) for row in frame])
    reaction = np.array([float(row['heat_reaction']) for row in frame])
    residual = np.zeros(88)
    for indices, stiffness, capacity in operators:
        residual[indices] += stiffness @ temperature[indices] + capacity @ (temperature[indices]-old[indices]) / 0.1
    columns.append((reaction-residual-baseline)[np.array(labels)-1])
    old = temperature
matrix = np.array(columns).T
with (root / (stem + '.csv')).open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['node'] + labels)
    for label, row in zip(labels, matrix):
        writer.writerow([label] + list(row))
print('maximum_asymmetry', np.max(abs(matrix-matrix.T)))
print('maximum_constant_temperature_residual', np.max(abs(matrix.sum(axis=1))))
print('eigenvalues', np.linalg.eigvalsh(matrix).tolist())
