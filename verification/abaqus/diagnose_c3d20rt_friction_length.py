"""Identify friction length from native fields; run from the repository root.

Check the existing force-transfer matrix against native normal forces before
using it to recover tangential tractions. No solver or input generation is used.
"""
import csv
from pathlib import Path
import re
import numpy as np

source = Path('src/core/cartesian3d_assembly.cpp').read_text()
begin = source.index('static const Matrix8 value', source.index('const Matrix8& abaqus_quad8_averaging'))
end = source.index('}};', begin)
averaging = np.array([
    float(value) for value in re.findall(r'-?\d+\.\d+(?:e[+-]?\d+)?', source[begin:end])
]).reshape(8, 8)
with open('verification/abaqus/c3d20rt_contact_friction_contact.csv') as stream:
    rows = list(csv.DictReader(stream))
labels = sorted({int(row['node']) for row in rows})
transfer = np.zeros((13, 13))
faces = [([57, 60, 64, 61, 68, 72, 76, 69], 1.2),
         ([60, 78, 80, 64, 83, 85, 88, 72], 0.8)]
for face, area in faces:
    for i, constraint in enumerate(face):
        fraction = 1.0 / 24.0 if i < 4 else 5.0 / 24.0
        for j, node in enumerate(face):
            transfer[labels.index(node), labels.index(constraint)] += area * fraction * averaging[i, j]

force = np.array([[float(row['t' + axis]) for axis in 'xyz'] for row in rows]).reshape(70, 13, 3)
normal = np.array([[float(row['n' + axis]) for axis in 'xyz'] for row in rows]).reshape(70, 13, 3)
slip = np.array([[float(row['slip' + axis]) for axis in 'xyz'] for row in rows]).reshape(70, 13, 3)
pressure = np.array([-1e11 * float(row['gap']) for row in rows]).reshape(70, 13)
normal_prediction = np.einsum('ij,tj->ti', transfer, pressure)
normal_error = np.max(abs(normal_prediction - normal[:, :, 0]))
print('native_normal_transfer_maximum_error', normal_error)
if normal_error >= 1e-6:
    raise RuntimeError('Native normal forces do not validate the transfer matrix')
reference = np.stack([np.linalg.solve(transfer, -value) for value in force])
identified = 0.3 * pressure[0] * np.linalg.norm(slip[0], axis=1) / np.linalg.norm(reference[0], axis=1)
length = 1e-5 * sum(np.sqrt(area) for _, area in faces) / len(faces)
print('geometric_elastic_slip_length', length)
print('identified_elastic_slip_length_minimum', identified.min())
print('identified_elastic_slip_length_maximum', identified.max())
old = np.zeros((13, 3))
last = np.zeros((13, 3))
predictions = []
for step in range(70):
    stiffness = 0.3 * pressure[step, :, None] / length
    trial = stiffness * (old + slip[step] - last)
    norm = np.linalg.norm(trial, axis=1)
    limit = 0.3 * pressure[step]
    scale = np.ones(13)
    sliding = norm > limit
    scale[sliding] = limit[sliding] / norm[sliding]
    traction = trial * scale[:, None]
    old = traction / stiffness
    last = slip[step]
    predictions.append(traction)
difference = np.array(predictions) - reference
print('native_traction_replay_relative_l2', np.linalg.norm(difference) / np.linalg.norm(reference))
print('native_traction_replay_maximum_absolute', np.max(abs(difference)))
