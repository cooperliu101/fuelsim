"""Replay the native state without solving or changing either input card.

The HEX20 geometry and bulk gradients are reconstructed independently in Python.
The neighborhood sampling values come from the existing identified QUAD8 rule;
this compares that rule at native temperatures, not an independent derivation of
its sampling constants. It does not assert arbitrary nonmatching-face agreement.
"""
import csv
import itertools
import json
from pathlib import Path
import sys

import netCDF4
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'abaqus'))
from diagnose_c3d20rt_finite_thermal_transfer import body_shapes
from diagnose_c3d20rt_contact_gap import FACES, shape
from check_h20_disk_transfer import samples


def main():
    with netCDF4.Dataset(ROOT / 'dc3d20_contact.e') as mesh:
        coordinates = np.array([mesh['coord' + axis][:] for axis in 'xyz']).T
        elements = np.vstack([mesh['connect1'][:], mesh['connect2'][:]]) - 1
    secondary = elements[0, FACES[1]]
    primary = elements[1, FACES[3]]
    primary = np.array([primary[np.argmin(np.linalg.norm(coordinates[primary] - coordinates[node], axis=1))]
                        for node in secondary])
    coordinate_error = float(np.max(abs(coordinates[secondary] - coordinates[primary])))
    if coordinate_error > 1e-14 or len(set(primary)) != 8:
        raise AssertionError('Diagnostic requires uniquely matching interface coordinates')
    with (ROOT / 'dc3d20_contact.nodes.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    temperature = np.array([float(row['temperature']) for row in rows])
    reaction = np.array([float(row['reaction']) for row in rows])
    body = np.zeros(len(temperature))
    points, weights = np.polynomial.legendre.leggauss(3)
    for nodes in elements:
        for indices in itertools.product(range(3), repeat=3):
            _, _, derivative = body_shapes(points[list(indices)])
            mapping = coordinates[nodes].T @ derivative
            gradient = derivative @ np.linalg.inv(mapping)
            measure = np.prod(weights[list(indices)]) * np.linalg.det(mapping)
            body[nodes] += 10 * measure * gradient @ (gradient.T @ temperature[nodes])
    outside = [node for node in range(len(temperature)) if node not in secondary and node not in primary]
    report = {'case': 'dc3d20_contact', 'matching_coordinate_maximum_m': coordinate_error,
              'noninterface_native_balance_maximum_W': float(np.max(abs(body[outside] - reaction[outside]))),
              'interface_native_balance_maximum_W': {}}
    for order in (2, 3, 5, 10):
        points, weights = np.polynomial.legendre.leggauss(order)
        matrix = np.zeros((8, 8))
        for i, j in itertools.product(range(order), repeat=2):
            values, dx, dy = shape(points[i], points[j])
            measure = np.linalg.norm(np.cross(dx @ coordinates[secondary], dy @ coordinates[secondary]))
            matrix += weights[i] * weights[j] * measure * np.outer(values, values)
        flux = 1000 * matrix @ (temperature[secondary] - temperature[primary])
        report['interface_native_balance_maximum_W'][f'pointwise_gauss_{order}'] = float(max(
            np.max(abs(body[secondary] + flux)), np.max(abs(body[primary] - flux))))
    matrix = np.zeros((8, 8))
    for constraint in range(8):
        area = 0.0
        moment = np.zeros(8)
        fraction = 1 / 24 if constraint < 4 else 5 / 24
        for a, b, weight in samples(constraint):
            values, dx, dy = shape(2 * a - 1, 2 * b - 1)
            measure = 4 * fraction * weight * np.linalg.norm(
                np.cross(dx @ coordinates[secondary], dy @ coordinates[secondary]))
            area += measure
            moment += measure * values
        matrix += np.outer(moment, moment) / area
    flux = 1000 * matrix @ (temperature[secondary] - temperature[primary])
    report['interface_native_balance_maximum_W']['neighborhood_average'] = float(max(
        np.max(abs(body[secondary] + flux)), np.max(abs(body[primary] - flux))))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
