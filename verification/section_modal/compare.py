"""Read production CSV/Exodus results; never assemble or solve a mechanics problem.

The homogeneous affine HEX20 fixture has a known 3x3x3 material-point rule.
Modal stresses (degree <= 5 in z) are interpolated from their six Gauss points
to the native solid Gauss points, retaining every solid element and point.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, chartostring


def names(data, key):
    return [str(s).rstrip('\x00') for s in chartostring(data[key][:])]


def summary(path):
    with path.open() as stream:
        return {row['metric']: float(row['value']) for row in csv.DictReader(stream)}


def relative(actual, reference, weights=None):
    difference = (actual - reference) ** 2
    squared = reference ** 2
    if weights is not None:
        difference = difference * weights
        squared = squared * weights
    denominator = np.sum(squared)
    if denominator == 0:
        raise ValueError('A zero reference must use an absolute metric')
    return float(np.sqrt(np.sum(difference) / denominator))


def read_solid(path):
    with Dataset(path) as data:
        xyz = np.column_stack([data['coord' + c][:] for c in 'xyz'])
        node_names = names(data, 'name_nod_var')
        u = np.column_stack([data[f'vals_nod_var{node_names.index("displacement_" + c) + 1}'][-1]
                             for c in 'xyz'])
        conn = np.asarray(data['connect1'][:], dtype=int) - 1
        fields = names(data, 'name_elem_var')
        stress = np.empty((len(conn), 27, 6))
        for q in range(27):
            for c, component in enumerate(('xx', 'yy', 'zz', 'xy', 'yz', 'xz')):
                index = fields.index(f'stress_{component}_q{q}') + 1
                stress[:, q, c] = data[f'vals_elem_var{index}eb1'][-1]
    g, w = np.polynomial.legendre.leggauss(3)
    natural = np.array([(x, y, z) for z in g for y in g for x in g])
    weights = np.array([x * y * z for z in w for y in w for x in w])
    lower, upper = xyz[conn[:, :8]].min(axis=1), xyz[conn[:, :8]].max(axis=1)
    positions = (lower[:, None] + upper[:, None]) / 2 + natural * (upper[:, None] - lower[:, None]) / 2
    volume = np.prod((upper - lower) / 2, axis=1)[:, None] * weights
    area = volume / (np.repeat(w, 9)[None, :] * (upper - lower)[:, 2, None] / 2)
    # Isotropic compliance, tensor shear convention; no eigenstrain in these cards.
    strain = stress * 1.3 / 7e10
    strain[:, :, :3] -= 0.3 / 7e10 * np.sum(stress[:, :, :3], axis=2, keepdims=True)
    energy = 0.5 * np.sum(volume[:, :, None] * stress * strain * np.array([1, 1, 1, 2, 2, 2]))
    return xyz, u, positions, stress, volume, area, float(energy), lower, upper


def read_modal(path, positions, mesh_path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if any(None in row for row in rows):
        raise ValueError('Modal CSV contains inconsistent column counts')
    nodes = [row for row in rows if row['kind'] == 'node']
    if [int(row['id']) for row in nodes] != list(range(len(nodes))):
        raise ValueError('Modal nodal coverage is incomplete')
    u = np.array([[float(row[c]) for c in ('ux', 'uy', 'uz')] for row in nodes])
    xyz = np.array([[float(row[c]) for c in 'xyz'] for row in nodes])
    points = [row for row in rows if row['kind'] == 'point']
    sections = [row for row in rows if row['kind'] == 'section']
    if not sections:
        raise ValueError('Production section-resultant output is missing')
    # Re-integrate the exported material-point stresses independently of the
    # production resultant fields. This fixture has area H*W = 4e-5 m^2.
    for section in sections:
        selected = [row for row in points if float(row['z']) == float(section['z'])]
        axial_weight = sum(float(row['weight']) for row in selected) / 4e-5
        force = sum(float(row['weight']) * float(row['szz']) for row in selected) / axial_weight
        moment_x = sum(float(row['weight']) * float(row['y']) * float(row['szz']) for row in selected) / axial_weight
        moment_y = -sum(float(row['weight']) * float(row['x']) * float(row['szz']) for row in selected) / axial_weight
        for field, value in [('axial_force', force), ('moment_x', moment_x), ('moment_y', moment_y)]:
            if abs(value - float(section[field])) > 1e-10 * abs(value) + 1e-12:
                raise ValueError('Production section resultant disagrees with material-point stresses')
    samples = {tuple(round(float(row[c]), 14) for c in 'xyz'):
               [float(row[c]) for c in ('sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz')] for row in points}
    if len(samples) != len(points):
        raise ValueError('Duplicate modal material-point coordinates')
    with Dataset(mesh_path) as mesh:
        coordinates = np.column_stack([mesh['coord' + c][:] for c in 'xyz'])
        connectivity = np.asarray(mesh['connect1'][:], dtype=int) - 1
        lower = coordinates[connectivity[:, :8]].min(axis=1)
        upper = coordinates[connectivity[:, :8]].max(axis=1)
    g3, _ = np.polynomial.legendre.leggauss(3)
    g6, _ = np.polynomial.legendre.leggauss(6)
    target = positions.reshape(-1, 3)
    interpolated = np.empty((len(target), 6))
    coverage = np.zeros(len(target), dtype=int)

    def lagrange(locations, values):
        weights = np.ones((len(values), len(locations)))
        for i in range(len(locations)):
            for j in range(len(locations)):
                if i != j:
                    weights[:, i] *= (values - locations[j]) / (locations[i] - locations[j])
        return weights

    for lo, hi in zip(lower, upper):
        # The comparison fixtures are affine extrusions. Their modal stress is
        # at most quadratic in each transverse coordinate and quintic in z.
        # Tensor interpolation therefore reconstructs the actual modal field,
        # including at every point of a separately refined 3D reference mesh.
        selected = np.all((target > lo) & (target < hi), axis=1)
        if not np.any(selected):
            continue
        positions_x, positions_y = lo[:2, None] + (hi - lo)[:2, None] * (g3 + 1) / 2
        positions_z = lo[2] + (hi[2] - lo[2]) * (g6 + 1) / 2
        values = np.array([samples[tuple(round(float(v), 14) for v in (x, y, z))]
                           for z in positions_z for y in positions_y for x in positions_x]).reshape(6, 3, 3, 6)
        natural = 2 * (target[selected] - lo) / (hi - lo) - 1
        nx = lagrange(g3, natural[:, 0])
        ny = lagrange(g3, natural[:, 1])
        nz = lagrange(g6, natural[:, 2])
        interpolated[selected] = np.einsum('pi,pj,pk,kjic->pc', nx, ny, nz, values, optimize=True)
        coverage[selected] += 1
    if not np.all(coverage == 1):
        raise ValueError('Reference point is outside the modal mesh or lies on a discontinuous cell trace')
    stress = interpolated.reshape(positions.shape[:2] + (6,))
    integration_energy = 0.0
    for row in points:
        strain = np.array([float(row[c]) for c in ('exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz')])
        sigma = np.array([float(row[c]) for c in ('sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz')])
        integration_energy += 0.5 * float(row['weight']) * np.dot(strain * [1, 1, 1, 2, 2, 2], sigma)
    return xyz, u, stress, integration_energy


def resultants(positions, stress, area):
    sections = {}
    for e in range(len(positions)):
        for q, p in enumerate(positions[e]):
            value = stress[e, q, 2] * area[e, q]
            sections.setdefault(round(float(p[2]), 12), np.zeros(3))[:] += value * np.array([1, p[1], -p[0]])
    return np.array([sections[z] for z in sorted(sections)])


def compare(directory, case, counts):
    xyz, solid_u, positions, solid_s, volume, area, solid_energy, lower, upper = read_solid(directory / f'{case}_solid.e')
    solid_g = resultants(positions, solid_s, area)
    records = []
    previous_energy = 0.0
    for count in counts:
        modal_xyz, modal_u, modal_s, energy = read_modal(directory / f'{case}_{count}_history.csv', positions, directory / 'plate.e')
        report = summary(directory / f'{case}_{count}_summary.csv')
        if not np.allclose(modal_xyz, xyz, rtol=0, atol=1e-14):
            raise ValueError('Solid and modal source nodes do not match')
        if not np.isfinite(modal_s).all() or not np.isfinite(modal_u).all():
            raise ValueError('Nonfinite production result')
        if abs(energy / report['strain_energy'] - 1) > 1e-10:
            raise ValueError('Material-point energy does not match the production integral')
        if abs(2 * energy / report['external_work'] - 1) > 1e-8:
            raise ValueError('Surface load work and material-point energy disagree')
        if energy < previous_energy * (1 - 1e-8):
            raise ValueError('Nested modal enrichment decreased compliance under fixed loads')
        previous_energy = energy
        if report['equilibrium_relative_residual'] > 1e-7 or report['constraint_maximum_error'] > 1e-9:
            raise ValueError('Production equilibrium or physical constraints failed')
        if np.max(np.abs(modal_u[xyz[:, 2] == 0])) > 1e-10:
            raise ValueError('Reconstructed fixed-end displacement is nonzero')
        if report['dof_count'] != 27 * count:
            raise ValueError('Cross-section mesh nodes entered the global unknown count')
        modal_g = resultants(positions, modal_s, area)
        record = dict(case=case, modes=count, modal_dofs=int(report['dof_count']), solid_displacement_dofs=solid_u.size,
                      displacement_L2=relative(modal_u, solid_u),
                      axial_stress_L2=relative(modal_s[:, :, 2], solid_s[:, :, 2], volume),
                      stress_tensor_L2=relative(modal_s, solid_s, volume[:, :, None] * [1, 1, 1, 2, 2, 2]),
                      bending_moment_L2=relative(modal_g[:, 1:], solid_g[:, 1:]) if case != 'axial' else '',
                      bending_moment_absolute=float(np.max(np.abs(modal_g[:, 1:] - solid_g[:, 1:]))),
                      axial_force_absolute=float(np.max(np.abs(modal_g[:, 0] - solid_g[:, 0]))),
                      energy_relative=abs(energy / solid_energy - 1), modal_energy=energy, solid_energy=solid_energy)
        records.append(record)
        print(', '.join(f'{key}={value:.8g}' if isinstance(value, float) else f'{key}={value}'
                        for key, value in record.items()))
    # These tests distinguish implementation contracts from model-accuracy evidence.
    # No fitted error gate is used to turn a nonconverged truncation into qualification.
    output = directory / f'{case}_comparison.csv'
    with output.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('case', choices=['axial', 'bending', 'transverse', 'nonuniform'])
    args = parser.parse_args()
    compare(args.directory, args.case, [3, 4, 6, 8, 12] if args.case == 'nonuniform' else [3, 12])
