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


def failed_metrics(record):
    fields = ('displacement_L2', 'axial_stress_L2', 'stress_tensor_L2', 'bending_moment_L2', 'energy_relative')
    failed = [name for name in fields if record[name] != '' and record[name] >= 0.01]
    if record['case'].startswith('axial') and record['bending_moment_absolute'] >= 1e-10:
        failed.append('zero_bending_moment_absolute')
    return failed


def lagrange(locations, values):
    weights = np.ones((len(values), len(locations)))
    for i in range(len(locations)):
        for j in range(len(locations)):
            if i != j:
                weights[:, i] *= (values - locations[j]) / (locations[i] - locations[j])
    return weights


def read_solid(path):
    with Dataset(path) as data:
        xyz = np.asarray(np.column_stack([data['coord' + c][:] for c in 'xyz']))
        node_names = names(data, 'name_nod_var')
        u = np.asarray(np.column_stack([data[f'vals_nod_var{node_names.index("displacement_" + c) + 1}'][-1]
                                        for c in 'xyz']))
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
    if not all(np.isfinite(field).all() for field in (xyz, u, stress, volume)) or not np.isfinite(energy):
        raise ValueError('Nonfinite solid reference field')
    return xyz, u, positions, stress, volume, area, float(energy), lower, upper


def read_modal(path, positions, mesh_path, target_nodes=None):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if any(None in row for row in rows):
        raise ValueError('Modal CSV contains inconsistent column counts')
    nodes = [row for row in rows if row['kind'] == 'node']
    if [int(row['id']) for row in nodes] != list(range(len(nodes))):
        raise ValueError('Modal nodal coverage is incomplete')
    u = np.array([[float(row[c]) for c in ('ux', 'uy', 'uz')] for row in nodes])
    xyz = np.array([[float(row[c]) for c in 'xyz'] for row in nodes])
    points = [row for row in rows if row['kind'] in ('point', 'solid_point')]
    sections = [row for row in rows if row['kind'] == 'section']
    points_by_z = {}
    for row in points:
        points_by_z.setdefault(float(row['z']), []).append(row)
    if not sections:
        raise ValueError('Production section-resultant output is missing')
    # Re-integrate the exported material-point stresses independently of the
    # production resultant fields. This fixture has area H*W = 4e-5 m^2.
    for section in sections:
        selected = points_by_z[float(section['z'])]
        axial_weight = sum(float(row['weight']) for row in selected) / 4e-5
        force = sum(float(row['weight']) * float(row['szz']) for row in selected) / axial_weight
        moment_x = sum(float(row['weight']) * float(row['y']) * float(row['szz']) for row in selected) / axial_weight
        moment_y = -sum(float(row['weight']) * float(row['x']) * float(row['szz']) for row in selected) / axial_weight
        for field, value in [('axial_force', force), ('moment_x', moment_x), ('moment_y', moment_y)]:
            if abs(value - float(section[field])) > 1e-10 * abs(value) + 1e-12:
                raise ValueError('Production section resultant disagrees with material-point stresses')
    samples = {tuple(round(float(row[c]), 14) for c in 'xyz'):
               [float(row[c]) for c in ('sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz', 'ux', 'uy', 'uz')] for row in points}
    if len(samples) != len(points):
        raise ValueError('Duplicate modal material-point coordinates')
    solid_samples = {tuple(round(float(row[c]), 14) for c in 'xyz')
                     for row in points if row['kind'] == 'solid_point'}
    with Dataset(mesh_path) as mesh:
        coordinates = np.asarray(np.column_stack([mesh['coord' + c][:] for c in 'xyz']))
        connectivity = np.asarray(mesh['connect1'][:], dtype=int) - 1
        lower = coordinates[connectivity[:, :8]].min(axis=1)
        upper = coordinates[connectivity[:, :8]].max(axis=1)
    g3, _ = np.polynomial.legendre.leggauss(3)
    g6, _ = np.polynomial.legendre.leggauss(6)
    target = positions.reshape(-1, 3)
    interpolated = np.empty((len(target), 6))
    coverage = np.zeros(len(target), dtype=int)
    displacement_positions = xyz if target_nodes is None else target_nodes
    interpolated_u = np.zeros((len(displacement_positions), 3))
    displacement_coverage = np.zeros(len(displacement_positions), dtype=int)
    # Restrict the coordinate search to one axial interval. This changes no
    # samples or interpolation, and avoids scanning the entire 3D reference
    # once for every transverse cell on a refined extrusion.
    axial_targets = {}
    axial_nodes = {}
    for lo, hi in zip(lower, upper):
        interval = (lo[2], hi[2])
        if interval not in axial_targets:
            axial_targets[interval] = np.flatnonzero((target[:, 2] > lo[2]) & (target[:, 2] < hi[2]))
            axial_nodes[interval] = np.flatnonzero((displacement_positions[:, 2] >= lo[2] - 1e-14)
                                                 & (displacement_positions[:, 2] <= hi[2] + 1e-14))

    for lo, hi in zip(lower, upper):
        # The comparison fixtures are affine extrusions. Their modal stress is
        # at most quadratic in each transverse coordinate and quintic in z.
        # Tensor interpolation therefore reconstructs the actual modal field,
        # including at every point of a separately refined 3D reference mesh.
        candidates = axial_targets[(lo[2], hi[2])]
        inside = np.all((target[candidates, :2] > lo[:2]) & (target[candidates, :2] < hi[:2]), axis=1)
        selected = candidates[inside]
        positions_x, positions_y = lo[:2, None] + (hi - lo)[:2, None] * (g3 + 1) / 2
        first_native = tuple(round(float(v), 14) for v in lo + (hi - lo) * (g3[0] + 1) / 2)
        axial_rule = g3 if first_native in solid_samples else g6
        positions_z = lo[2] + (hi[2] - lo[2]) * (axial_rule + 1) / 2
        values = np.array([samples[tuple(round(float(v), 14) for v in (x, y, z))]
                           for z in positions_z for y in positions_y for x in positions_x]).reshape(len(axial_rule), 3, 3, 9)
        natural = 2 * (target[selected] - lo) / (hi - lo) - 1
        nx = lagrange(g3, natural[:, 0])
        ny = lagrange(g3, natural[:, 1])
        nz = lagrange(axial_rule, natural[:, 2])
        interpolated[selected] = np.einsum('pi,pj,pk,kjic->pc', nx, ny, nz, values[:, :, :, :6], optimize=True)
        coverage[selected] += 1
        candidates = axial_nodes[(lo[2], hi[2])]
        inside = np.all((displacement_positions[candidates, :2] >= lo[:2] - 1e-14)
                        & (displacement_positions[candidates, :2] <= hi[:2] + 1e-14), axis=1)
        selected = candidates[inside]
        natural = 2 * (displacement_positions[selected] - lo) / (hi - lo) - 1
        value = np.einsum('pi,pj,pk,kjic->pc', lagrange(g3, natural[:, 0]), lagrange(g3, natural[:, 1]),
                          lagrange(axial_rule, natural[:, 2]), values[:, :, :, 6:], optimize=True)
        previous = displacement_coverage[selected] > 0
        if not np.allclose(value[previous], interpolated_u[selected[previous]], rtol=1e-9, atol=1e-13):
            raise ValueError('Reconstructed displacement is discontinuous across a modal cell trace')
        interpolated_u[selected] = value
        displacement_coverage[selected] += 1
    if not np.all(coverage == 1):
        raise ValueError('Reference point is outside the modal mesh or lies on a discontinuous cell trace')
    if not np.all(displacement_coverage > 0):
        raise ValueError('A reference displacement node is outside the modal mesh')
    # The independently exported nodal values check the polynomial displacement
    # reconstruction wherever these nested verification meshes share a node.
    node_ids = {tuple(np.round(p, 14)): i for i, p in enumerate(displacement_positions)}
    common = [(i, node_ids[tuple(np.round(p, 14))]) for i, p in enumerate(xyz)
              if tuple(np.round(p, 14)) in node_ids]
    if not common or not np.allclose(u[[a for a, _ in common]], interpolated_u[[b for _, b in common]],
                                    rtol=1e-9, atol=1e-13):
        raise ValueError('Point and nodal displacement outputs disagree')
    stress = interpolated.reshape(positions.shape[:2] + (6,))
    integration_energy = 0.0
    for row in points:
        strain = np.array([float(row[c]) for c in ('exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz')])
        sigma = np.array([float(row[c]) for c in ('sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz')])
        integration_energy += 0.5 * float(row['weight']) * np.dot(strain * [1, 1, 1, 2, 2, 2], sigma)
    return displacement_positions, interpolated_u, stress, integration_energy


def resultants(positions, stress, area):
    points = np.asarray(positions).reshape(-1, 3)
    _, section = np.unique(np.round(points[:, 2], 12), return_inverse=True)
    force = (stress[:, :, 2] * area).reshape(-1)
    return np.column_stack([np.bincount(section, weights=force),
                            np.bincount(section, weights=force * points[:, 1]),
                            np.bincount(section, weights=-force * points[:, 0])])


def stress_diagnostics(path, positions, actual, reference, volume):
    """Localize errors without removing any point from the full-domain acceptance."""
    weights = volume[:, :, None] * [1, 1, 1, 2, 2, 2]
    difference = weights * (actual - reference) ** 2
    squared = weights * reference ** 2
    global_error, global_reference = np.sum(difference), np.sum(squared)
    z = positions[:, :, 2]
    rows = []
    for family, selections in [
            ('region', [('root_0_6mm', z < 0.006), ('tip_194_200mm', z > 0.194),
                        ('transition_6_30mm', ((z >= 0.006) & (z < 0.03)) | ((z > 0.17) & (z <= 0.194))),
                        ('interior_30_170mm', (z >= 0.03) & (z <= 0.17))]),
            ('component', [(name, i) for i, name in enumerate(('xx', 'yy', 'zz', 'xy', 'yz', 'xz'))])]:
        for name, selection in selections:
            error = np.sum(difference[selection] if family == 'region' else difference[:, :, selection])
            norm = np.sum(squared[selection] if family == 'region' else squared[:, :, selection])
            rows.append(dict(family=family, name=name,
                             squared_error_share=float(error / global_error) if global_error else 0.0,
                             error_over_global_reference=float(np.sqrt(error / global_reference)),
                             relative_L2=float(np.sqrt(error / norm)) if norm else '',
                             absolute_weighted_error=float(np.sqrt(error))))
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check_point_support(xyz, displacement):
    selections = [((0, 0, 0), (0, 1, 2)), ((0.001, 0, 0), (1, 2)), ((0, 0.01, 0), (2,))]
    rigid_rows = []
    for position, components in selections:
        nodes = np.flatnonzero(np.all(np.abs(xyz - position) < 1e-14, axis=1))
        if len(nodes) != 1:
            raise ValueError('A three-point support location is missing or duplicated')
        if np.max(np.abs(displacement[nodes[0], list(components)])) > 1e-10:
            raise ValueError('A prescribed three-point support displacement is nonzero')
        x, y, z = xyz[nodes[0]]
        rotation = np.array([[0, z, -y], [-z, 0, x], [y, -x, 0]])
        rigid_rows.extend(np.column_stack([np.eye(3), rotation])[list(components)])
    if np.linalg.matrix_rank(np.array(rigid_rows)) != 6:
        raise ValueError('The six scalar constraints do not remove all rigid motions')
    root_displacement = float(np.max(np.abs(displacement[xyz[:, 2] == 0])))
    if root_displacement <= 1e-15:
        raise ValueError('This loaded point-supported fixture unexpectedly retained a fully fixed root')
    return root_displacement


def displacement_diagnostics(path, xyz, actual, reference):
    """Explain displacement error; never align fields used by the acceptance metrics."""
    difference = actual - reference
    position = xyz - np.mean(xyz, axis=0)
    translation = np.mean(difference, axis=0)
    inertia = np.eye(3) * np.sum(position * position) - position.T @ position
    rotation = np.linalg.solve(inertia, np.sum(np.cross(position, difference - translation), axis=0))
    remainder = difference - translation - np.cross(rotation, position)
    error = np.sum(difference * difference)
    row = dict(rigid_fit_squared_error_fraction=1 - np.sum(remainder * remainder) / error if error else 0.0,
               remaining_displacement_L2=relative(remainder + reference, reference),
               translation_x=translation[0], translation_y=translation[1], translation_z=translation[2],
               rotation_x=rotation[0], rotation_y=rotation[1], rotation_z=rotation[2])
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=row, lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)


def compare(directory, case, counts, mesh_name='plate.e', accuracy_count=None, reference_path=None,
            diagnostics=False, point_support=False):
    reference_path = directory / f'{case}_solid.e' if reference_path is None else reference_path
    xyz, solid_u, positions, solid_s, volume, area, solid_energy, lower, upper = read_solid(reference_path)
    solid_g = resultants(positions, solid_s, area)
    records = []
    previous_energy = 0.0
    previous_count = 0
    for count in counts:
        modal_xyz, modal_u, modal_s, energy = read_modal(
            directory / f'{case}_{count}_history.csv', positions, directory / mesh_name, xyz)
        report = summary(directory / f'{case}_{count}_summary.csv')
        if not np.allclose(modal_xyz, xyz, rtol=0, atol=1e-14):
            raise ValueError('Solid and modal source nodes do not match')
        if not np.isfinite(modal_s).all() or not np.isfinite(modal_u).all():
            raise ValueError('Nonfinite production result')
        if abs(energy / report['strain_energy'] - 1) > 1e-10:
            raise ValueError('Material-point energy does not match the production integral')
        if abs(2 * energy / report['external_work'] - 1) > 1e-8:
            raise ValueError('Surface load work and material-point energy disagree')
        if count >= 6 and previous_count >= 6 and energy < previous_energy * (1 - 1e-8):
            raise ValueError('Modal enrichment decreased compliance in this fixed-load convergence check')
        previous_energy = energy
        previous_count = count
        if report['equilibrium_relative_residual'] > 1e-7 or report['constraint_maximum_error'] > 1e-9:
            raise ValueError('Production equilibrium or physical constraints failed')
        if point_support:
            modal_root_displacement = check_point_support(modal_xyz, modal_u)
            solid_root_displacement = check_point_support(xyz, solid_u)
            if report['constraint_count'] != 6:
                raise ValueError('The modal problem did not impose exactly six independent physical constraints')
        elif np.max(np.abs(modal_u[xyz[:, 2] == 0])) > 1e-10:
            raise ValueError('Reconstructed fixed-end displacement is nonzero')
        with Dataset(directory / mesh_name) as mesh:
            corner_ids = np.asarray(mesh['connect1'][:, :8], dtype=int) - 1
            axial_coordinates = np.unique(mesh['coordz'][:][corner_ids])
            axial_nodes = len(axial_coordinates)
            source_z = np.asarray(mesh['coordz'][:])
        if report['axial_nodes'] != axial_nodes:
            raise ValueError('Axial mesh does not match the production extrusion')
        with (directory / f'{case}_{count}_history.csv').open() as stream:
            amplitudes = [row for row in csv.DictReader(stream)
                          if row['kind'] in ('amplitude', 'local_amplitude')]
        amplitude_ids = {(int(row['id']), int(row['mode'])) for row in amplitudes}
        if len(amplitude_ids) != len(amplitudes):
            raise ValueError('Duplicate axial amplitude output')
        solid_dofs = 0
        expected_nodes = set(range(axial_nodes))
        if report.get('solid_elements', 0):
            lo, hi = report['solid_lower_interface'], report['solid_upper_interface']
            expected_nodes = set(np.flatnonzero((axial_coordinates >= lo) & (axial_coordinates <= hi)))
            solid_dofs = int(3 * np.sum((source_z < lo) | (source_z > hi)))
            if report['solid_displacement_dofs'] != solid_dofs:
                raise ValueError('Independent solid displacement count disagrees with the extrusion')
        if {node for node, _ in amplitude_ids} != expected_nodes:
            raise ValueError('Axial amplitude node coverage disagrees with the mesh')
        for node in sorted(expected_nodes):
            node_rows = [row for row in amplitudes if int(row['id']) == node]
            if not node_rows:
                raise ValueError('Missing axial amplitude node')
            expected = max(int(row['mode']) for row in node_rows) + 1
            if (node_rows[0]['kind'] == 'local_amplitude'
                    and not count < expected <= report['section_basis_size']):
                raise ValueError('Local amplitude count is outside the enriched section basis')
            if node_rows[0]['kind'] == 'amplitude' and expected != count:
                raise ValueError('Retained axial node changed the full-length mode count')
            if ({int(row['mode']) for row in node_rows} != set(range(expected))
                    or len({row['kind'] for row in node_rows}) != 1):
                raise ValueError('Axial amplitude coverage does not match its section basis')
        retained = 3 * sum(row['kind'] == 'amplitude' for row in amplitudes)
        local = 3 * sum(row['kind'] == 'local_amplitude' for row in amplitudes)
        if (report['dof_count'] != retained or report['recovered_local_dofs'] != local + solid_dofs
                or report['active_axial_dofs'] != retained + local):
            raise ValueError('Global and recovered local amplitude counts disagree with production output')
        modal_g = resultants(positions, modal_s, area)
        record = dict(case=case, modes=count, modal_dofs=int(report['dof_count']),
                      local_dofs=local + solid_dofs, section_basis_size=int(report['section_basis_size']),
                      solid_displacement_dofs=solid_u.size,
                      displacement_L2=relative(modal_u, solid_u),
                      axial_stress_L2=relative(modal_s[:, :, 2], solid_s[:, :, 2], volume),
                      stress_tensor_L2=relative(modal_s, solid_s, volume[:, :, None] * [1, 1, 1, 2, 2, 2]),
                      bending_moment_L2=relative(modal_g[:, 1:], solid_g[:, 1:]) if not case.startswith('axial') else '',
                      bending_moment_absolute=float(np.max(np.abs(modal_g[:, 1:] - solid_g[:, 1:]))),
                      axial_force_absolute=float(np.max(np.abs(modal_g[:, 0] - solid_g[:, 0]))),
                      energy_relative=abs(energy / solid_energy - 1), modal_energy=energy, solid_energy=solid_energy)
        if point_support:
            record.update(modal_root_displacement=modal_root_displacement,
                          solid_root_displacement=solid_root_displacement)
        records.append(record)
        if diagnostics:
            stress_diagnostics(directory / f'{case}_{count}_diagnostics.csv', positions, modal_s, solid_s, volume)
            if point_support:
                displacement_diagnostics(directory / f'{case}_{count}_displacement_diagnostics.csv', xyz, modal_u, solid_u)
        print(', '.join(f'{key}={value:.8g}' if isinstance(value, float) else f'{key}={value}'
                        for key, value in record.items()))
    output = directory / f'{case}_comparison.csv'
    with output.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    if accuracy_count is not None:
        selected = [row for row in records if row['modes'] == accuracy_count]
        if len(selected) != 1:
            raise ValueError('The requested accuracy case was not evaluated')
        for field in ('displacement_L2', 'axial_stress_L2', 'stress_tensor_L2', 'bending_moment_L2', 'energy_relative'):
            value = selected[0][field]
            if value != '' and (not np.isfinite(value) or value >= 0.01):
                raise ValueError(f'{case}: {field}={value} does not satisfy the strict 1.0% acceptance limit')
        if case.startswith('axial') and selected[0]['bending_moment_absolute'] >= 1e-10:
            raise ValueError('Axial tension acquired a nonzero bending moment')
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('case', choices=['axial', 'bending', 'transverse', 'nonuniform'])
    args = parser.parse_args()
    compare(args.directory, args.case, [3, 4, 6, 8, 12] if args.case == 'nonuniform' else [3, 12])
