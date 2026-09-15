"""Automatic and manual verification driver: copy complete tracked cards unchanged and run the production CLI."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from compare import compare, summary


def stage(source, destination, name, mesh_name):
    destination.mkdir(parents=True, exist_ok=True)
    for filename in (mesh_name, name + '.fsi'):
        shutil.copyfile(source / filename, destination / filename)
        if (source / filename).read_bytes() != (destination / filename).read_bytes():
            raise ValueError('Staged production input differs from its tracked source')
    for suffix in ('_history.csv', '_summary.csv', '.e'):
        (destination / (name + suffix)).unlink(missing_ok=True)


def run(executable, source, destination, name, launcher, mesh_name='plate.e'):
    stage(source, destination, name, mesh_name)
    command = launcher + [str(executable), '-i', str(destination / (name + '.fsi'))]
    log = destination / (name + '.log')
    with log.open('w') as stream:
        process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                 env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    if process.returncode:
        raise RuntimeError(f'{name} failed with {process.returncode}\n{log.read_text()}')


def mpi_check(serial, parallel, name='nonuniform_12'):
    left, right = summary(serial / (name + '_summary.csv')), summary(parallel / (name + '_summary.csv'))
    for field in ('mode_count', 'axial_nodes', 'dof_count', 'constraint_count', 'section_basis_size',
                  'recovered_local_dofs', 'active_axial_dofs', 'global_system_size', 'local_system_size',
                  'independent_axial_warping_modes', 'shear_free_distortion_modes', 'transverse_corrector_modes'):
        if left[field] != right[field]:
            raise ValueError(f'MPI changed {field}')
    for field in ('strain_energy', 'external_work'):
        if not np.isfinite([left[field], right[field]]).all() or abs(right[field] / left[field] - 1) > 1e-8:
            raise ValueError(f'MPI changed {field}')
    with (serial / (name + '_history.csv')).open() as stream:
        a = list(csv.DictReader(stream))
    with (parallel / (name + '_history.csv')).open() as stream:
        b = list(csv.DictReader(stream))
    if len(a) != len(b):
        raise ValueError('MPI changed result coverage')
    for left_row, right_row in zip(a, b):
        if any(left_row[field] != right_row[field] for field in ('kind', 'id', 'x', 'y', 'z')):
            raise ValueError('MPI changed physical result locations or ordering')
    fields_by_kind = [('node', ['ux', 'uy', 'uz']), ('point', ['sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz'])]
    if any(row['kind'] == 'solid_point' for row in a):
        fields_by_kind.append(('solid_point', ['sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz']))
    for kind, fields in fields_by_kind:
        x = np.array([[float(row[field]) for field in fields] for row in a if row['kind'] == kind])
        y = np.array([[float(row[field]) for field in fields] for row in b if row['kind'] == kind])
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError('MPI result contains a nonfinite field')
        error = np.linalg.norm(x - y) / np.linalg.norm(x)
        print(f'MPI {kind} relative L2 = {error:.12g}')
        if error > 1e-8:
            raise ValueError('MPI changed physical fields')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('case', choices=['axial', 'bending', 'transverse', 'nonuniform', 'prescribed'])
    parser.add_argument('--mpiexec', type=Path)
    parser.add_argument('--serial', type=Path)
    parser.add_argument('--accuracy', action='store_true')
    parser.add_argument('--solid-ends', action='store_true')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    if args.solid_ends:
        name = args.case + '_hybrid_12'
        if args.mpiexec:
            if args.serial is None:
                args.serial = args.directory / 'serial'
                run(args.executable, source, args.serial, name, [], 'plate_solid_ends_128_32.e')
            run(args.executable, source, args.directory, name, [str(args.mpiexec), '-n', '2'], 'plate_solid_ends_128_32.e')
            mpi_check(args.serial, args.directory, name)
        else:
            reference = args.case + '_support_solid'
            run(args.executable, source, args.directory, reference, [], 'plate_support_reference.e')
            run(args.executable, source, args.directory, name, [], 'plate_solid_ends_128_32.e')
            compare(args.directory, args.case + '_hybrid', [12], 'plate_solid_ends_128_32.e', accuracy_count=12,
                    reference_path=args.directory / (reference + '.e'), diagnostics=True, point_support=True)
    elif args.accuracy and args.mpiexec:
        if args.case != 'nonuniform':
            raise ValueError('Refined MPI comparison requires the nonuniform case')
        name = 'nonuniform_local_12'
        if args.serial is None:
            args.serial = args.directory / 'serial'
            run(args.executable, source, args.serial, name, [], 'plate_local.e')
        run(args.executable, source, args.directory, name, [str(args.mpiexec), '-n', '2'], 'plate_local.e')
        mpi_check(args.serial, args.directory, name)
    elif args.accuracy:
        if args.case == 'prescribed':
            raise ValueError('Accuracy comparison requires one of the four mechanical loads')
        case = args.case + '_local'
        reference = args.case + '_refined_solid'
        run(args.executable, source, args.directory, reference, [], 'plate_reference.e')
        run(args.executable, source, args.directory, case + '_12', [], 'plate_local.e')
        compare(args.directory, case, [12], 'plate_local.e', accuracy_count=12,
                reference_path=args.directory / (reference + '.e'))
    elif args.mpiexec:
        run(args.executable, source, args.directory, 'nonuniform_12', [str(args.mpiexec), '-n', '2'])
        mpi_check(args.serial, args.directory)
    elif args.case == 'prescribed':
        run(args.executable, source, args.directory, 'prescribed_8', [])
        report = summary(args.directory / 'prescribed_8_summary.csv')
        if report['strain_energy'] <= 0 or abs(2 * report['strain_energy'] / report['external_work'] - 1) > 1e-8:
            raise ValueError('Nonzero prescribed displacement reaction work failed')
        with (args.directory / 'prescribed_8_history.csv').open() as stream:
            for row in csv.DictReader(stream):
                if row['kind'] != 'node':
                    continue
                if float(row['z']) == 0 and max(abs(float(row[c])) for c in ('ux', 'uy', 'uz')) > 1e-12:
                    raise ValueError('Fixed-end physical displacement failed')
                if float(row['z']) == 0.2 and abs(float(row['uz']) - 1e-6) > 1e-12:
                    raise ValueError('Nonzero end displacement failed')
        print('Nonzero physical displacement constraints and reaction work passed')
    else:
        counts = [3, 4, 6, 8, 12] if args.case == 'nonuniform' else [3, 12]
        for suffix in ['solid'] + [str(count) for count in counts]:
            run(args.executable, source, args.directory, args.case + '_' + suffix, [])
        compare(args.directory, args.case, counts)
