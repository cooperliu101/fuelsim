"""CTest driver: copy complete tracked cards unchanged and run the production CLI."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from compare import compare, summary


def run(executable, source, destination, name, launcher):
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ('plate.e', name + '.fsi'):
        shutil.copyfile(source / filename, destination / filename)
        if (source / filename).read_bytes() != (destination / filename).read_bytes():
            raise ValueError('Staged production input differs from its tracked source')
    for suffix in ('_history.csv', '_summary.csv', '.e'):
        (destination / (name + suffix)).unlink(missing_ok=True)
    command = launcher + [str(executable), '-i', str(destination / (name + '.fsi'))]
    process = subprocess.run(command, capture_output=True, text=True,
                             env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    (destination / (name + '.log')).write_text(process.stdout + process.stderr)
    if process.returncode:
        raise RuntimeError(f'{name} failed with {process.returncode}\n{process.stdout}\n{process.stderr}')


def mpi_check(serial, parallel):
    name = 'nonuniform_12'
    left, right = summary(serial / (name + '_summary.csv')), summary(parallel / (name + '_summary.csv'))
    for field in ('mode_count', 'axial_nodes', 'dof_count', 'constraint_count'):
        if left[field] != right[field]:
            raise ValueError(f'MPI changed {field}')
    for field in ('strain_energy', 'external_work'):
        if abs(right[field] / left[field] - 1) > 1e-8:
            raise ValueError(f'MPI changed {field}')
    with (serial / (name + '_history.csv')).open() as stream:
        a = list(csv.DictReader(stream))
    with (parallel / (name + '_history.csv')).open() as stream:
        b = list(csv.DictReader(stream))
    if len(a) != len(b):
        raise ValueError('MPI changed result coverage')
    for kind, fields in [('node', ['ux', 'uy', 'uz']), ('point', ['sxx', 'syy', 'szz', 'sxy', 'syz', 'sxz'])]:
        x = np.array([[float(row[field]) for field in fields] for row in a if row['kind'] == kind])
        y = np.array([[float(row[field]) for field in fields] for row in b if row['kind'] == kind])
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
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    if args.mpiexec:
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
