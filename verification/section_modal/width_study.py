"""Run complete production cards; failed 1% comparisons remain explicit research results."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import subprocess

from compare import compare, summary, failed_metrics
from run import stage, run, mpi_check


def run_reference(executable, directory, cpu):
    source = Path(__file__).resolve().parent
    executable, directory = executable.resolve(), directory.resolve()
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    records = []
    for case in ('nonuniform', 'axial', 'bending', 'transverse'):
        name = case + '_support_solid'
        target = directory / case
        stage(source, target, name, 'plate_support_reference.e')
        timing = target / 'timing.csv'
        command = ['/usr/bin/time', '-o', str(timing), '-f', 'elapsed_seconds,%e\npeak_rss_kib,%M',
                   'taskset', '-c', str(cpu), str(executable), '-i', str(target / (name + '.fsi'))]
        with (target / (name + '.log')).open('w') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment)
        if completed.returncode:
            raise RuntimeError(f'Point-supported solid reference failed; preserve {target}')
        with timing.open() as stream:
            times = {key: float(value) for key, value in csv.reader(stream)}
        records.append(dict(case=case, cpu=cpu, **times,
                            executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                            reference_sha256=hashlib.sha256((target / (name + '.e')).read_bytes()).hexdigest()))
        with (directory / 'reference_runs.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
        print(f'Point-supported solid {case}: {times}', flush=True)


def study(executable, directory, reference_root, cpu):
    executable, directory, reference_root = executable.resolve(), directory.resolve(), reference_root.resolve()
    source = Path(__file__).resolve().parent
    if cpu not in os.sched_getaffinity(0):
        raise ValueError('The selected benchmark CPU is unavailable')
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    digest = hashlib.sha256(executable.read_bytes()).hexdigest()
    records = []
    jobs = [(case, 1, f'width{strips}', count) for case in ('axial', 'bending', 'transverse')
            for strips, count in ((2, 38), (4, 64))]
    order = [('support', 12), ('width1', 25), ('width2', 38), ('width4', 64)]
    for repeat in (1, 2):
        jobs.extend(('nonuniform', repeat, kind, count) for kind, count in (order if repeat == 1 else order[::-1]))
    for case, repeat, kind, count in jobs:
        label = case + '_' + kind
        name = label + '_' + str(count)
        target = directory / case / str(repeat) / kind
        mesh = 'plate_width.e'
        stage(source, target, name, mesh)
        timing = target / 'timing.csv'
        command = ['/usr/bin/time', '-o', str(timing), '-f', 'elapsed_seconds,%e\npeak_rss_kib,%M',
                   'taskset', '-c', str(cpu), str(executable), '-i', str(target / (name + '.fsi'))]
        with (target / (name + '.log')).open('w') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment)
        if completed.returncode:
            raise RuntimeError(f'Production solve failed; preserve {target}')
        with timing.open() as stream:
            times = {key: float(value) for key, value in csv.reader(stream)}
        reference = reference_root / case / (case + '_support_solid.e')
        measured = compare(target, label, [count], mesh, reference_path=reference,
                           diagnostics=True, point_support=True)[0]
        failures = failed_metrics(measured)
        report = summary(target / (name + '_summary.csv'))
        record = dict(repeat=repeat, model=kind, cpu=cpu, executable_sha256=digest,
                      reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
                      global_system_size=int(report['global_system_size']),
                      local_system_size=int(report['local_system_size']),
                      accuracy_status='FAIL' if failures else 'PASS', failed_metrics=';'.join(failures),
                      **times, **measured)
        records.append(record)
        directory.mkdir(parents=True, exist_ok=True)
        with (directory / 'study.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
        print(f"{case} {kind} repeat={repeat}: {times}, accuracy={record['accuracy_status']}, "
              f"failed={record['failed_metrics']}", flush=True)
    for case in ('axial', 'bending', 'transverse', 'nonuniform'):
        selected = [row for row in records if row['repeat'] == 1 and row['case'].startswith(case + '_width')]
        selected.sort(key=lambda row: row['modes'])
        for left, right in zip(selected, selected[1:]):
            if right['modal_energy'] < left['modal_energy'] * (1 - 1e-8):
                raise ValueError('Nested width enrichment reduced fixed-load compliance')
    return all(row['accuracy_status'] == 'PASS' for row in records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('reference_root', type=Path)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--mpiexec', type=Path)
    parser.add_argument('--reference-only', action='store_true')
    args = parser.parse_args()
    if args.reference_only:
        if args.mpiexec:
            raise ValueError('Reference timing requires one process')
        run_reference(args.executable, args.reference_root, args.cpu)
    elif args.mpiexec:
        source = Path(__file__).resolve().parent
        parallel = args.directory / 'mpi'
        run(args.executable.resolve(), source, parallel.resolve(), 'nonuniform_width4_64',
            [str(args.mpiexec), '-n', '2'], 'plate_width.e')
        mpi_check(args.directory / 'nonuniform/1/width4', parallel, 'nonuniform_width4_64')
    else:
        raise SystemExit(0 if study(args.executable, args.directory, args.reference_root, args.cpu) else 1)
