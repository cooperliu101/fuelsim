"""Time unchanged production cards on one CPU; retain the cost of all preprocessing and output."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import subprocess

from compare import compare
from run import stage


def benchmark(executable, directory, cpu, repeats, baseline_executable=None, solid_ends=False):
    source = Path(__file__).resolve().parent
    executable, directory = executable.resolve(), directory.resolve()
    if baseline_executable is not None:
        baseline_executable = baseline_executable.resolve()
    binaries = {'solid': executable, 'modal': executable}
    if baseline_executable is not None:
        binaries['baseline'] = baseline_executable
    hashes = {kind: hashlib.sha256(path.read_bytes()).hexdigest() for kind, path in binaries.items()}
    if cpu not in os.sched_getaffinity(0) or repeats < 1:
        raise ValueError('The benchmark requires an available CPU and a positive repetition count')
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    records = []
    for repeat in range(1, repeats + 1):
        # Alternate order to avoid attributing a fixed warm-cache ordering to the method.
        order = ['solid', 'baseline', 'modal'] if baseline_executable is not None else ['solid', 'modal']
        if repeat % 2 == 0:
            order.reverse()
        for kind in order:
            name, mesh = (('transverse_refined_solid', 'plate_reference.e') if kind == 'solid'
                          else ('transverse_local_12', 'plate_local.e'))
            if solid_ends:
                name, mesh = (('transverse_support_solid', 'plate_support_reference.e') if kind == 'solid'
                              else ('transverse_hybrid_12', 'plate_solid_ends_128_32.e'))
            target = directory / str(repeat) / kind
            stage(source, target, name, mesh)
            timing = target / 'timing.csv'
            command = ['/usr/bin/time', '-o', str(timing), '-f',
                       'elapsed_seconds,%e\npeak_rss_kib,%M\nuser_seconds,%U\nsystem_seconds,%S',
                       'taskset', '-c', str(cpu), str(binaries[kind]), '-i', str(target / (name + '.fsi'))]
            with (target / (name + '.log')).open('w') as log:
                process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment)
            if process.returncode:
                raise RuntimeError(f'Production benchmark failed: {target}')
            with timing.open() as stream:
                values = {key: float(value) for key, value in csv.reader(stream)}
            records.append(dict(repeat=repeat, model=kind, cpu=cpu, executable_sha256=hashes[kind], **values))
            print(records[-1], flush=True)
        # Accuracy is checked outside timed production execution, over every reference point.
        if solid_ends:
            compare(directory / str(repeat) / 'modal', 'transverse_hybrid', [12], 'plate_solid_ends_128_32.e',
                    accuracy_count=12, reference_path=directory / str(repeat) / 'solid' / 'transverse_support_solid.e',
                    diagnostics=True, point_support=True)
        else:
            compare(directory / str(repeat) / 'modal', 'transverse_local', [12], 'plate_local.e', accuracy_count=12,
                    reference_path=directory / str(repeat) / 'solid' / 'transverse_refined_solid.e')
        if baseline_executable is not None:
            compare(directory / str(repeat) / 'baseline', 'transverse_local', [12], 'plate_local.e', accuracy_count=12,
                    reference_path=directory / str(repeat) / 'solid' / 'transverse_refined_solid.e')
    with (directory / 'benchmark.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--cpu', type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--baseline-executable', type=Path)
    parser.add_argument('--solid-ends', action='store_true')
    args = parser.parse_args()
    if args.solid_ends and args.baseline_executable:
        parser.error('The native solid end study compares the same executable in its two production paths')
    benchmark(args.executable, args.directory, args.cpu, args.repeats, args.baseline_executable, args.solid_ends)
