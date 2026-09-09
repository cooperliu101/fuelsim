#!/usr/bin/env python3
"""Compare two Fuelsim binaries using one unchanged, output-disabled RZ card."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
CARD = ROOT / 'verification/fuelsim/steady_rz_performance_medium_timing.fsi'


def run(simd_on, simd_off, results, cpu):
    results.mkdir(parents=True, exist_ok=True)
    binaries = {'on': simd_on, 'off': simd_off}
    env = os.environ.copy()
    threads = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    env.update(threads)
    env.pop('PETSC_OPTIONS', None)
    files = [CARD, Path(__file__).resolve(), simd_on, simd_off]
    provenance = dict(
        git_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        cpu=cpu, threads=threads,
        sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (results / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    records = []
    # Warm both variants. Reverse the order in the second measured pair.
    for repeat, order in enumerate((('on', 'off'), ('on', 'off'), ('off', 'on'))):
        for variant in order:
            log = results / f'simd_{variant}_run{repeat}.log'
            print(f'Starting SIMD {variant}, run {repeat}', flush=True)
            start = time.perf_counter()
            with log.open('w') as output:
                subprocess.run(['taskset', '-c', str(cpu), str(binaries[variant]), '-i', str(CARD)],
                               cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT, check=True)
            seconds = time.perf_counter() - start
            values = dict(re.findall(r'^([^=\n]+)=([^\n]+)$', log.read_text(), re.M))
            for key, expected in (('completed', 'true'), ('load_steps_completed', '20'), ('mpi_ranks', '1')):
                if values.get(key) != expected:
                    raise RuntimeError(f'{variant}: unexpected {key}')
            records.append(dict(simd=variant, run=repeat, warmup=repeat == 0,
                                external_wall_seconds=seconds,
                                solver_internal_seconds=float(values['total_seconds']),
                                nonlinear_iterations=int(values['nonlinear_iterations_total']),
                                jacobian_evaluations=int(values['jacobian_evaluations_total']),
                                load_predictor_fallbacks=int(values['load_predictor_fallbacks'])))
            (results / 'timing.json').write_text(json.dumps(records, indent=2) + '\n')
            print(f'Completed: {seconds:.6f} seconds', flush=True)
    for name, expected in provenance['sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Measured file changed during timing: ' + name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--simd-on', type=Path, required=True)
    parser.add_argument('--simd-off', type=Path, default=ROOT / 'build/fuelsim')
    parser.add_argument('--results-directory', type=Path,
                        default=ROOT / 'verification/abaqus/rz_performance/medium_simd_off')
    parser.add_argument('--cpu', type=int, default=0)
    args = parser.parse_args()
    run(args.simd_on.resolve(), args.simd_off.resolve(), args.results_directory.resolve(), args.cpu)
