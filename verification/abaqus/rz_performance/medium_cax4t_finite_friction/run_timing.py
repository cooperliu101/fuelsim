"""Time the committed friction input using one CPU and one MUMPS process."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--executable', type=Path, required=True)
parser.add_argument('--results', type=Path, required=True)
parser.add_argument('--element', choices=('cax4t', 'cax4rt', 'cax8t', 'cax8rt'), default='cax4t')
args = parser.parse_args()
root = args.root.resolve()
executable = args.executable.resolve()
results = args.results.resolve()
results.mkdir(parents=True, exist_ok=True)
card = root / f'verification/fuelsim/quasistatic_rz_performance_medium_{args.element}_finite_friction_timing.fsi'
env = os.environ.copy()
env.pop('PETSC_OPTIONS', None)
env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
records = []
expected = None
for run in range(3):
    log = results / ('fuelsim_run%d.log' % run)
    print('START', run, flush=True)
    start = time.perf_counter()
    with log.open('w') as stream:
        subprocess.run(['taskset', '-c', '0', str(executable), '-i', str(card)],
                       cwd=root, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
    elapsed = time.perf_counter() - start
    values = dict(re.findall(r'^([^=\n]+)=([^\n]+)$', log.read_text(), re.M))
    assert values['completed'] == 'true' and values['accepted_steps'] == '20'
    assert values['rejected_steps'] == '0' and values['mpi_ranks'] == '1'
    fields = {key: value for key, value in values.items()
              if key.startswith(('region.', 'contact.')) or key == 'nonlinear_iterations_total'}
    if expected is None:
        expected = fields
    else:
        assert fields == expected
    records.append(dict(run=run, warmup=run == 0, external_seconds=elapsed,
                        internal_seconds=float(values['total_seconds']),
                        nonlinear_iterations=int(values['nonlinear_iterations_total'])))
    (results / 'fuelsim_timing.json').write_text(json.dumps(records, indent=2) + '\n')
    print('COMPLETE', elapsed, flush=True)
(results / 'timing_provenance.json').write_text(json.dumps({
    'cpu': 0, 'threads': 1, 'sha256': {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (executable, card)}}, indent=2) + '\n')
