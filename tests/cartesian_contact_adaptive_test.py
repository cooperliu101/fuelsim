"""Run an unchanged complete input card and inspect production progress/summary."""
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

executable, source, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
(work / 'verification/fuelsim').mkdir(parents=True, exist_ok=True)
(work / 'verification/meshes').mkdir(parents=True, exist_ok=True)
card = work / 'verification/fuelsim/transient_cartesian_contact_adaptive.fsi'
shutil.copyfile(source / 'verification/fuelsim' / card.name, card)
shutil.copyfile(source / 'verification/meshes/b526_friction_reversal.e',
                work / 'verification/meshes/b526_friction_reversal.e')
run = subprocess.run([str(executable), '-i', str(card)], capture_output=True, text=True, timeout=90)
(work / 'production.log').write_text(run.stdout + run.stderr)
assert run.returncode == 0, run.stdout + run.stderr
with card.with_name('transient_cartesian_contact_adaptive_summary.csv').open() as stream:
    summary = dict(list(csv.reader(stream))[1:])
assert int(summary['time_error_rejections']) == 1
assert int(summary['rejected_steps']) == 1
assert summary['last_rejected.failure_category'] == 'time_discretization'
prefix = 'last_rejected.time_error.'
components = {k[len(prefix):]: float(v) for k, v in summary.items() if k.startswith(prefix)}
assert math.isfinite(components['contact_friction']) and components['contact_friction'] > 1
assert all(math.isfinite(v) and v < 1 for k, v in components.items() if k != 'contact_friction'), components
steps = []
for line in run.stdout.splitlines():
    if '=' not in line:
        continue
    key, value = line.split('=', 1)
    if key == 'progress.accepted_steps':
        steps.append({})
    if steps and key.startswith('progress.'):
        steps[-1][key] = value
assert len(steps) >= 2
first = steps[0]
assert float(first['progress.time_step']) < float(summary['last_rejected.time_step'])
assert float(first['progress.time_error.contact_friction']) < 1
assert all(float(step['progress.time_error_estimate']) <= 1 for step in steps)
assert float(steps[-1]['progress.time']) > float(first['progress.time'])
assert abs(float(summary['committed_time']) - 0.4) < 1e-12
report = {'rejected_components': components, 'accepted_progress': steps, 'committed_time': summary['committed_time']}
(work / 'contact_adaptive_evidence.json').write_text(json.dumps(report, indent=2) + '\n')
print('Contact error alone rejected the coarse step; smaller steps were accepted and completed the run.')
