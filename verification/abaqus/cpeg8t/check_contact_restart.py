"""Checkpoint and failed-step recovery with contact and both inelastic histories."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from check_contact_supplement import run
from check_restart import equivalent


def output_path(work, case):
    log = (work/(case+'.log')).read_text()
    return Path(next(line.split('=',1)[1] for line in log.splitlines() if line.startswith('results_file=')))


def check(source, work, executable, mpiexec=None):
    for case in ['coupled_contact','coupled_contact_split']:
        if not run(case, source, work, executable):
            raise AssertionError(case+' did not complete')
    checkpoint = work/'coupled_contact_split.chk'
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    case = 'coupled_contact_failure'
    if run(case, source, work, executable):
        raise AssertionError('The one-iteration failure probe unexpectedly completed')
    log = (work/(case+'.log')).read_text()
    for text in ['committed_time=5.000000000000e-01','accepted_steps=0','rejected_steps=1',
                 'last_rejected.attempted_end_time=5.625000000000e-01',
                 'last_rejected.convergence_reason=DIVERGED_MAX_IT']:
        if text not in log:
            raise AssertionError('Missing failed-step evidence: '+text)
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest() != digest:
        raise AssertionError('Failed continuation changed the input checkpoint')
    equivalent(output_path(work,case),output_path(work,'coupled_contact'))
    if not run('coupled_contact_restart',source,work,executable):
        raise AssertionError('Continuation did not complete')
    equivalent(output_path(work,'coupled_contact_restart'),output_path(work,'coupled_contact'),True)
    equivalent(output_path(work,'coupled_contact_split'),output_path(work,'coupled_contact'))
    report = dict(status='passed',failed_trial_end=.5625,restored_time=.5,
                  all_exported_nodal_material_and_contact_fields=True,
                  mpi_ranks=1,checkpoint_sha256=digest)
    (work/'restart_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    if mpiexec:
        parallel = work/'parallel'
        parallel.mkdir(exist_ok=True)
        for name in ['coupled_contact_split.chk','contact_law.e','coupled_contact_restart.fsi']:
            shutil.copyfile(work/name,parallel/name)
        command = [str(mpiexec),'-n','4',str(executable.resolve()),'-i','coupled_contact_restart.fsi']
        r = subprocess.run(command,cwd=parallel,capture_output=True,text=True)
        (parallel/'coupled_contact_restart.log').write_text(r.stdout+r.stderr)
        if r.returncode or 'completed=true' not in r.stdout:
            raise AssertionError('Four-process continuation failed')
        try:
            equivalent(output_path(parallel,'coupled_contact_restart'),output_path(work,'coupled_contact'),True)
        except AssertionError as error:
            (work/'restart_mpi_comparison.json').write_text(json.dumps(
                dict(status='failed',mpi_ranks=4,reason=str(error)),indent=2)+'\n')
            raise
        report['mpi_ranks'] = 4
        (work/'restart_mpi_comparison.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--mpiexec',type=Path)
    args = p.parse_args()
    check(Path(__file__).resolve().parent,args.work,args.executable,args.mpiexec)
    print('Contact transition, failed trial recovery and complete checkpoint histories agree')
