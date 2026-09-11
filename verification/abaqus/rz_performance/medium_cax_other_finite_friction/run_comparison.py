"""Run committed finite-strain friction cards; validate accuracy before timing."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
p.add_argument('--executable', type=Path, required=True)
p.add_argument('--results', type=Path, required=True)
p.add_argument('--elements', nargs='+', choices=('cax4rt', 'cax8t', 'cax8rt'), default=['cax4rt', 'cax8t', 'cax8rt'])
p.add_argument('--phase', choices=('accuracy', 'timing', 'all'), default='all')
a = p.parse_args()
root, executable, results = a.root.resolve(), a.executable.resolve(), a.results.resolve()
ab = root / 'verification/abaqus/rz_performance'
env = os.environ.copy()
env.pop('PETSC_OPTIONS', None)
env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
           NETCDF_LIBRARY='/home/cooper/miniforge/envs/moose/lib/libnetcdf.so')
python = '/home/cooper/miniforge/envs/moose/bin/python'


def run(command, log):
    print('START', str(log), flush=True)
    with log.open('w') as stream:
        subprocess.run(command, env=env, cwd=root, stdout=stream, stderr=subprocess.STDOUT, check=True)
    print('COMPLETE', str(log), flush=True)


def windows(path):
    return r'\\wsl.localhost\Ubuntu' + str(path).replace('/', '\\')


def abaqus(model, output, timing=False):
    output.mkdir(parents=True, exist_ok=True)
    command = ['powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', windows(ab / 'run.ps1'),
               '-SourceDirectory', windows(ab), '-Size', 'medium', '-Element', model, '-Strain', 'finite',
               '-Friction', '-Runs', '3' if timing else '1', '-ResultsDirectory', windows(output)]
    if timing:
        command += ['-Timing']
    run(command, output / 'abaqus_launch.log')


if a.phase in ('accuracy', 'all'):
    for model in a.elements:
        work = results / model / 'accuracy'
        (work / 'fuelsim').mkdir(parents=True, exist_ok=True)
        (work / 'meshes').mkdir(exist_ok=True)
        card = root / 'verification/fuelsim' / f'quasistatic_rz_performance_medium_{model}_finite_friction.fsi'
        mesh = (card.parent / re.search(r'^  file = (.+)$', card.read_text(), re.M)[1]).resolve()
        copied = work / 'fuelsim' / card.name
        shutil.copyfile(card, copied)
        shutil.copyfile(mesh, work / 'meshes' / mesh.name)
        assert copied.read_bytes() == card.read_bytes()
        run(['taskset', '-c', '0', str(executable), '-i', str(copied)], work / 'fuelsim.log')
        abaqus(model, work / 'abaqus')
        output = work / 'fuelsim' / f'quasistatic_rz_performance_medium_{model}_finite_friction_results.e'
        prefix = work / 'abaqus' / f'rz_performance_medium_{model}_finite_friction'
        run([python, str(ab / 'compare.py'), str(output), str(prefix), '--element', model, '--incremental',
             '--friction', '--report', str(work / 'accuracy.json')], work / 'comparison.log')

if a.phase in ('timing', 'all'):
    for model in a.elements:
        metrics = json.loads((results / model / 'accuracy/accuracy.json').read_text())
        if not all(value['passed'] for value in metrics.values()):
            raise RuntimeError('Accuracy must pass before timing ' + model)
        work = results / model / 'timing'
        work.mkdir(parents=True, exist_ok=True)
        run([python, str(ab / 'medium_cax4t_finite_friction/run_timing.py'), '--root', str(root),
             '--executable', str(executable), '--results', str(work), '--element', model], work / 'fuelsim_launch.log')
        abaqus(model, work / 'abaqus', True)
