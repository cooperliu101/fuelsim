"""Run identical, unchanged production cards with one and two MPI ranks."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import netCDF4
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--fuelsim',type=Path,required=True)
p.add_argument('--mpiexec',required=True)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--work',type=Path,required=True)
a=p.parse_args()
outputs=[]
for ranks in [1,2]:
    folder=a.work/str(ranks)
    folder.mkdir(parents=True,exist_ok=True)
    for name in ['dc3d20_transient.fsi','dc3d20.e']:
        shutil.copyfile(a.source/name,folder/name)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    run=subprocess.run([a.mpiexec,'-n',str(ranks),str(a.fuelsim.resolve()),'-i',str((folder/'dc3d20_transient.fsi').resolve())],capture_output=True,text=True,env=env)
    (folder/'run.log').write_text(run.stdout+run.stderr)
    if run.returncode:
        raise RuntimeError(run.stdout+run.stderr)
    with netCDF4.Dataset(folder/'dc3d20_transient_results.e') as data:
        outputs.append({name:np.asarray(value[:]) for name,value in data.variables.items() if name.startswith('vals_') or name=='time_whole'})
if outputs[0].keys()!=outputs[1].keys():raise AssertionError('MPI output fields differ')
for name in outputs[0]:
    if not np.allclose(outputs[0][name],outputs[1][name],rtol=1e-11,atol=1e-9,equal_nan=True):
        raise AssertionError('MPI field differs: '+name)
print('All native production output fields match between one and two ranks.')
