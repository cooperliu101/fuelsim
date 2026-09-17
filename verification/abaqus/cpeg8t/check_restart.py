"""Run static complete cards; compare restarted production output at matching times."""
import argparse
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from check_thermal_mass import check


def equivalent(actual_path, reference_path, restart=False):
    with netCDF4.Dataset(actual_path) as actual, netCDF4.Dataset(reference_path) as reference:
        times=np.asarray(actual['time_whole'][:])
        ref_times=np.asarray(reference['time_whole'][:])
        indices=[int(np.argmin(abs(ref_times-t))) for t in times]
        np.testing.assert_allclose(ref_times[indices],times,rtol=0,atol=1e-14)
        if restart and (times[0]!=.5 or times[-1]!=1):
            raise AssertionError('Restart must continue from the middle of the history')
        for name in reference.variables:
            if not name.startswith(('vals_nod_var','vals_elem_var','vals_glo_var')):
                continue
            a,r=np.asarray(actual[name][:]),np.asarray(reference[name][indices])
            np.testing.assert_array_equal(np.isnan(a),np.isnan(r))
            valid=np.isfinite(r)
            np.testing.assert_allclose(a[valid],r[valid],rtol=1e-12,atol=1e-12,err_msg=name)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--mpiexec',type=Path)
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    args.work.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(source/'distorted.e',args.work/'distorted.e')
    for case in ['thermal_mass','thermal_mass_split','thermal_mass_restart']:
        shutil.copyfile(source/(case+'.fsi'),args.work/(case+'.fsi'))
        assert (source/(case+'.fsi')).read_bytes()==(args.work/(case+'.fsi')).read_bytes()
        command=[str(args.executable.resolve()),'-i',case+'.fsi']
        # A serial checkpoint must also be restartable with four MPI processes.
        if args.mpiexec and case=='thermal_mass_restart':
            command=[str(args.mpiexec),'-n','4']+command
        run=subprocess.run(command,cwd=args.work,capture_output=True,text=True)
        (args.work/(case+'.log')).write_text(run.stdout+run.stderr)
        if run.returncode or 'completed=true' not in run.stdout:
            raise RuntimeError(run.stdout+run.stderr)
    check(args.work)
    equivalent(args.work/'thermal_mass_restart_results.part1.e',args.work/'thermal_mass_results.e',True)
    equivalent(args.work/'thermal_mass_split_results.e',args.work/'thermal_mass_results.e')
    print('All nodal fields, material tensors, scalar histories and section results match uninterrupted production')
