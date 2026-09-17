"""Manual spatial and temporal refinement study using complete static cards."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from check_thermal_mass import check


def run(case,mesh,work,source,executable):
    work.mkdir(parents=True,exist_ok=True)
    for name in [case+'.fsi',mesh]:
        shutil.copyfile(source/name,work/name)
        assert (source/name).read_bytes()==(work/name).read_bytes()
    r=subprocess.run([str(executable.resolve()),'-i',case+'.fsi'],cwd=work,capture_output=True,text=True)
    (work/'production.log').write_text(r.stdout+r.stderr)
    if r.returncode or 'completed=true' not in r.stdout:
        raise RuntimeError(r.stdout+r.stderr)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    args=p.parse_args()
    source=Path(__file__).resolve().parent
    spatial=[]
    for count in [2,4,8]:
        case=f'mesh_convergence_{count}'
        work=args.work/case
        run(case,f'convergence_{count}.e',work,source,args.executable)
        with netCDF4.Dataset(work/(case+'_results.e')) as d:
            names=list(netCDF4.chartostring(d['name_nod_var'][:]))
            temperature=np.asarray(d[f'vals_nod_var{names.index("temperature")+1}'][-1])
            y=np.asarray(d['coordy'][:])
            valid=np.isfinite(temperature)
            # k=10, Q=1000, T(0)=300, k*T'(1)=10*(350-T(1)).
            exact=300+100*y-50*y*y
            np.testing.assert_allclose(temperature[valid],exact[valid],rtol=1e-11,atol=1e-9)
            elements=list(netCDF4.chartostring(d['name_elem_var'][:]))
            a,r=[],[]
            for j in range(3):
                eta=[-np.sqrt(3/5),0,np.sqrt(3/5)][j]
                yy=(np.repeat(np.arange(count),count)+.5*(1+eta))/count
                for i in range(3):
                    var=elements.index(f'stress_xx_q{3*j+i}')+1
                    a.extend(d[f'vals_elem_var{var}eb1'][-1])
                    r.extend(-20*(100*yy-50*yy*yy))
            a,r=np.asarray(a),np.asarray(r)
            spatial.append(dict(elements=count**2,maximum_absolute=float(max(abs(a-r))),
                                relative_l2=float(np.linalg.norm(a-r)/np.linalg.norm(r)),
                                maximum_pointwise_relative=float(max(abs((a-r)/r)))))
    for coarse,fine in zip(spatial[:-1],spatial[1:]):
        np.testing.assert_allclose(coarse['maximum_absolute']/fine['maximum_absolute'],4,rtol=1e-7,atol=1e-7)
    temporal=[]
    for case in ['thermal_mass','time_convergence','thermal_mass_adaptive']:
        work=args.work/case
        run(case,'distorted.e',work,source,args.executable)
        check(work,case,case=='thermal_mass_adaptive')
        with netCDF4.Dataset(work/(case+'_results.e')) as d:
            times=np.asarray(d['time_whole'][:])
            names=list(netCDF4.chartostring(d['name_nod_var'][:]))
            temperature=np.asarray(d[f'vals_nod_var{names.index("temperature")+1}'][-1])[:4]
            # Integral cp(T)dT=Q*t/rho gives the continuous exact solution.
            exact=300+2000/(100+np.sqrt(10200))
            error=float(max(abs(temperature-exact)))
            effective_step=float(sum(np.diff(times)**2)/(2 if case=='thermal_mass_adaptive' else 1))
            temporal.append(dict(case=case,accepted_steps=len(times)-1,
                                 effective_step=effective_step,temperature_absolute_error=error))
    rate=np.log(temporal[0]['temperature_absolute_error']/temporal[1]['temperature_absolute_error'])/np.log(4)
    if not .95<rate<1.05:
        raise AssertionError(f'Backward Euler temporal convergence is not first order: {rate}')
    report=dict(spatial=spatial,temporal=temporal,temporal_order=float(rate))
    (args.work/'convergence.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
