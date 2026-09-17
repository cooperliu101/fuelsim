"""Compare a production curved interface with independently solved parabola geometry."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from compare import compare, verify_references


def check(directory, source):
    reports=[]
    with netCDF4.Dataset(directory/'curved_contact_results.e') as d:
        g=dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][-1]))
        names=list(netCDF4.chartostring(d['name_nod_var'][:]))
        total_heat=0.
        for q,(xi,w) in enumerate(zip([-np.sqrt(3/5),0,np.sqrt(3/5)],[5/9,8/9,5/9])):
            x=.008*xi
            y=10*x*x-.0001
            roots=np.roots([200,0,1-20*y,-x])
            real=[r.real for r in roots if abs(r.imag)<1e-14 and -.01<=r.real<=.01]
            projected=min(real,key=lambda r:(r-x)**2+(10*r*r-y)**2)
            normal=np.array([-20*projected,1.])/np.hypot(20*projected,1)
            gap=np.dot([x-projected,y-10*projected**2],normal)
            area=.008*w*.1*np.hypot(20*x,1)
            pressure=-1e9*gap
            heat=1000*100*area
            total_heat+=heat
            for name,expected,atol in [('gap',gap,1e-13),('pressure',0.,1e-7),('area',area,1e-13),
                                       ('force',0.,1e-8),('heat_rate',heat,1e-8)]:
                np.testing.assert_allclose(g[f'contact_0_q{q+3}_{name}'],expected,rtol=1e-11,atol=atol)
            reports.append(dict(gap=float(gap),area=float(area),pressure=float(pressure),
                                projected_x=float(projected)))
        fy=np.asarray(d[f'vals_nod_var{names.index("reaction_y")+1}'][-1])
        fx=np.asarray(d[f'vals_nod_var{names.index("reaction_x")+1}'][-1])
        heat=np.asarray(d[f'vals_nod_var{names.index("heat_reaction")+1}'][-1])
        np.testing.assert_allclose([sum(fy),sum(fx)], [0,0],rtol=0,atol=1e-8)
        np.testing.assert_allclose([np.nansum(heat[:8]),np.nansum(heat[8:])],[-total_heat,total_heat],
                                   rtol=1e-10,atol=1e-8)
        native=json.loads((source/'curved_contact_fields.json').read_text())['steps']['CONTACT'][-1]['fields']
        native_report,failures={},[]
        # All zero masks follow the prescribed rigid translations, reflection
        # symmetry, and the list of nodes actually lying on each interface.
        for field,key,component,active,zeros in [
            ('reaction_x','RF',0,range(1,17),lambda n:n not in (3,4,9,10)),
            ('reaction_y','RF',1,range(1,17),lambda n:n not in (3,4,7,9,10,13)),
            ('heat_reaction','RFL11',None,[1,2,3,4,9,10,11,12],lambda n:n not in (3,4,9,10))]:
            ref={v['nodeLabel']:v['data'] for v in native[key]['values']}
            expected=[ref[n] if component is None else ref[n][component] for n in active]
            actual=np.asarray(d[f'vals_nod_var{names.index(field)+1}'][-1])[np.array(list(active))-1]
            try:
                compare(actual,expected,[zeros(n) for n in active],1e-8,field,native_report)
            except AssertionError as error:
                failures.append(str(error))
        result=dict(analytical_status='passed',analytical=reports,
                    native_status='failed' if failures else 'passed',native=native_report,native_failures=failures)
    (directory/'curved_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Curved thermal projection and heat rates, total force and heat conservation passed analytical checks')
    print(f'Native curved contact comparison: {len(failures)} failed fields')
    return failures


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--require-native',action='store_true')
    args=p.parse_args()
    source=Path(__file__).resolve().parent
    verify_references(source, 'extended_reference.sha256')
    args.work.mkdir(parents=True,exist_ok=True)
    for name in ['curved_contact.fsi','curved.e']:
        shutil.copyfile(source/name,args.work/name)
        assert (source/name).read_bytes()==(args.work/name).read_bytes()
    r=subprocess.run([str(args.executable.resolve()),'-i','curved_contact.fsi'],cwd=args.work,capture_output=True,text=True)
    (args.work/'production.log').write_text(r.stdout+r.stderr)
    if r.returncode or 'completed=true' not in r.stdout:
        raise RuntimeError(r.stdout+r.stderr)
    failures=check(args.work,source)
    if args.require_native and failures:
        raise AssertionError('\n'.join(failures))
