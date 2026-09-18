"""Repeat identical checkpoint continuations and report all differing fields."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import netCDF4
import numpy as np

ROOT=Path(__file__).resolve().parent


def compare(actual,reference):
    report={};failed=[]
    with netCDF4.Dataset(actual) as a,netCDF4.Dataset(reference) as r:
        times=np.asarray(a['time_whole'][:]);rt=np.asarray(r['time_whole'][:])
        indices=[int(np.argmin(abs(rt-t))) for t in times]
        np.testing.assert_allclose(rt[indices],times,atol=1e-14,rtol=0)
        names={family:list(netCDF4.chartostring(r['name_'+family+'_var'][:])) for family in ('nod','elem','glo')}
        for key in r.variables:
            if not key.startswith(('vals_nod_var','vals_elem_var','vals_glo_var')):continue
            av=np.asarray(a[key][:]);rv=np.asarray(r[key][indices]);valid=np.isfinite(rv)
            assert np.array_equal(np.isnan(av),np.isnan(rv))
            if key=='vals_glo_var':
                columns=[(name,av[:,i],rv[:,i],valid[:,i]) for i,name in enumerate(names['glo'])]
            else:
                family='nod' if key.startswith('vals_nod_var') else 'elem'
                index=int(key.split('_var')[1].split('eb')[0])-1
                columns=[(names[family][index]+('/block'+key.split('eb')[1] if family=='elem' else ''),av,rv,valid)]
            for name,x,y,mask in columns:
                if not np.any(mask):continue
                difference=abs(x-y);limit=1e-12+1e-12*abs(y)
                bad=mask&(difference>limit)
                report[name]=dict(maximum_absolute=float(np.max(difference[mask])),failed_count=int(sum(bad.flatten())))
                if np.any(bad):
                    where=np.argwhere(bad)
                    failed.append(dict(field=name,count=len(where),examples=[dict(
                        time=float(times[row[0]]),index=row.tolist(),actual=float(x[tuple(row)]),
                        reference=float(y[tuple(row)]),difference=float(difference[tuple(row)])) for row in where[:5]]))
    return dict(status='failed' if failed else 'passed',fields=report,failures=failed)


def run(args):
    source=ROOT.parent;args.work.mkdir(parents=True,exist_ok=True)
    outputs={}
    for ranks,repeat in [(1,0),(1,1),(2,0),(4,0),(4,1)]:
        work=args.work/f'rank{ranks}_repeat{repeat}';work.mkdir(exist_ok=True)
        for name in ['contact_law.e','coupled_contact_restart.fsi']:
            shutil.copyfile(source/name,work/name)
            assert (source/name).read_bytes()==(work/name).read_bytes()
        shutil.copyfile(args.checkpoint,work/'coupled_contact_split.chk')
        command=[str(args.executable.resolve()),'-i','coupled_contact_restart.fsi']
        if ranks>1:command=[str(args.mpiexec),'-n',str(ranks)]+command
        result=subprocess.run(command,cwd=work,capture_output=True,text=True)
        (work/'production.log').write_text(result.stdout+result.stderr)
        assert result.returncode==0 and 'completed=true' in result.stdout,result.stdout+result.stderr
        outputs[f'rank{ranks}_repeat{repeat}']=Path(next(l.split('=',1)[1] for l in result.stdout.splitlines()
                                                      if l.startswith('results_file=')))
    comparisons={}
    for name,path in outputs.items():
        comparisons[name+'_vs_serial']=compare(path,outputs['rank1_repeat0'])
    comparisons['rank4_repeatability']=compare(outputs['rank4_repeat1'],outputs['rank4_repeat0'])
    report=dict(executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),
                checkpoint_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),comparisons=comparisons)
    (ROOT/'mpi_differences.json').write_text(json.dumps(report,indent=2)+'\n')
    for name,r in comparisons.items():print(name,r['status'],len(r['failures']))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--mpiexec',type=Path,required=True)
    p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    run(p.parse_args())
