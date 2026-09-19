"""Manual matched-input validation and alternating single-core wall-time study."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import time
import netCDF4
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from compare import metric


def identity():
    for case in ['steady','transient','interface']:
        a=(ROOT/(case+'.fsi')).read_text().split('[Outputs]')[0]
        b=(ROOT/(case+'_timing.fsi')).read_text().split('[Outputs]')[0]
        assert a == b,case
        a=(ROOT/(case+'.inp')).read_text().split('*OUTPUT')[0]
        b=(ROOT/(case+'_timing.inp')).read_text().split('*OUTPUT')[0]
        assert a == b,case


def windows(path):
    return subprocess.check_output(['wslpath','-w',str(path.resolve())],text=True).strip()


def native(case,folder,extract):
    folder.mkdir(parents=True,exist_ok=True)
    command=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe',
             '-NoProfile','-ExecutionPolicy','Bypass','-File',windows(ROOT/'run_native.ps1'),
             '-SourceDirectory',windows(ROOT),'-DestinationDirectory',windows(folder),'-CaseName',case]
    if extract:
        command.append('-Extract')
    with (folder/'runner.log').open('w') as stream:
        subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
    result=json.loads((folder/'timing.json').read_text())
    data=(folder/(case+'.dat')).read_text(encoding='latin1')
    summary=data[data.rindex('JOB TIME SUMMARY'):]
    for key,label in [('analysis_cpu_seconds','TOTAL CPU TIME'),('analysis_wall_seconds','WALLCLOCK TIME')]:
        result[key]=float(re.search(label+r' \(SEC\)\s*=\s*([\d.E+-]+)',summary)[1])
    message=(folder/(case+'.msg')).read_text(encoding='latin1')
    for key,label in [('increments','INCREMENTS'),('cutbacks','CUTBACKS IN AUTOMATIC INCREMENTATION'),
                      ('iterations','ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT'),
                      ('factorizations','INVOLVE MATRIX DECOMPOSITION, INCLUDING')]:
        result[key]=int(re.search(r'(\d+)\s+'+label,message[message.rindex('ANALYSIS SUMMARY:'):])[1])
    assert result['increments']==(1 if case.startswith('steady') else 20)
    assert result['cutbacks']==0
    (folder/'timing.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def fuelsim(executable,case,folder):
    folder.mkdir(parents=True,exist_ok=True)
    for previous in folder.glob(case+'_results*.e'):
        previous.unlink()
    for name in [case+'.fsi', 'interface.e' if case.startswith('interface') else 'bulk.e']:
        shutil.copyfile(ROOT/name,folder/name)
    env=os.environ.copy()
    env.update({key:'1' for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})
    with (folder/'fuelsim.log').open('w') as stream:
        start=time.perf_counter()
        subprocess.run(['taskset','-c','0',str(executable),'-i',case+'.fsi'],cwd=folder,env=env,
                       stdout=stream,stderr=subprocess.STDOUT,check=True)
        elapsed=time.perf_counter()-start
    log=(folder/'fuelsim.log').read_text()
    assert 'completed=true' in log
    result={'external_wall_seconds':elapsed}
    for key,value in re.findall(r'^([^=\n]+)=([^\n]+)$',log,re.M):
        try:
            result[key]=float(value)
        except ValueError:
            pass
    assert result.get('accepted_steps',result.get('load_steps_completed'))==(1 if case.startswith('steady') else 20)
    assert result.get('rejected_steps',0)==0 and result.get('rejected_load_steps',0)==0
    if case.endswith('_timing'):
        assert not list(folder.glob('*results*')) and not list(folder.glob('*.checkpoint'))
    (folder/'timing.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def compare(case,folder):
    ref=np.load(folder/'abaqus'/(case+'.npz'))
    with netCDF4.Dataset(folder/'fuelsim'/(case+'_results.e')) as f:
        times=np.asarray(f['time_whole'][:]);indices=np.flatnonzero(times>0)
        assert np.array_equal(times[indices].astype(np.float32),ref['time'].astype(np.float32))
        names=netCDF4.chartostring(f['name_nod_var'][:]).tolist()
        t=np.asarray(f['vals_nod_var'+str(names.index('temperature')+1)][:])[indices]
        r=np.asarray(f['vals_nod_var'+str(names.index('heat_reaction')+1)][:])[indices]
        enames=netCDF4.chartostring(f['name_elem_var'][:]).tolist()
        def field(name):
            k=enames.index(name)+1
            return np.concatenate([np.asarray(f[f'vals_elem_var{k}eb{b+1}'][:])[indices]
                                   for b in range(len(f.dimensions['num_el_blk']))],axis=1)
        xyz=np.stack([np.stack([field(f'point_{c}_q{q}')[0] for c in 'xyz'],axis=-1)
                      for q in range(8)],axis=1)
        distances=np.linalg.norm(xyz[:,:,None,:]-ref['coordinates'][:,None,:,:],axis=-1)
        mapping=np.argmin(distances,axis=2)
        assert np.max(np.min(distances,axis=2))<1e-8
        assert np.all(np.sort(mapping,axis=1)==np.arange(8))
        fields={'temperature':metric(t,ref['temperature']),
                'temperature_rise':metric(t-300,ref['temperature']-300),
                'heat_reaction':metric(r,ref['reaction'])}
        for d,c in enumerate('xyz'):
            actual=np.stack([field(f'heat_flux_{c}_q{q}') for q in range(8)],axis=2)
            reference=np.take_along_axis(ref['flux'][:,:,:,d],mapping[None,:,:],axis=2)
            fields['heat_flux_'+c]=metric(actual,reference,zero=c!='x')
    result={'passed':all(v['passed'] for v in fields.values()),'fields':fields,
            'frames':len(indices),'nodes_per_frame':t.shape[1],'points_per_frame':int(np.prod(xyz.shape[:2]))}
    (folder/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:
        raise AssertionError(json.dumps(result))
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fuelsim',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mode',choices=['validate','check','timing'],required=True)
    p.add_argument('--cases',nargs='+',default=['steady','transient','interface'])
    args=p.parse_args();args.work=args.work.resolve();args.fuelsim=args.fuelsim.resolve()
    identity()
    meshes={'interface.e' if case=='interface' else 'bulk.e' for case in args.cases}
    if any(not (ROOT/name).is_file() for name in meshes):
        p.error('Generate local geometry first: python verification/thermal/performance/mesh.py')
    if args.mode=='check' and any(not (ROOT/'references'/(case+'.npz')).is_file() for case in args.cases):
        p.error('Full native references are local-only; use --mode validate with Abaqus on a fresh checkout')
    input_hashes={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in ROOT.iterdir()
                  if path.suffix in ['.fsi','.inp','.e']}
    executable_hash=hashlib.sha256(args.fuelsim.read_bytes()).hexdigest()
    if args.mode=='timing':
        validated=json.loads((args.work/'validate_provenance.json').read_text())
        assert validated['inputs']==input_hashes, 'Inputs changed since native validation'
        assert validated['executable_sha256']==executable_hash, 'Executable changed since native validation'
    results={}
    for case in args.cases:
        if args.mode in ['validate','check']:
            folder=args.work/'validation'/case
            print('validate '+case,flush=True)
            fuelsim(args.fuelsim,case,folder/'fuelsim')
            if args.mode=='validate':
                native(case,folder/'abaqus',True)
            else:
                (folder/'abaqus').mkdir(parents=True,exist_ok=True)
                shutil.copyfile(ROOT/'references'/(case+'.npz'),folder/'abaqus'/(case+'.npz'))
            results[case]=compare(case,folder)
            print(case+' accuracy passed',flush=True)
        else:
            assert json.loads((args.work/'validation'/case/'comparison.json').read_text())['passed']
            samples=[]
            for i in range(4):
                # Sample zero warms caches and is explicitly excluded from statistics.
                folder=args.work/'timing'/case/str(i)
                sample={}
                order=['fuelsim','abaqus'] if i%2==0 else ['abaqus','fuelsim']
                for solver in order:
                    print(f'{case} sample {i} {solver}',flush=True)
                    sample[solver]=(fuelsim(args.fuelsim,case+'_timing',folder/solver) if solver=='fuelsim'
                                    else native(case+'_timing',folder/solver,False))
                samples.append(sample)
            results[case]={'warmup':samples[0],'samples':samples[1:]}
            for solver,keys in [('fuelsim',['nonlinear_iterations_total','linear_iterations_total',
                                          'residual_evaluations_total','jacobian_evaluations_total']),
                                ('abaqus',['increments','cutbacks','iterations','factorizations'])]:
                for key in keys:
                    assert len({s[solver][key] for s in samples})==1,(case,solver,key)
            for solver in ['fuelsim','abaqus']:
                values=[s[solver]['external_wall_seconds'] for s in samples[1:]]
                results[case][solver+'_external']={'median':statistics.median(values),'min':min(values),'max':max(values)}
            results[case]['abaqus_over_fuelsim_external_ratio']=(results[case]['abaqus_external']['median']/
                                                               results[case]['fuelsim_external']['median'])
        (args.work/(args.mode+'_summary.json')).write_text(json.dumps(results,indent=2)+'\n')
    provenance={'executable':str(args.fuelsim),'executable_sha256':executable_hash,
                'linux':subprocess.check_output(['uname','-a'],text=True).strip(),
                'cpu':subprocess.check_output(['lscpu'],text=True),
                'inputs':input_hashes}
    (args.work/(args.mode+'_provenance.json')).write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__':
    main()
