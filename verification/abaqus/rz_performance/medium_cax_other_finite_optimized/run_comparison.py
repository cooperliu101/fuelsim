"""Compare existing CAX inputs using isolated, byte-identical cards and meshes."""
import argparse, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
root=a.root.resolve();out=a.results.resolve();out.mkdir(parents=True,exist_ok=True)
bins={'baseline':a.baseline.resolve(),'candidate':a.candidate.resolve()}
env=os.environ.copy();env.pop('PETSC_OPTIONS',None);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',NETCDF_LIBRARY='/home/cooper/miniforge/envs/moose/lib/libnetcdf.so')
models=['cax8t','cax8rt','cax4rt'];records=[];hashes={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in bins.values()}
for model in models:
 for timing in [True,False]:
  card=root/'verification/fuelsim'/f'quasistatic_rz_performance_medium_{model}_finite{"_timing" if timing else ""}.fsi'
  hashes[str(card)]=hashlib.sha256(card.read_bytes()).hexdigest()
  # Read the existing mesh path without changing the card.
  mesh=(card.parent/re.search(r'^  file = (.+)$',card.read_text(),re.M)[1]).resolve()
  hashes[str(mesh)]=hashlib.sha256(mesh.read_bytes()).hexdigest()
  order=[('baseline',True),('candidate',True),('baseline',False),('candidate',False),('candidate',False),('baseline',False)] if timing else [('baseline',False),('candidate',False)]
  expected=None
  for sequence,(variant,warmup) in enumerate(order):
   work=out/model/('timing' if timing else 'accuracy')/variant
   (work/'fuelsim').mkdir(parents=True,exist_ok=True);(work/'meshes').mkdir(exist_ok=True)
   copied=work/'fuelsim'/card.name;shutil.copyfile(card,copied);shutil.copyfile(mesh,work/'meshes'/mesh.name)
   assert copied.read_bytes()==card.read_bytes()
   log=out/f'{model}_{"timing" if timing else "accuracy"}_{sequence}_{variant}.log'
   print('START',log.name,flush=True);start=time.perf_counter()
   with log.open('w') as f:subprocess.run(['taskset','-c','0',str(bins[variant]),'-i',str(copied)],env=env,cwd=work,stdout=f,stderr=subprocess.STDOUT,check=True)
   elapsed=time.perf_counter()-start
   v=dict(re.findall(r'^([^=\n]+)=([^\n]+)$',log.read_text(),re.M));assert v['completed']=='true' and v['mpi_ranks']=='1' and v['accepted_steps']=='20' and v['rejected_steps']=='0'
   engineering={k:v for k,v in v.items() if k.startswith(('region.','contact.'))}
   if expected is None:expected=engineering
   else:assert engineering==expected
   records.append(dict(element=model,timing=timing,variant=variant,warmup=warmup,sequence=sequence,external_seconds=elapsed,internal_seconds=float(v['total_seconds']),iterations=int(v['nonlinear_iterations_total']),log=log.name))
   (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n');print('COMPLETE',elapsed,flush=True)
  if not timing:
   old=out/model/'accuracy/baseline/fuelsim'/f'quasistatic_rz_performance_medium_{model}_finite_results.e';new=out/model/'accuracy/candidate/fuelsim'/old.name
   for file in [old,new]:hashes[str(file)]=hashlib.sha256(file.read_bytes()).hexdigest()
   check=root/'verification/abaqus/rz_performance/medium_cax4t_finite_optimized/check_equivalence.py'
   subprocess.run(['/home/cooper/miniforge/envs/moose/bin/python',str(check),str(old),str(new),'--report',str(out/f'{model}_equivalence.json')],env=env,check=True)
   prefix=root/'verification/abaqus/rz_performance'/f'medium_{model}_finite_incremental'/f'rz_performance_medium_{model}_finite'
   subprocess.run(['/home/cooper/miniforge/envs/moose/bin/python',str(check.parent.parent/'compare.py'),str(new),str(prefix),'--element',model,'--incremental','--report',str(out/f'{model}_accuracy.json')],env=env,check=True)
 (out/'provenance.json').write_text(json.dumps({'baseline_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'sha256':hashes,'cpu':0,'threads':1},indent=2)+'\n')
