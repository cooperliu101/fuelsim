"""Alternate baseline/candidate production runs; inputs are never modified."""
import argparse,hashlib,json,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
a.results.mkdir(parents=True,exist_ok=True)
binaries={'baseline':a.baseline.resolve(),'candidate':a.candidate.resolve()}
env=os.environ.copy();env.pop('PETSC_OPTIONS',None);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
cards={'finite':a.root/'verification/fuelsim/quasistatic_rz_performance_medium_cax4t_finite_timing.fsi','small':a.root/'verification/fuelsim/steady_rz_performance_medium_timing.fsi'}
paths=list(binaries.values())+list(cards.values())
(a.results/'provenance.json').write_text(json.dumps({'git_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.root,text=True).strip(),'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'order':'warm baseline, warm candidate, baseline, candidate, candidate, baseline','cpu':0,'threads':1},indent=2)+'\n')
records=[]
for strain,card in cards.items():
 expected=None
 for sequence,(variant,warmup) in enumerate([('baseline',True),('candidate',True),('baseline',False),('candidate',False),('candidate',False),('baseline',False)]):
  path=a.results/f'{strain}_{sequence}_{variant}.log';print('START',strain,variant,sequence,flush=True)
  start=time.perf_counter()
  with path.open('w') as log:subprocess.run(['taskset','-c','0',str(binaries[variant]),'-i',str(card)],env=env,cwd=a.root,stdout=log,stderr=subprocess.STDOUT,check=True)
  elapsed=time.perf_counter()-start
  values=dict(re.findall(r'^([^=\n]+)=([^\n]+)$',path.read_text(),re.M))
  assert values['completed']=='true' and values['mpi_ranks']=='1'
  assert values['accepted_steps' if strain=='finite' else 'load_steps_completed']=='20'
  if strain=='finite':assert values['rejected_steps']=='0'
  engineering={k:v for k,v in values.items() if k.startswith(('region.','contact.'))}
  if expected is None:expected=engineering
  else:assert engineering==expected
  records.append({'strain':strain,'variant':variant,'warmup':warmup,'sequence':sequence,'external_seconds':elapsed,'internal_seconds':float(values['total_seconds']),'nonlinear_iterations':int(values['nonlinear_iterations_total']),'log':path.name})
  (a.results/'timing.json').write_text(json.dumps(records,indent=2)+'\n');print('COMPLETE',elapsed,flush=True)
