from pathlib import Path
import os,subprocess,shutil,csv,json,re,statistics
import argparse
parser=argparse.ArgumentParser(description="Alternating pinned single-thread Fuelsim timing; one warmup and five measured runs")
parser.add_argument("baseline", type=Path)
parser.add_argument("current", type=Path)
parser.add_argument("work", type=Path)
parser.add_argument("report", type=Path)
args=parser.parse_args()
repo=Path(__file__).resolve().parents[2]
base=args.work.resolve(); base.mkdir(parents=True, exist_ok=True)
binaries={"baseline":args.baseline.resolve(), "hardened":args.current.resolve()}
cards=['verification/thermal/dc3d20_transient.fsi','verification/fuelsim/transient_b73_rz_finite_creep.fsi','verification/fuelsim/transient_b527_engineering_contact.fsi']
cpu=min(os.sched_getaffinity(0)); env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
rows=[]
for card in cards:
 for variant in binaries:
  root=base/variant/Path(card).stem; local=root/card; local.parent.mkdir(parents=True,exist_ok=True)
  shutil.copy2(repo/card,local)
  for directory in ['meshes','moose','abaqus']:
   link=root/'verification'/directory
   if not link.exists(): link.symlink_to(repo/'verification'/directory,target_is_directory=True)
  if '/thermal/' in card:
   mesh_block=re.search(r'\[Mesh\](.*?)(?:\[\])',local.read_text(),re.S)[1]
   mesh_name=re.search(r'^\s*file\s*=\s*(\S+)',mesh_block,re.M)[1]
   link=local.parent/mesh_name
   if not link.exists(): link.symlink_to((repo/Path(card).parent/mesh_name).resolve())
 for repeat in range(6):
  for variant in (list(binaries) if repeat%2==0 else list(reversed(binaries))):
   local=base/variant/Path(card).stem/card; log=local.parent/f'run{repeat}.log'; timing=local.parent/f'time{repeat}.txt'
   with log.open('w') as out:
    result=subprocess.run(['/usr/bin/time','-f','%e %M','-o',str(timing),'taskset','-c',str(cpu),str(binaries[variant]),'-i',str(local)],stdout=out,stderr=subprocess.STDOUT,env=env)
   if result.returncode: raise RuntimeError(str(log))
   metrics={}
   for row in log.read_text().splitlines():
    if '=' in row:
     k,v=row.split('=',1); metrics[k.strip()]=v.strip()
   for summary in local.parent.glob('*summary.csv'):
    metrics.update(dict(list(csv.reader(summary.open()))[1:]))
   elapsed,rss=timing.read_text().split()
   rows.append({'case':Path(card).stem,'variant':variant,'repeat':repeat,'cpu':cpu,'wall_seconds':float(elapsed),'peak_rss_kib':int(rss), 'metrics':{k:v for k,v in metrics.items() if 'seconds' in k or 'workspace' in k or 'evaluations' in k or k=='accepted_steps'}})
args.report.write_text(json.dumps(rows,indent=2))
for card in cards:
 case=Path(card).stem
 for variant in binaries:
  selected=[r for r in rows if r['case']==case and r['variant']==variant and r['repeat']>0]
  print(case,variant,'wall median',statistics.median(r['wall_seconds'] for r in selected),'range',min(r['wall_seconds'] for r in selected),max(r['wall_seconds'] for r in selected),'rss median',statistics.median(r['peak_rss_kib'] for r in selected),flush=True)
