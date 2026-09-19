from pathlib import Path
import subprocess,shutil,re,os,json,numpy as np
from netCDF4 import Dataset
import argparse
parser=argparse.ArgumentParser(description="Cross-read existing full input-card checkpoints with two Fuelsim executables")
parser.add_argument("baseline", type=Path)
parser.add_argument("current", type=Path)
parser.add_argument("work", type=Path)
parser.add_argument("report", type=Path)
args=parser.parse_args()
repo=Path(__file__).resolve().parents[2]
root=args.work.resolve()
bin={"old":args.baseline.resolve(), "new":args.current.resolve()}
pairs=[('thermal','verification/thermal','study_adaptive_first.fsi','study_adaptive_resume.fsi'),('radial','verification/fuelsim','transient_gps_contact_checkpoint.fsi','transient_gps_contact_restart.fsi'),('rz4','verification/fuelsim','transient_large_sliding_contact_checkpoint.fsi','transient_large_sliding_contact_restart.fsi'),('rz8','verification/fuelsim','transient_b125_cax8rt_finite_coupled_split.fsi','transient_b125_cax8rt_finite_coupled_restart.fsi'),('plane','verification/abaqus/cpeg8t','coupled_contact_split.fsi','coupled_contact_restart.fsi'),('hex8','verification/fuelsim','transient_hex8_sts_friction_checkpoint.fsi','transient_hex8_sts_friction_restart.fsi'),('hex20','verification/fuelsim','transient_hex20_coupled_checkpoint.fsi','transient_hex20_coupled_restart.fsi')]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
records=[]
for name,directory,split,restart in pairs:
 outputs={}
 for writer,reader in [('old','old'),('old','new'),('new','old'),('new','new')]:
  stage=root/name/(writer+'-'+reader)/directory; stage.mkdir(parents=True,exist_ok=True)
  for card in [split,restart]:
   source=repo/directory/card; target=stage/card; shutil.copy2(source,target)
   match=re.search(r'\[Mesh\](.*?)(?:\[\])',source.read_text(),re.S)
   mesh=re.search(r'^\s*file\s*=\s*(\S+)',match[1],re.M)[1]
   dest=(stage/mesh).resolve(); dest.parent.mkdir(parents=True,exist_ok=True)
   if not dest.exists(): dest.symlink_to((source.parent/mesh).resolve())
  for variant,card in [(writer,split),(reader,restart)]:
   with (stage/(card+'.log')).open('w') as log:
    run=subprocess.run([str(bin[variant]),'-i',str(stage/card)],stdout=log,stderr=subprocess.STDOUT,env=env)
   if run.returncode: raise RuntimeError(str(stage/card))
  expected={}
  for result in stage.glob('*'):
   if result.suffix=='.e' and not result.is_symlink():
    with Dataset(result) as ds:
     for field in ds.variables:
      values=np.asarray(ds[field][:]);
      if values.dtype.kind in 'fiu': expected[result.name+'/'+field]=values
   elif result.suffix in ['.checkpoint','.chk']: expected[result.name]=result.read_bytes()
   elif result.suffix=='.csv' and 'history' in result.name: expected[result.name]=result.read_bytes()
  outputs[writer+'-'+reader]=expected
 reference=outputs['old-old']; differences=[]; values=0
 for variant,results in outputs.items():
  if results.keys()!=reference.keys(): differences.append([variant,'output names differ']); continue
  for key,x in reference.items():
   y=results[key]
   if isinstance(x,np.ndarray):
    values+=x.size
    if x.shape!=y.shape or not np.array_equal(x,y,equal_nan=True): differences.append([variant,key])
   elif x!=y: differences.append([variant,key])
 records.append({'case':name,'runs':8,'compared_numeric_values':values,'differences':differences})
 print(json.dumps(records[-1]),flush=True)
args.report.write_text(json.dumps(records,indent=2))
if any(row['differences'] for row in records): raise RuntimeError('Compatibility results differ')
