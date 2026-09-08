from pathlib import Path
import csv, hashlib, json, subprocess, sys
import numpy as np
from netCDF4 import Dataset, chartostring
root=Path('verification/abaqus/b13_small_larger_steps')
tag=sys.argv[1]; dt=float(sys.argv[2]); times=np.array([0.,.5]+[.5+dt*i for i in range(1,int(19.5/dt)+1)]+[20.]); count=len(times)-1; d=root/tag; summary={}
for el,nq in [('cax4t',4),('cax4rt',1),('cax8t',9),('cax8rt',4)]:
 c='b13_small_'+el; result=d/(el+'_results.e')
 meta=dict(x.split('=',1) for x in (d/(c+'_run.txt')).read_text().splitlines())
 assert hashlib.sha256((d/(c+'.inp')).read_bytes()).hexdigest()==meta['input_sha256'].lower()
 rows=list(csv.DictReader((d/(c+'_points.csv')).open()))
 assert sorted({float(x['time']) for x in rows})==list(times[1:])
 cmd=['build/fuelsim_production_pcmi_abaqus',str(result),str(d/(c+'_nodes.csv')),str(d/(c+'_points.csv')),str(d/(c+'_contact.csv')),'small',str(count)]
 run=subprocess.run(cmd,capture_output=True,text=True); (d/(el+'_metrics.txt')).write_text(run.stdout+run.stderr)
 vals=dict(x.split('=',1) for x in run.stdout.splitlines() if '=' in x)
 metrics={k:float(v) for k,v in vals.items() if k.endswith(('_relative_l2','_relative_absolute_peak','_maximum_pointwise_relative'))}
 with Dataset(result) as f:
  assert np.allclose(f['time_whole'][:],times,rtol=0,atol=1e-12)
  names=list(chartostring(f['name_elem_var'][:]))
  def field(name):return np.array(f[f'vals_elem_var{names.index(name)+1}eb2'][:])
  p=np.stack([field('equiv_plastic_q'+str(q)) for q in range(nq)],axis=2)
  cr=np.stack([field('equiv_creep_q'+str(q)) for q in range(nq)],axis=2)
  s=np.stack([np.stack([field('stress_'+k+'_q'+str(q)) for k in ['rr','zz','hoop','rz']],axis=-1) for q in range(nq)],axis=2)
 native=[x for x in rows if int(x['element'])>=25]
 ap=np.array([float(x['equiv_plastic']) for x in native]).reshape(count,10,nq)
 ac=np.array([float(x['equiv_creep']) for x in native]).reshape(count,10,nq)
 ass=np.array([[float(x['stress_'+k]) for k in ['rr','zz','hoop','rz']] for x in native]).reshape(count,10,nq,4)
 def audit(plastic,creep,stress,initial):
  dp=np.diff(plastic,axis=0) if initial else np.diff(plastic,axis=0,prepend=np.zeros((1,10,nq)))
  dc=np.diff(creep,axis=0) if initial else np.diff(creep,axis=0,prepend=np.zeros((1,10,nq)))
  st=stress[1:] if initial else stress
  q=np.sqrt(((st[...,0]-st[...,1])**2+(st[...,1]-st[...,2])**2+(st[...,2]-st[...,0])**2)/2+3*st[...,3]**2)
  be=np.diff(times)[:,None,None]*1e-5*(q/5e6)**3
  return {'simultaneous':int(((dp>1e-12)&(dc>1e-12)).sum()),'total':dp.size,'min_plastic_increment':float(dp.min()),'min_creep_increment':float(dc.min()),'max_backward_euler_relative_residual':float(np.max(abs(dc-be)/be))}
 a=audit(p,cr,s,True);b=audit(ap,ac,ass,False)
 summary[el]={'comparison_passed':run.returncode==0,'maximum_relative_percent':max(metrics.values())*100,'maximum_metric':max(metrics,key=metrics.get),'active_contact_samples':int(vals['active_contact_samples']),'fuelsim':a,'abaqus':b}
 print(el,summary[el],flush=True)
(d/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
assert all(x['comparison_passed'] and all(x[k]['simultaneous']==x[k]['total'] and x[k]['max_backward_euler_relative_residual']<1e-3 for k in ['fuelsim','abaqus']) for x in summary.values())
