"""Read production outputs and compare coupled material fields and activation history."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('result', type=Path)
p.add_argument('prefix', type=Path)
p.add_argument('--report', type=Path, required=True)
a = p.parse_args()
spec = importlib.util.spec_from_file_location('comparison', Path(__file__).resolve().parents[1] / 'compare.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
metrics = m.compare(a.result, a.prefix, 'cax4t', True, True)
with gzip.open(str(a.prefix)+'_points.csv.gz', 'rt') as f:
    points = list(csv.DictReader(f))
with gzip.open(str(a.prefix)+'_inelastic_history.csv.gz', 'rt') as f:
    history = list(csv.DictReader(f))
with m.result_database(a.result) as db:
    times = db.variables['time_whole']
    fields = {name: np.concatenate([db.variables['vals_elem_var%deb%d' % (i,b)] for b in range(1,db.dimensions['num_el_blk']+1)],axis=1)
              for i,name in enumerate(m.names(db.variables['name_elem_var']),1)}
qmap = [0,1,3,2]

def values(field, rows, frame=-1):
    return np.array([fields[field+'_q'+str(qmap[int(r['point'])-1])][frame,int(r['element'])-1] for r in rows])

for tensor in ('elastic','plastic','creep'):
    actual = np.stack([values(tensor+'_'+c,points) for c in ('rr','zz','hoop','rz')],axis=1)
    reference = np.array([[float(r[tensor+'_'+c]) for c in ('rr','zz','hoop','rz')] for r in points])
    actual[:,3] *= np.sqrt(2); reference[:,3] *= np.sqrt(2)
    metrics[tensor+'_strain_tensor'] = m.metric(actual,reference,1e-12)
for scalar in ('plastic','creep'):
    metrics['equivalent_'+scalar] = m.metric(values('equiv_'+scalar,points),[float(r['equiv_'+scalar]) for r in points],1e-12)

steps=[]
previous_a={};previous_f={}
if len(times) != 21 or not np.array_equal(times, np.arange(21)):
    raise ValueError('Expected the initial state and all 20 one-second increments')
for frame,time in enumerate(times):
    if frame==0:continue
    rows=[r for r in history if float(r['time'])==time]
    if not rows:raise ValueError('Missing material history frame')
    pa=np.array([float(r['equiv_plastic']) for r in rows]);ca=np.array([float(r['equiv_creep']) for r in rows]);qa=np.array([float(r['mises']) for r in rows])
    pf=values('equiv_plastic',rows,frame);cf=values('equiv_creep',rows,frame)
    rr,zz,hoop,rz=[values('stress_'+c,rows,frame) for c in ('rr','zz','hoop','rz')]
    qf=np.sqrt(.5*((rr-zz)**2+(zz-hoop)**2+(hoop-rr)**2)+3*rz**2)
    keys=[(int(r['element']),int(r['point'])) for r in rows]
    if len(keys) != 4096 or len(set(keys)) != 4096 or (previous_a and set(keys) != set(previous_a)):
        raise ValueError('Expected the same 4096 unique cladding material points in every increment')
    old_a=np.array([previous_a.get(k,(0,0,0)) for k in keys]);old_f=np.array([previous_f.get(k,(0,0,0)) for k in keys])
    dpa,dca=pa-old_a[:,0],ca-old_a[:,1];dpf,dcf=pf-old_f[:,0],cf-old_f[:,1]
    dt=float(time-times[frame-1]);threshold=1e-14
    def integration(delta,q,oldq):
        mask=delta>threshold
        if not mask.any():return {'active_samples':0}
        be=dt*1e-5*(q[mask]/5e6)**3;fe=dt*1e-5*(oldq[mask]/5e6)**3
        return {'active_samples':int(mask.sum()),'backward_euler_relative_l2':float(np.linalg.norm(delta[mask]-be)/np.linalg.norm(delta[mask])),
                'forward_euler_relative_l2':float(np.linalg.norm(delta[mask]-fe)/np.linalg.norm(delta[mask]))}
    steps.append({'time':float(time),'fuelsim_simultaneous_points':int(((dpf>threshold)&(dcf>threshold)).sum()),
                  'abaqus_simultaneous_points':int(((dpa>threshold)&(dca>threshold)).sum()),
                  'shared_simultaneous_points':int(((dpa>threshold)&(dca>threshold)&(dpf>threshold)&(dcf>threshold)).sum()),
                  'equivalent_plastic_comparison':m.metric(pf,pa,1e-12),'equivalent_creep_comparison':m.metric(cf,ca,1e-12),
                  'fuelsim_integration':integration(dcf,qf,old_f[:,2]),'abaqus_integration':integration(dca,qa,old_a[:,2])})
    for name, actual, reference in [('plastic',pf,pa),('creep',cf,ca)]:
        nonzero=np.flatnonzero(reference != 0)
        if nonzero.size:
            worst=nonzero[np.argmax(abs(actual[nonzero]-reference[nonzero])/abs(reference[nonzero]))]
            steps[-1]['equivalent_'+name+'_comparison']['worst_point']={
                'element':keys[worst][0], 'abaqus_point':keys[worst][1],
                'fuelsim':float(actual[worst]), 'abaqus':float(reference[worst]),
                'absolute_difference':float(abs(actual[worst]-reference[worst]))}
    previous_a=dict(zip(keys,zip(pa,ca,qa)));previous_f=dict(zip(keys,zip(pf,cf,qf)))
activation=any(r['shared_simultaneous_points']>0 for r in steps)
final_passed=all(v['passed'] for v in metrics.values())
history_project_passed=all(
    v['zero_reference_absolute_error'] <= v['zero_absolute_tolerance'] and
    all(v.get(k,0) < .005 for k in ('relative_l2','relative_absolute_peak','maximum_pointwise_relative'))
    for r in steps for v in (r['equivalent_plastic_comparison'],r['equivalent_creep_comparison']))
passed=activation and final_passed and history_project_passed
report={'all_final_metrics_passed':final_passed,'simultaneous_activation_verified':activation,'activation_increment_threshold':1e-14,
        'history_project_0_5_percent_passed':history_project_passed,
        'metrics':metrics,'material_history':steps,'history_metrics_all_passed':all(r['equivalent_plastic_comparison']['passed'] and r['equivalent_creep_comparison']['passed'] for r in steps)}
a.report.write_text(json.dumps(report,indent=2)+'\n')
for key,value in metrics.items():print(key,'PASS' if value['passed'] else 'FAIL',value.get('maximum_pointwise_relative',0)*100)
print('simultaneous_activation_verified',activation)
print('history_project_0_5_percent_passed',history_project_passed)
raise SystemExit(0 if passed else 1)
