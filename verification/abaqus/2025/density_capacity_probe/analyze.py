"""Identify native density and capacity from prescribed fields and independent quadrature."""
import argparse
import collections
import csv
import itertools
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
HEX=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
QUAD=HEX[:4,:2]

def geometry(x):
    signs=HEX if x.shape[1]==3 else QUAD
    values=[];weights=[]
    for xi in signs/np.sqrt(3):
        factors=1+signs*xi
        n=np.prod(factors,axis=1)/len(signs)
        dn=np.column_stack([signs[:,d]*np.prod(factors[:,[k for k in range(x.shape[1]) if k!=d]],axis=1)/len(signs) for d in range(x.shape[1])])
        w=np.linalg.det(x.T@dn)
        if x.shape[1]==2:w*=2*np.pi*(n@x[:,0])
        assert w>0
        values.append(n);weights.append(w)
    return np.array(values),np.array(weights)

def average_gradient(x):
    signs=HEX if x.shape[1]==3 else QUAD
    gradients=[]
    for xi in signs/np.sqrt(3):
        factors=1+signs*xi
        dn=np.column_stack([signs[:,d]*np.prod(factors[:,[k for k in range(x.shape[1]) if k!=d]],axis=1)/len(signs) for d in range(x.shape[1])])
        gradients.append(dn@np.linalg.inv(x.T@dn))
    _,w=geometry(x)
    return sum(g*weight for g,weight in zip(gradients,w))/sum(w)

def rho(t):return np.interp(t,[300,360,420,480],[2000,1900,1400,1300])
def write(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def validate_input(layout, native):
    coordinates={}; connectivity={}; initial={}; section=''
    for line in (HERE/'density_capacity.inp').read_text().splitlines():
        if line.startswith('**') or not line.strip():
            continue
        if line.startswith('*'):
            section=line.split(',')[0].lower()
            continue
        fields=[value.strip() for value in line.split(',')]
        if section=='*node':coordinates[int(fields[0])]=[float(v) for v in fields[1:]]
        elif section=='*element':connectivity[int(fields[0])]=[int(v) for v in fields[1:]]
        elif section=='*initial conditions':initial[int(fields[0])]=float(fields[1])
    assert len(coordinates)==960 and len(connectivity)==144 and len(initial)==960
    assert len(native)==4800
    increment=layout.get('temperature_increment',1)
    for n in layout['nodes']:
        assert coordinates[n['id']]==n['x'] and initial[n['id']]==n['t0']
        for k,stage in enumerate(['REFERENCE','DEFORM_A','HOLD_A','DEFORM_B','RETURN']):
            expected=n['t0']+((k+1)*increment if n['active'] else 0)
            assert abs(float(native[stage,n['id']]['temperature'])-expected)<1e-10
    for e in layout['elements']:
        assert connectivity[e['id']]==e['nodes']

def main():
    layout=json.loads((HERE/'probe_layout.json').read_text());rate=layout.get('temperature_increment',1)/layout.get('time_increment',1);ns={n['id']:n for n in layout['nodes']}
    native={(r['step'],int(r['node'])):r for r in csv.DictReader((HERE/'nodes.csv').open())}
    validate_input(layout,native)
    rows=[];base=[];identified=[];inactive=0
    for e in layout['elements']:
        ids=e['nodes'];i=e['active'];dim=len(ns[ids[0]]['x']);n=len(ids)
        x=np.array([ns[k]['x'] for k in ids]);t0=np.array([ns[k]['t0'] for k in ids]);sh,w0=geometry(x);v0=sum(w0);m0=sh.T@w0
        rhos={'mean_temperature':rho(t0.mean()),'volume_mean_temperature':rho((sh@t0)@w0/v0),'mean_density':rho(t0).mean(),'volume_mean_density':(sh@rho(t0))@w0/v0,'current_node':rho(float(native['REFERENCE',ids[i]]['temperature'])),'paired_node':rho(t0[i]),'paired_gauss':rho((sh@t0)[i])}
        if e['material']=='CONSTANT':rhos={k:2000. for k in rhos}
        for stage in ['REFERENCE','DEFORM_A','HOLD_A','DEFORM_B','RETURN']:
            obs=[native[stage,k] for k in ids];current=x+np.array([[float(r['u'+str(d+1)]) for d in range(dim)] for r in obs]);_,wc=geometry(current);vc=sum(wc);mc=sh.T@wc
            cp=3000 if e['material']!='DENSITY_CP' else 3000+4*(float(obs[i]['temperature'])-300)
            actual=float(obs[i]['rfl']);inactive=max(inactive,max(abs(float(obs[j]['rfl'])) for j in range(n) if j!=i))
            if stage=='REFERENCE':
                for rname,r in rhos.items():
                    for wname,w in [('paired_gauss',w0[i]),('row_sum',m0[i])]:
                        pred=r*cp*w*rate;base.append(dict(type=e['type'],geometry=e['geometry'],material=e['material'],element=e['id'],node=i+1,density_rule=rname,capacity_rule=wname,predicted=pred,actual=actual,relative_error=abs(pred/actual-1)))
            # Divide by the observed undeformed capacity to isolate geometry, independent of density mapping.
            ref=native['REFERENCE',ids[i]];cpref=3000 if e['material']!='DENSITY_CP' else 3000+4*(float(ref['temperature'])-300)
            measured_ratio=actual/float(ref['rfl'])*cpref/cp
            center_gradient=(HEX if dim==3 else QUAD)/n
            center_ratio=np.linalg.det(current.T@center_gradient)/np.linalg.det(x.T@center_gradient)
            if dim==2:center_ratio*=current[:,0].mean()/x[:,0].mean()
            mean_f_ratio=np.linalg.det(current.T@average_gradient(x))
            if dim==2:mean_f_ratio*=sum(sh@current[:,0]*w0)/sum(sh@x[:,0]*w0)
            rules={'rows_over_mean_f':mc[i]/m0[i]/mean_f_ratio,'rows_over_center':mc[i]/m0[i]/center_ratio,'fixed_initial_mass':1.,'global_volume_ratio':vc/v0,'local_volume_ratio':wc[i]/w0[i], 'global_over_local':vc/v0*w0[i]/wc[i], 'current_rows_over_initial_rows':mc[i]/m0[i], 'rows_over_local':mc[i]/m0[i]*w0[i]/wc[i], 'rows_over_global':mc[i]/m0[i]*v0/vc, 'integrated_density':sum(sh[:,i]*wc*w0/wc)/m0[i]}
            if e['type']=='C3D8T':
                capacity=rhos['mean_temperature']*cp*w0[i]*(vc/v0)/(wc[i]/w0[i])
            elif e['type']=='CAX4T':
                capacity=rhos['mean_temperature']*cp*mc[i]/(vc/v0)
            else:
                current_rho=2000 if e['material']=='CONSTANT' else rho(float(obs[i]['temperature']))
                jacobian=mean_f_ratio if e['type']=='C3D8RT' else vc/v0
                capacity=current_rho*cp*mc[i]/jacobian
            prediction=capacity*rate
            identified.append(dict(type=e['type'],geometry=e['geometry'],material=e['material'],step=stage,
                                   element=e['id'],node=i+1,prediction_w=prediction,native_rfl_w=actual,
                                   relative_error=abs(prediction/actual-1),absolute_error_w=abs(prediction-actual)))
            for rule,pred in list(rules.items()):
                if e['material']!='CONSTANT':rules[rule+'_current_density']=pred*rho(float(obs[i]['temperature']))/rho(float(ref['temperature']))
                else:rules[rule+'_current_density']=pred
            for rule,pred in rules.items():rows.append(dict(type=e['type'],geometry=e['geometry'],material=e['material'],step=stage,element=e['id'],node=i+1,rule=rule,prediction=pred,measured_ratio=measured_ratio,relative_error=abs(pred/measured_ratio-1)))
    write('initial_density_candidates.tsv',base);write('geometry_candidates.tsv',rows)
    summary=[]
    for typ,mat in itertools.product(['C3D8T','C3D8RT','CAX4T','CAX4RT'],['CONSTANT','DENSITY','DENSITY_CP']):
        selected=[r for r in base if r['type']==typ and r['material']==mat]
        ranks=[]
        for d,w in itertools.product(rhos,['paired_gauss','row_sum']):
            err=max(r['relative_error'] for r in selected if r['density_rule']==d and r['capacity_rule']==w);ranks.append((err,d,w))
        print('INITIAL',typ,mat,sorted(ranks)[:2])
    for typ in ['C3D8T','C3D8RT','CAX4T','CAX4RT']:
        for rule in rules:
            selected=[r for r in rows if r['type']==typ and r['rule']==rule]
            summary.append(dict(type=typ,rule=rule,maximum_relative_error=max(r['relative_error'] for r in selected)))
        print('GEOMETRY',typ,sorted([(r['maximum_relative_error'],r['rule']) for r in summary if r['type']==typ]))
    write('geometry_summary.tsv',summary)
    write('identified_predictions.tsv',identified)
    direct=[]
    for typ in ['C3D8T','C3D8RT','CAX4T','CAX4RT']:
        selected=[r for r in identified if r['type']==typ]
        direct.append(dict(type=typ,samples=len(selected),maximum_relative_error=max(r['relative_error'] for r in selected),
                           maximum_absolute_error_w=max(r['absolute_error_w'] for r in selected)))
    write('identified_summary.tsv',direct)
    mass_rows=[]
    element_by_id={e['id']:e for e in layout['elements']}
    grouped=collections.defaultdict(list)
    for row in identified:
        if row['material']=='CONSTANT':grouped[row['type'],row['geometry'],row['step']].append(row)
    for (typ,kind,stage),samples in grouped.items():
        e=element_by_id[samples[0]['element']]
        initial_volume=sum(geometry(np.array([ns[k]['x'] for k in e['nodes']]))[1])
        effective_mass=sum(r['native_rfl_w'] for r in samples)/(3000*rate)
        mass_rows.append(dict(type=typ,geometry=kind,step=stage,initial_mass_kg=2000*initial_volume,
                              inferred_capacity_mass_kg=effective_mass,
                              relative_mass_change=effective_mass/(2000*initial_volume)-1))
    write('capacity_mass.tsv',mass_rows)
    print('IDENTIFIED',direct)
    print('MAX_INACTIVE_RFL',inactive)
    assert len(identified)==720
    assert max(r['relative_error'] for r in identified)<1e-10
    assert inactive<1e-12
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=HERE)
    args=parser.parse_args()
    HERE=args.directory.resolve()
    main()
