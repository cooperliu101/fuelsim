"""Output-only diagnosis of the four formal small-strain PCMI histories.

No problem generation, FE solve, or production-library access. Norton and J2
identities are evaluated independently from public material-point outputs.
"""
import csv
import json
from collections import Counter
from pathlib import Path

from analyze_b15_time import Exodus, np, ROOT
import hashlib

OUT = ROOT / 'verification/abaqus/small_inelastic_diagnosis'
OUT.mkdir(exist_ok=True)
REF = ROOT / 'verification/abaqus'
COMPONENTS = ['rr', 'zz', 'hoop', 'rz']
G = 7.5e10 / 2.6
A = 8e-26
H = 2e9


def dev(s):
    result = np.array(s, dtype=float)
    result[:3] -= np.mean(result[:3])
    return result


def norm(s):
    return float(np.sqrt(np.dot(s[:3], s[:3]) + 2*s[3]**2))


def mises(s):
    return np.sqrt(1.5) * norm(dev(s))


def root(trial, dt):
    # Independent monotone scalar equation; diagnostic, not a FE re-solve.
    lo, hi = 0., trial
    for _ in range(80):
        mid = .5*(lo+hi)
        if mid + 3*G*dt*A*mid**3 > trial:
            hi = mid
        else:
            lo = mid
    return .5*(lo+hi)


def read(element):
    case = 'b13_small_' + element
    run = dict(line.split('=',1) for line in (REF/(case+'_run.txt')).read_text().splitlines())
    assert hashlib.sha256((REF/(case+'.inp')).read_bytes()).hexdigest().upper() == run['input_sha256']
    assert hashlib.sha256((REF/'extract_b13_pcmi.py').read_bytes()).hexdigest().upper() == run['extractor_sha256']
    path = ROOT / ('build/blackbox/'+case+'/verification/fuelsim/transient_'+case+'_results.e')
    file = Exodus(path)
    times = file.get('time_whole')
    assert np.max(abs(times - np.arange(21))) < 1e-12
    names = file.names('name_elem_var')
    count = {'cax4t':4, 'cax4rt':1, 'cax8t':9, 'cax8rt':4}[element]
    order = list(range(9)) if count == 9 else [0, 1, 3, 2][:count]
    fields = [p+'_'+c for p in ['stress','elastic','plastic','creep'] for c in COMPONENTS]
    fields += ['equiv_creep','equiv_plastic']
    cache = {(field,q):file.get('vals_elem_var%deb2' % (names.index(field+'_q'+str(q))+1))
             for field in fields for q in order}
    file.close()
    actual, reference = {}, {}
    for t in range(21):
        for e in range(25,35):
            for p,q in enumerate(order,1):
                actual[t,e,p] = {field:float(cache[field,q][t,e-25]) for field in fields}
                if t == 0:
                    reference[t,e,p] = {field:0. for field in fields}
    with (REF/(case+'_points.csv')).open() as f:
        for r in csv.DictReader(f):
            t,e,p = round(float(r['time'])),int(r['element']),int(r['point'])
            if e >= 25:
                reference[t,e,p] = {field:float(r[field]) for field in fields}
    assert actual.keys() == reference.keys()
    assert all(abs(v) < 1e-30 for key,state in actual.items() if key[0] == 0 for v in state.values())
    return actual,reference,path


def tensor(r, prefix):
    return np.array([r[prefix+'_'+c] for c in COMPONENTS])


def diagnose(old,new,dt=1.):
    s0,s1 = tensor(old,'stress'),tensor(new,'stress')
    q0,q1 = mises(s0),mises(s1)
    dc = new['equiv_creep']-old['equiv_creep']
    dp = new['equiv_plastic']-old['equiv_plastic']
    forward,backward = dt*A*q0**3,dt*A*q1**3
    ef = abs(dc-forward)/max(abs(dc),forward,1e-300)
    eb = abs(dc-backward)/max(abs(dc),backward,1e-300)
    # Threshold is only for identifying meaningful activity, never acceptance.
    mode = ('inactive' if max(abs(dc),forward,backward) < 1e-12 else
            'both' if ef < 1e-7 and eb < 1e-7 else
            'explicit' if ef < 1e-7 else 'implicit' if eb < 1e-7 else 'other')
    dct = tensor(new,'creep')-tensor(old,'creep')
    dpt = tensor(new,'plastic')-tensor(old,'plastic')
    trial = 2*G*dev(tensor(new,'elastic')+dct+dpt)
    qc = root(mises(trial),dt)
    yold = 4e6+H*old['equiv_plastic']
    ynew = 4e6+H*new['equiv_plastic']
    direction0 = 1.5*dev(s0)/q0 if q0 else np.zeros(4)
    direction1 = 1.5*dev(s1)/q1 if q1 else np.zeros(4)
    return dict(q_old=q0,q_new=q1,creep_increment=dc,plastic_increment=dp,
                forward_prediction=forward,backward_prediction=backward,
                forward_relative_difference=ef,backward_relative_difference=eb,mode=mode,
                creep_forward_tensor_abs=norm(dct-forward*direction0),
                creep_backward_tensor_abs=norm(dct-backward*direction1),
                trial_q=mises(trial),creep_only_q=qc,old_yield=yold,
                implicit_activation_margin=qc-yold,end_yield_difference=q1-ynew,
                plastic_direction_abs=norm(dpt-dp*direction1))


def write(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    allrows, differences, summaries,provenance = [],[],{},[]
    focus_rows,plastic_mismatches=[] ,[]
    for element in ['cax4t','cax4rt','cax8t','cax8rt']:
        actual,reference,path = read(element)
        provenance.append(dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        case='b13_small_'+element
        for p in [REF/(case+s) for s in ['.inp','_points.csv','_run.txt']]+[
                ROOT/('verification/fuelsim/transient_'+case+'.fsi'),REF/'extract_b13_pcmi.py']:
            provenance.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        rows=[]
        for t,e,p in sorted(actual):
            if not t: continue
            oldkey=(t-1,e,p)
            for solver,states in [('fuelsim',actual),('abaqus',reference)]:
                r=dict(element_type=element,solver=solver,time=t,element=e,point=p,
                       **diagnose(states[oldkey],states[t,e,p]))
                rows.append(r)
            a,b=actual[t,e,p],reference[t,e,p]
            differences.append(dict(element_type=element,time=t,element=e,point=p,
                stress_difference_pa=norm(tensor(a,'stress')-tensor(b,'stress')),
                creep_difference=abs(a['equiv_creep']-b['equiv_creep']),
                plastic_difference=abs(a['equiv_plastic']-b['equiv_plastic']),
                fuelsim_creep=a['equiv_creep'],abaqus_creep=b['equiv_creep'],
                fuelsim_plastic=a['equiv_plastic'],abaqus_plastic=b['equiv_plastic']))
        allrows+=rows
        summary={}
        for solver in ['fuelsim','abaqus']:
            selected=[r for r in rows if r['solver']==solver]
            summary[solver]={'modes':dict(Counter(r['mode'] for r in selected))}
            active=[r for r in selected if abs(r['creep_increment'])>=1e-12]
            plastic=[r for r in selected if r['plastic_increment']>=1e-12]
            summary[solver]['first_creep']=active[0] if active else None
            summary[solver]['first_plastic']=plastic[0] if plastic else None
            summary[solver]['maximum_active_yield_difference_pa']=max(abs(r['end_yield_difference']) for r in plastic)
            summary[solver]['plastic_integration_modes']=dict(Counter(r['mode'] for r in plastic))
            matched=[r for r in selected if r['mode'] in ['explicit','implicit','both']]
            summary[solver]['maximum_matching_creep_tensor_absolute_difference']=max(
                r['creep_forward_tensor_abs'] if r['mode']=='explicit' else r['creep_backward_tensor_abs'] for r in matched)
            summary[solver]['maximum_plastic_direction_absolute_difference']=max(r['plastic_direction_abs'] for r in selected)
            summary[solver]['implicit_activation_disagreement_count']=int(sum(
                (r['plastic_increment']>=1e-12) != (r['implicit_activation_margin']>1.)
                for r in selected if abs(r['implicit_activation_margin'])>1.))
        ds=[r for r in differences if r['element_type']==element]
        for field,threshold in [('stress_difference_pa',1.),('creep_difference',1e-12),('plastic_difference',1e-12)]:
            summary['first_'+field]=next((r for r in ds if r[field]>threshold),None)
        significant=[r for r in ds if r['creep_difference']>1e-8]
        first_time=min(r['time'] for r in significant)
        focus=max((r for r in significant if r['time']==first_time),key=lambda r:r['creep_difference'])
        summary['first_creep_difference_above_1e_minus8']=focus
        focus_rows.extend(r for r in rows if all(r[k]==focus[k] for k in ['time','element','point']))
        bykey={(r['solver'],r['time'],r['element'],r['point']):r for r in rows}
        for a in rows:
            if a['solver']!='fuelsim':continue
            b=bykey['abaqus',a['time'],a['element'],a['point']]
            if (a['plastic_increment']>1e-12)!=(b['plastic_increment']>1e-12):
                plastic_mismatches.extend([a,b])
        summaries[element]=summary
    write('point_increments.csv',allrows)
    write('differences.csv',differences)
    write('first_significant_creep.csv',focus_rows)
    write('plastic_activation_differences.csv',plastic_mismatches)
    (OUT/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (OUT/'fuelsim_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    for element,s in summaries.items():
        print(element)
        for solver in ['fuelsim','abaqus']:
            v=s[solver]
            print(solver,v['modes'],'first creep',[(v['first_creep'] or {}).get(k) for k in ['time','element','point']],
                  'first plastic',[(v['first_plastic'] or {}).get(k) for k in ['time','element','point']],
                  'activation disagreement',v['implicit_activation_disagreement_count'])
        for f in ['stress_difference_pa','creep_difference','plastic_difference']:
            print('first',f,s['first_'+f])


if __name__=='__main__':main()
