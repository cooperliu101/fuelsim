"""Compare native contact matrices with independent surface integrals."""
import json
from pathlib import Path
import netCDF4
import numpy as np
from analyze_projection import curve,coefficients,project,cuts,PRIMARY,SECONDARY
from create_thermal_probe import CORNERS

ROOT=Path(__file__).resolve().parent


def averaged_matrix(nodes,thermal_dual=False):
    accum={}
    for sid in SECONDARY:
        s=nodes[np.array(sid)-1]
        for local,node in enumerate(sid[:2] if thermal_dual else sid):
            accum.setdefault(node,[0.,np.zeros(12)])
            for a,b,n in [(-1,-2/3,3),(-2/3,2/3,5),(2/3,1,3)]:
                for j in range(n):
                    xi=.5*(a+b)+(b-a)*(j-(n-1)/2)/np.sqrt(n*n-1)
                    test=([1 if xi<0 else .5 if xi==0 else 0,
                           1 if xi>0 else .5 if xi==0 else 0][local] if thermal_dual else
                          [max(-xi,0),max(xi,0),1-abs(xi)][local])
                    position,tangent=curve(s,xi);length=np.linalg.norm(tangent);tangent/=length
                    radius=(1-1/np.sqrt(2))/12*np.linalg.norm(s[1]-s[0])
                    for ids in PRIMARY:
                        p=nodes[np.array(ids)-1]
                        left=(p[1]-position)@tangent;right=(p[0]-position)@tangent
                        fraction=max(0,min(right,radius)-max(left,-radius))/(2*radius)
                        if fraction==0:continue
                        xp=((position-p[2])@tangent)/(.5*(p[1]-p[0])@tangent)
                        shape=np.zeros(12)
                        for i in (0,1):
                            shape[CORNERS.index(sid[i])]+=.5*(1+(-1 if i==0 else 1)*xi)
                            shape[CORNERS.index(ids[i])]-=.5*(1+(-1 if i==0 else 1)*xp)
                        area=.1*(b-a)/n*test*length*fraction
                        accum[node][0]+=area;accum[node][1]+=area*shape
    return sum(1000*np.outer(v,v)/area for area,v in accum.values())


def native_matrix(at):
    fields=[]
    for enabled in ('on','off'):
        n=json.loads((ROOT/f'thermal_matrix_{enabled}_fields.json').read_text())['steps']
        columns=[]
        for c in range(1,13):
            f=n[f'F{at}_C{c}'][-1]['fields']['RFL11']['values']
            v={r['nodeLabel']:r['data'] for r in f}
            columns.append([v[k] for k in CORNERS])
        fields.append(np.array(columns).T)
    return fields[0]-fields[1]


def matrix(nodes,rule='segmented',chord=False,smooth=False):
    result=np.zeros((12,12))
    for sid in SECONDARY:
        s=nodes[np.array(sid)-1].copy()
        if chord:s[2]=.5*(s[0]+s[1])
        intervals=[(-1.,1.)]
        if rule=='segmented':
            all_cuts=[-1.,1.]
            for ids in PRIMARY:all_cuts+=cuts(s,nodes[np.array(ids)-1])[1:-1]
            all_cuts=sorted(set(all_cuts));intervals=list(zip(all_cuts[:-1],all_cuts[1:]))
        samples=[]
        if rule=='equal11':
            for left,right,count in [(-1,-2/3,3),(-2/3,2/3,5),(2/3,1,3)]:
                samples += [(.5*(left+right)+(right-left)*(j-(count-1)/2)/np.sqrt(count*count-1),
                             (right-left)/count) for j in range(count)]
        else:
            x,w=np.polynomial.legendre.leggauss(3 if rule in ('segmented','gauss3') else 64)
            for a,b in intervals:samples += [(.5*((1-xi)*a+(1+xi)*b),.5*(b-a)*wi) for xi,wi in zip(x,w)]
        for xi,weight in samples:
            position,tangent=curve(s,xi);length=np.linalg.norm(tangent);tangent=tangent/length
            radius=(1-1/np.sqrt(2))/12*np.linalg.norm(s[1]-s[0])
            for ids in PRIMARY:
                p=nodes[np.array(ids)-1]
                xp=project(p,position)
                fraction=1.
                if smooth:
                    left=(p[1]-position)@tangent;right=(p[0]-position)@tangent
                    fraction=max(0,min(right,radius)-max(left,-radius))/(2*radius)
                    # Flat primary: normal-to-secondary projection can be
                    # evaluated even slightly beyond an element endpoint.
                    xp=((position-p[2])@tangent)/(.5*(p[1]-p[0])@tangent)
                if xp is None or fraction==0:continue
                shape=np.zeros(12)
                for i in (0,1):
                    shape[CORNERS.index(sid[i])]+=.5*(1+(-1 if i==0 else 1)*xi)
                    shape[CORNERS.index(ids[i])]-=.5*(1+(-1 if i==0 else 1)*xp)
                result+=1000*.1*weight*length*fraction*np.outer(shape,shape)
    return result


def analyze():
    with netCDF4.Dataset(ROOT.parent/'deforming_contact.e') as d:
        original=np.column_stack([d['coordx'][:],d['coordy'][:]])
    frames=json.loads((ROOT.parent/'deforming_fixed_primary_fields.json').read_text())['steps']['HISTORY']
    report={}
    for at in (1,13):
        u={v['nodeLabel']:v['data'] for v in frames[at]['fields']['U']['values']}
        nodes=original+np.array([u[n][:2] for n in range(1,27)])
        native=native_matrix(at);metrics={}
        m=averaged_matrix(nodes)
        metrics['mechanical_neighborhood_average']=dict(
            relative_frobenius=float(np.linalg.norm(m-native)/np.linalg.norm(native)),
            maximum_absolute=float(np.max(abs(m-native))))
        m=averaged_matrix(nodes,True)
        metrics['thermal_corner_dual_average']=dict(
            relative_frobenius=float(np.linalg.norm(m-native)/np.linalg.norm(native)),
            maximum_absolute=float(np.max(abs(m-native))))
        for rule in ['segmented','gauss3','gauss64','equal11']:
            for chord in (False,True):
                for smooth in (False,True):
                    m=matrix(nodes,rule,chord,smooth)
                    metrics[f'{rule}/chord={chord}/smooth={smooth}']=dict(
                        relative_frobenius=float(np.linalg.norm(m-native)/np.linalg.norm(native)),
                        maximum_absolute=float(np.max(abs(m-native))))
        report[str(at)]=dict(time=frames[at]['time'],native_matrix=native.tolist(),metrics=metrics,
                            native_singular_values=np.linalg.svd(native,compute_uv=False).tolist(),
                            native_symmetry=float(np.max(abs(native-native.T))),
                            native_conservation=float(np.max(abs(native.sum(axis=0)))))
    (ROOT/'thermal_matrix_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    r=analyze()
    for at,row in r.items():
        print(at,sorted(row['metrics'].items(),key=lambda v:v[1]['relative_frobenius'])[:5])
