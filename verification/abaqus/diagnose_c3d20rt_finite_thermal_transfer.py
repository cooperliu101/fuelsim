"""Diagnostic replay at the first native state; never solves or changes an input.

Compare the existing single-face thermal transfer and the mechanically identified
disk transfer. The latter is a hypothesis for thermal contact, not a native
thermal-operator identification or an end-to-end acceptance result.
"""
import csv
import itertools
from pathlib import Path
import netCDF4 as nc
import numpy as np
from diagnose_c3d20rt_contact_gap import FACES, shape, point, project
from check_h20_disk_transfer import samples

ROOT = Path(__file__).resolve().parents[2]
SIGNS = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                  [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
EDGES = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),
         (4,5),(5,6),(6,7),(7,4)]


def body_shapes(q):
    factors = 1 + SIGNS*q
    n = np.prod(factors, axis=1)/8
    d = np.array([SIGNS[:,i]*np.prod(np.delete(factors,i,axis=1),axis=1)/8 for i in range(3)]).T
    d20 = np.empty((20,3))
    for i in range(8):
        d20[i] = d[i]*(SIGNS[i]@q-2)+n[i]*SIGNS[i]
    for e,(i,j) in enumerate(EDGES):
        s = (SIGNS[i]+SIGNS[j])/2
        axis = int(np.where(s == 0)[0][0]); other = [k for k in range(3) if k != axis]
        d20[8+e,axis] = -.5*q[axis]*np.prod(1+s[other]*q[other])
        for k in other:
            l = next(v for v in other if v != k)
            d20[8+e,k] = .25*(1-q[axis]**2)*s[k]*(1+s[l]*q[l])
    return n,d,d20


def circle_polygon(vertices):
    area = 0.
    for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
        edge=b-a; length=edge@edge
        if length == 0: continue
        linear=2*a@edge; disc=linear*linear-4*length*(a@a-1)
        cuts=[0.,1.]
        if disc > 0:
            cuts += [t for t in ((-linear-np.sqrt(disc))/(2*length),
                                 (-linear+np.sqrt(disc))/(2*length)) if 0<t<1]
        cuts.sort()
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            u=a+lo*edge; v=a+hi*edge; middle=(u+v)/2
            cross=u[0]*v[1]-u[1]*v[0]
            area += .5*(cross if middle@middle<1 else np.arctan2(cross,u@v))
    return abs(area)/np.pi


def extrapolate(position, face):
    natural=np.zeros(2)
    for _ in range(30):
        p,t,_,_=point(face,*natural)
        delta=np.linalg.solve(t.T@t,t.T@(position-p)); natural+=delta
        if np.linalg.norm(delta)<1e-13: return natural
    raise RuntimeError('Unconverged diagnostic projection')


def main():
    with nc.Dataset(ROOT/'verification/meshes/c3d20rt_finite_contact.e') as dataset:
        x=np.array([dataset['coord'+a][:] for a in 'xyz']).T
        elements=np.vstack([dataset['connect1'][:],dataset['connect2'][:]])-1
        names=nc.chartostring(dataset['ss_names'][:]).tolist()
        def faces(name):
            i=names.index(name)+1
            return [elements[int(e)-1,FACES[int(s)-1]] for e,s in
                    zip(dataset['elem_ss'+str(i)][:],dataset['side_ss'+str(i)][:])]
        primary=faces('primary_contact'); secondary=faces('secondary_contact')
    with (Path(__file__).parent/'c3d20rt_finite_contact_frictionless_nodes.csv').open() as f:
        rows=list(itertools.islice(csv.DictReader(f),88))
    u=np.array([[float(r['u'+a]) for a in 'xyz'] for r in rows]); current=x+u
    temperature=np.array([float(r['temperature']) for r in rows])
    old=np.where(np.arange(88)<56,300.,400.)
    native=np.array([float(r['heat_reaction']) for r in rows]); body=np.zeros(88)
    for ids in elements:
        for q in itertools.product([-1/np.sqrt(3),1/np.sqrt(3)],repeat=3):
            n,d,d20=body_shapes(np.array(q))
            measure=np.linalg.det(current[ids].T@d20)
            gradient=d@np.linalg.inv((x[ids]+u[ids]/2).T@d20)
            body[ids[:8]] += measure*(10*gradient@(gradient.T@temperature[ids[:8]]) +
                n*(n@(temperature[ids[:8]]-old[ids[:8]]))/.025)
    reference=native-body
    r0=.028276366456418445
    for mode in ('single_face','disk'):
        moments={}
        for ids in secondary:
            face=current[ids]
            radius=r0*min(np.linalg.norm(face[(i+1)%4]-face[i]) for i in range(4))
            for local in range(8):
                for a,b,w in samples(local):
                    location,tangents,_,measure=point(face,2*a-1,2*b-1)
                    ns=np.array([(1-a)*(1-b),a*(1-b),a*b,(1-a)*b])
                    transfer=np.zeros(88)
                    support=[]
                    if mode == 'single_face':
                        candidates=[(abs(q[1]),p,q[2]) for p in primary
                                    if (q:=project(location,current[p])) is not None]
                        _,p,natural=min(candidates,key=lambda entry:entry[0])
                        support=[(p,1.,natural)]
                    else:
                        first=tangents[:,0]/np.linalg.norm(tangents[:,0])
                        second=tangents[:,1]-first*(first@tangents[:,1]); second/=np.linalg.norm(second)
                        basis=np.array([first,second]).T
                        for p in primary:
                            polygon=((current[p]-location)@basis/radius)[[0,4,1,5,2,6,3,7]]
                            fraction=circle_polygon(polygon)
                            if fraction<1e-12: continue
                            support.append((p,fraction,extrapolate(location,current[p])))
                    for p,fraction,(xi,eta) in support:
                        transfer[p[:4]]+=fraction*.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),
                                                               (1+xi)*(1+eta),(1-xi)*(1+eta)])
                    transfer/=transfer.sum()
                    transfer[ids[:4]]-=ns
                    patches=[i for i in range(4) if
                             (a>=.5-1e-12 if i in (1,2) else a<=.5+1e-12) and
                             (b>=.5-1e-12 if i in (2,3) else b<=.5+1e-12)]
                    weight=4*measure*w*(1 if local<4 else 5)/24/len(patches)
                    for i in patches:
                        key=ids[i]
                        if key not in moments: moments[key]=[0.,np.zeros(88)]
                        moments[key][0]+=weight; moments[key][1]+=weight*transfer
        result=np.zeros(88)
        for area,coefficients in moments.values():
            result+=100*coefficients*(coefficients@temperature)/area
        error=result-reference
        print(mode,'relative_l2',np.linalg.norm(error)/np.linalg.norm(reference),
              'maximum_absolute_w',max(abs(error)),'net_heat_w',sum(result))


if __name__ == '__main__':
    main()
