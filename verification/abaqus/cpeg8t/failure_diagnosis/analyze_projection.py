"""Reconstruct projection intervals from native coordinates and displacements."""
import json
from pathlib import Path
import netCDF4
import numpy as np

ROOT = Path(__file__).resolve().parent
PRIMARY = [[3,4,7],[10,3,13]]
SECONDARY = [[14,15,18],[15,22,24]]


def coefficients(nodes):
    return np.array([nodes[2],.5*(nodes[1]-nodes[0]),.5*(nodes[0]+nodes[1])-nodes[2]])


def curve(nodes,x):
    c,b,a=coefficients(nodes)
    return c+b*x+a*x*x,b+2*a*x


def roots(coeff):
    coeff=np.polynomial.polynomial.polytrim(coeff,tol=1e-14*np.max(np.abs(coeff)))
    return sorted(float(x.real) for x in np.polynomial.polynomial.polyroots(coeff)
                  if abs(x.imag)<1e-9 and -1+1e-12<x.real<1-1e-12)


def cuts(s,p):
    result=[-1.,1.]
    for end in (-1,1):
        position,tangent=curve(p,end)
        coeff=coefficients(s).copy();coeff[0]-=position
        result+=roots(coeff@tangent)
    return sorted(result)


def project(p,position):
    c,b,a=coefficients(p);c=c-position
    polynomial=[c@b,b@b+2*a@c,3*a@b,2*a@a]
    trimmed=np.polynomial.polynomial.polytrim(polynomial,tol=1e-14*max(abs(v) for v in polynomial))
    solutions=[float(z.real) for z in np.polynomial.polynomial.polyroots(trimmed)
               if abs(z.imag)<1e-9 and -1-1e-12<=z.real<=1+1e-12 and
               np.polynomial.polynomial.polyval(z.real,[polynomial[1],2*polynomial[2],3*polynomial[3]])>0]
    if not solutions:
        return None
    return min(solutions,key=lambda z:np.linalg.norm(curve(p,z)[0]-position))


def analyze():
    with netCDF4.Dataset(ROOT.parent/'deforming_contact.e') as d:
        original=np.column_stack([d['coordx'][:],d['coordy'][:]])
    frames=json.loads((ROOT.parent/'deforming_contact_fields.json').read_text())['steps']['HISTORY']
    report=[]
    for frame in frames[1:]:
        u={v['nodeLabel']:v['data'] for v in frame['fields']['U']['values']}
        nodes=original+np.array([u[n][:2] for n in range(1,27)])
        primary=[nodes[np.array(ids)-1] for ids in PRIMARY]
        t0=curve(primary[0],-1)[1];t1=curve(primary[1],1)[1]
        angle=float(np.arctan2(np.linalg.det(np.array([t0,t1])),t0@t1))
        for edge,ids in enumerate(SECONDARY):
            secondary=nodes[np.array(ids)-1]
            intervals=[]
            for i,p in enumerate(primary):
                cc=cuts(secondary,p)
                for left,right in zip(cc[:-1],cc[1:]):
                    if project(p,curve(secondary,.5*(left+right))[0]) is not None:
                        intervals.append((left,right,i))
            intervals.sort()
            overlap=0.;uncovered=0.;physical_overlap=0.
            for before,after in zip(intervals[:-1],intervals[1:]):
                overlap=max(overlap,before[1]-after[0])
                uncovered=max(uncovered,after[0]-before[1])
                if before[1]>after[0]:
                    x,w=np.polynomial.legendre.leggauss(8)
                    physical_overlap+=.5*(before[1]-after[0])*sum(
                        wi*np.linalg.norm(curve(secondary,.5*((1-xi)*after[0]+(1+xi)*before[1]))[1])
                        for xi,wi in zip(x,w))
            # The shared endpoint has one closest projection on this secondary
            # curve, independent of which adjacent primary tangent is used.
            shared_coordinate=project(secondary,nodes[2])
            report.append(dict(time=frame['time'],secondary_edge=edge,intervals=intervals,
                               natural_overlap=overlap,natural_uncovered=uncovered,
                               physical_overlap_m=physical_overlap,primary_tangent_jump_radians=angle,
                               shared_endpoint_secondary_normal_coordinate=shared_coordinate))
    result=dict(rows=report,max_overlap_m=max(r['physical_overlap_m'] for r in report),
                affected_rows=sum(r['natural_overlap']>1e-10 or r['natural_uncovered']>1e-10 for r in report))
    (ROOT/'projection_intervals.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    r=analyze();print(json.dumps(dict(first=r['rows'][:2],max_overlap_m=r['max_overlap_m'],
                                     affected_rows=r['affected_rows']),indent=2))
