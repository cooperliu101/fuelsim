"""Diagnostic only: compare geometric averaging variants at native nodal states.

Reads the existing identified parent-face transfer rule from C++ source. This is
not independent validation of that rule and does not alter production inputs.
"""
from pathlib import Path
import argparse
import re
import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FACES = np.array([[0,1,5,4,8,13,16,12], [1,2,6,5,9,14,17,13],
                  [2,3,7,6,10,15,18,14], [3,0,4,7,11,12,19,15],
                  [0,3,2,1,11,10,9,8], [4,5,6,7,16,17,18,19]])
LOCATIONS = [(-.75,-.75),(.75,-.75),(.75,.75),(-.75,.75),(0,-.5),(.5,0),(0,.5),(-.5,0)]

def shape(x, y):
    signs = np.array([[-1,-1],[1,-1],[1,1],[-1,1]])
    sx, sy = signs.T
    n = np.empty(8); dx = np.empty(8); dy = np.empty(8)
    n[:4] = .25*(1+sx*x)*(1+sy*y)*(sx*x+sy*y-1)
    dx[:4] = .25*sx*(1+sy*y)*(2*sx*x+sy*y)
    dy[:4] = .25*sy*(1+sx*x)*(sx*x+2*sy*y)
    n[4:] = [.5*(1-x*x)*(1-y),.5*(1+x)*(1-y*y),.5*(1-x*x)*(1+y),.5*(1-x)*(1-y*y)]
    dx[4:] = [-x*(1-y),.5*(1-y*y),-x*(1+y),-.5*(1-y*y)]
    dy[4:] = [-.5*(1-x*x),-(1+x)*y,.5*(1-x*x),-(1-x)*y]
    return n, dx, dy

def point(face, x, y):
    n, dx, dy = shape(x,y)
    tangents = np.array([dx@face,dy@face]).T
    normal = np.cross(tangents[:,0],tangents[:,1]); measure = np.linalg.norm(normal)
    return n@face, tangents, normal/measure, measure

def project(position, face):
    natural = np.zeros(2)
    for _ in range(25):
        p, tangent, normal, _ = point(face,*natural)
        delta = np.linalg.solve(tangent.T@tangent,tangent.T@(position-p))
        natural += delta
        if np.linalg.norm(delta) < 1e-13: break
    else: raise RuntimeError('Projection did not converge')
    if np.max(abs(natural)) > 1+1e-8: return None
    p, _, normal, _ = point(face,*natural)
    return p, float(normal@(position-p)), natural

def rules():
    source = (ROOT/'elements/src/quad8_face.cpp').read_text()
    result = []
    for name in ['corner','edge']:
        body = source.split('abaqus_quad8_'+name+'_transfer_rule() {',1)[1].split('return value;',1)[0]
        values = np.array([[float(v) for v in match] for match in re.findall(r'\{\s*([\d.e+-]+),\s*([\d.e+-]+),\s*([\d.e+-]+)\s*\}',body)])
        assert len(values) == (12 if name == 'corner' else 40)
        result.append(values)
    return result

def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--frame',type=int,default=1)
    parser.add_argument('--reference-prefix', default='c3d20rt_finite_contact_friction')
    parser.add_argument('--disk-radius', choices=['current','reference'])
    parser.add_argument('--disk-plane', choices=['point','constraint'], default='point')
    args = parser.parse_args()
    nodes = np.loadtxt(ROOT/'verification/abaqus'/(args.reference_prefix+'_nodes.csv'),delimiter=',',skiprows=1).reshape(-1,88,10)[args.frame-1]
    contact = np.loadtxt(ROOT/'verification/abaqus'/(args.reference_prefix+'_contact.csv'),delimiter=',',skiprows=1).reshape(-1,13,14)[args.frame-1]
    with nc.Dataset(ROOT/'verification/meshes/c3d20rt_finite_contact.e') as d:
        reference_coordinates = np.array([d['coord'+a][:] for a in 'xyz']).T
        coordinates = reference_coordinates + nodes[:,3:6]
        elements = np.vstack([d['connect1'][:],d['connect2'][:]])-1
        names = nc.chartostring(d['ss_names'][:]).tolist()
        def faces(name):
            index = names.index(name)+1
            return [elements[int(e)-1,FACES[int(s)-1]] for e,s in zip(d['elem_ss'+str(index)][:],d['side_ss'+str(index)][:])]
        primary_ids = faces('primary_contact')
        primary = [coordinates[f] for f in primary_ids]
        primary_reference = [reference_coordinates[f] for f in primary_ids]
        secondary = faces('secondary_contact')
    samples = rules(); accum = {}; area_sum = {}
    for ids in secondary:
        face = coordinates[ids]
        gauss, weights = np.polynomial.legendre.leggauss(3)
        area = sum(wx*wy*point(face,x,y)[3] for x,wx in zip(gauss,weights) for y,wy in zip(gauss,weights))
        for local in range(8):
            transfer = samples[0 if local<4 else 1].copy()
            for v in transfer:
                x,y = v[:2]
                if local==1: v[0]=1-x
                elif local==2: v[:2]=[1-x,1-y]
                elif local==3: v[1]=1-y
                elif local==5: v[:2]=[1-y,x]
                elif local==6: v[1]=1-y
                elif local==7: v[:2]=[y,x]
            transfer[:,2] /= transfer[:,2].sum()
            normal = point(face,*LOCATIONS[local])[2]
            values = np.zeros(7); measure_sum = 0.
            for x,y,w in transfer:
                p,_,pn,j = point(face,2*x-1,2*y-1)
                candidates = [v for f in primary if (v:=project(p,f)) is not None]
                if not candidates: raise RuntimeError('Transfer point lost all primary projections')
                projected, pgap, _ = min(candidates,key=lambda v:abs(v[1]))
                if args.disk_radius:
                    from diagnose_c3d20rt_finite_thermal_transfer import circle_polygon, extrapolate
                    radius_face=face if args.disk_radius=='current' else reference_coordinates[ids]
                    radius=.028276366456418445*min(np.linalg.norm(radius_face[(i+1)%4]-radius_face[i]) for i in range(4))
                    _,basis,_,_=point(face,*(LOCATIONS[local] if args.disk_plane=='constraint' else (2*x-1,2*y-1)))
                    first=basis[:,0]/np.linalg.norm(basis[:,0])
                    second=basis[:,1]-first*(first@basis[:,1]); second/=np.linalg.norm(second)
                    basis=np.array([first,second]).T
                    total=0.; projected=np.zeros(3)
                    for f in primary:
                        polygon=((f-p)@basis/radius)[[0,4,1,5,2,6,3,7]]
                        fraction=circle_polygon(polygon)
                        if fraction<1e-12: continue
                        projected+=fraction*point(f,*extrapolate(p,f))[0];total+=fraction
                    projected/=total
                separation = projected-p
                refp,_,refnormal,_ = point(reference_coordinates[ids],2*x-1,2*y-1)
                ref_candidates = [(i,v) for i,f in enumerate(primary_reference) if (v:=project(refp,f)) is not None]
                primary_index, (_,_,natural) = min(ref_candidates,key=lambda v:abs(v[1][1]))
                fixed_projected = point(primary[primary_index],*natural)[0]
                fixed_separation = fixed_projected-p
                values += w*np.array([normal@separation,pn@separation,j*(normal@separation),j*(pn@separation),pgap,
                                      normal@fixed_separation, refnormal@fixed_separation])
                measure_sum += w*j
            values[2:4] /= measure_sum
            factor = area*(1 if local<4 else 5)/24
            node = int(ids[local])+1
            accum[node] = accum.get(node,np.zeros(7)) + factor*values
            area_sum[node] = area_sum.get(node,0.) + factor
    print('diagnostic_only,time='+str(nodes[0,0]))
    print('node,abaqus_gap,fixed_normal,point_normal,area_fixed_normal,area_point_normal,primary_normal,fixed_projection_current_normal,fixed_projection_reference_normal')
    for row in contact:
        node = int(row[1]); values = accum[node]/area_sum[node]
        print(','.join([str(node),format(row[2],'.12e')]+[format(v,'.12e') for v in values]))

if __name__ == '__main__': main()
