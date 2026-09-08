"""Compare secondary normal-force rules at prescribed native geometry and gaps.

This isolates force transfer from the nonlinear solution and pressure recovery.
No candidate rule is fitted to the force values.
"""
import argparse
from pathlib import Path
import numpy as np
import netCDF4 as nc
from diagnose_c3d20rt_contact_gap import FACES, LOCATIONS, point, shape
from check_h20_disk_transfer import samples


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--frame',type=int,default=160)
    args=parser.parse_args(); root=Path(__file__).resolve().parents[2]
    stem=root/'verification/abaqus/c3d20rt_finite_contact_frictionless'
    nodes=np.loadtxt(str(stem)+'_nodes.csv',delimiter=',',skiprows=1).reshape(-1,88,10)[args.frame-1]
    contact=np.loadtxt(str(stem)+'_contact.csv',delimiter=',',skiprows=1).reshape(-1,13,14)[args.frame-1]
    with nc.Dataset(root/'verification/meshes/c3d20rt_finite_contact.e') as d:
        x=np.array([d['coord'+a][:] for a in 'xyz']).T+nodes[:,3:6]
        elements=np.vstack([d['connect1'][:],d['connect2'][:]])-1
        i=nc.chartostring(d['ss_names'][:]).tolist().index('secondary_contact')+1
        faces=[elements[int(e)-1,FACES[int(s)-1]] for e,s in zip(d['elem_ss'+str(i)][:],d['side_ss'+str(i)][:])]
    labels=contact[:,1].astype(int)-1; index={int(label):i for i,label in enumerate(labels)}
    coefficients={mode:np.zeros((13,13,3)) for mode in ('local','sampled','area_sampled','raw_area_sampled','common')}
    scalar=np.zeros((13,13)); normals=np.zeros((13,3)); areas=np.zeros(13)
    gauss,weights=np.polynomial.legendre.leggauss(3)
    for ids in faces:
        face=x[ids]; area=sum(wx*wy*point(face,a,b)[3] for a,wx in zip(gauss,weights) for b,wy in zip(gauss,weights))
        target=[index[int(i)] for i in ids]
        for i in range(8):
            row=target[i]; local_area=area*(1 if i<4 else 5)/24
            normal=point(face,*LOCATIONS[i])[2]
            averaging=np.zeros(8); sampled=np.zeros((8,3)); area_sampled=np.zeros((8,3)); weight_sum=0.
            for a,b,w in samples(i):
                n=shape(2*a-1,2*b-1)[0]; _,_,pn,measure=point(face,2*a-1,2*b-1)
                averaging+=w*n; sampled+=w*n[:,None]*pn
                area_sampled+=w*measure*n[:,None]*pn; weight_sum+=w*measure
            coefficients['local'][row,target]+=local_area*averaging[:,None]*normal
            coefficients['sampled'][row,target]+=local_area*sampled
            coefficients['area_sampled'][row,target]+=local_area*area_sampled/weight_sum
            coefficients['raw_area_sampled'][row,target]+=4*(1 if i<4 else 5)/24*area_sampled
            scalar[row,target]+=local_area*averaging; normals[row]+=local_area*normal; areas[row]+=local_area
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    coefficients['common']=scalar[:,:,None]*normals[:,None,:]
    pressure=-1e11*contact[:,2]; reference=contact[:,4:7]
    print('diagnostic_time',nodes[0,0])
    for mode,matrix in coefficients.items():
        predicted=-np.einsum('i,ijk->jk',pressure,matrix); error=predicted-reference
        print(mode,'relative_l2',np.linalg.norm(error)/np.linalg.norm(reference),
              'component_relative_l2',np.linalg.norm(error,axis=0)/np.linalg.norm(reference,axis=0),
              'component_maximum_absolute_n',np.max(abs(error),axis=0))
    frames=np.loadtxt(str(stem)+'_frames.csv',delimiter=',',skiprows=1).reshape(-1,13,10)[args.frame-1]
    if not np.array_equal(frames[:,1],contact[:,1]): raise RuntimeError('Native contact frame ordering mismatch')
    native_normal=np.cross(frames[:,2:5],frames[:,5:8])
    averaged_normal=coefficients['raw_area_sampled'].sum(axis=1)
    averaged_normal/=np.linalg.norm(averaged_normal,axis=1)[:,None]
    print('maximum_native_frame_normal_difference',
          np.max(np.minimum(np.linalg.norm(averaged_normal-native_normal,axis=1),
                            np.linalg.norm(averaged_normal+native_normal,axis=1))))


if __name__=='__main__': main()
