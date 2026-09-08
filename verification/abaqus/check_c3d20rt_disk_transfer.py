"""Independent six-element flat-geometry check of disk-averaged extrapolation."""
from pathlib import Path
import numpy as np
import netCDF4 as nc
from diagnose_c3d20rt_contact_gap import ROOT, FACES, shape
from check_h20_disk_transfer import disk_rectangle, samples

root=Path(__file__).resolve().parent
native=np.loadtxt(root/'c3d20rt_finite_gap_matrix_native.csv',delimiter=',',skiprows=1)
steps=np.loadtxt(root/'c3d20rt_finite_gap_matrix_steps.csv',delimiter=',',skiprows=1)
columns=steps[1::2,1].astype(int)
labels=np.array([57,60,61,64,68,69,72,76,78,80,83,85,88])
with nc.Dataset(ROOT/'verification/meshes/c3d20rt_finite_contact.e') as d:
 coordinates=np.array([d['coord'+a][:] for a in 'xyz']).T
 elements=np.vstack([d['connect1'][:],d['connect2'][:]])-1
 names=nc.chartostring(d['ss_names'][:]).tolist()
 def faces(name):
  k=names.index(name)+1
  return [elements[int(e)-1,FACES[int(s)-1]] for e,s in zip(d['elem_ss'+str(k)][:],d['side_ss'+str(k)][:])]
 primary=faces('primary_contact');secondary=faces('secondary_contact')
r0=min(min(np.min(samples(i)[:,:2]),np.min(1-samples(i)[:,:2])) for i in range(8))
lengths=[np.array([np.linalg.norm(coordinates[ids[1]]-coordinates[ids[0]]),np.linalg.norm(coordinates[ids[3]]-coordinates[ids[0]])]) for ids in secondary]
for mode in ['local_minimum','global_minimum']:
 matrix=np.zeros((88,88));areas=np.zeros(88)
 for ids,length in zip(secondary,lengths):
  radius=r0*(min(length) if mode=='local_minimum' else min(min(l) for l in lengths) if mode=='global_minimum' else np.sqrt(np.prod(length)))
  face=coordinates[ids]
  for local in range(8):
   row=np.zeros(88)
   for x,y,w in samples(local):
    n=shape(2*x-1,2*y-1)[0];pos=n@face;row[ids]+=w*n
    total=0.
    for primary_ids in primary:
     q=coordinates[primary_ids];lo=q[:4,1:].min(axis=0);hi=q[:4,1:].max(axis=0)
     f=disk_rectangle(lo[0]-pos[1],hi[0]-pos[1],lo[1]-pos[2],hi[1]-pos[2],radius)
     if f<=0:continue
     total+=f
     natural=2*(pos[1:]-lo)/(hi-lo)-1
     row[primary_ids]-=w*f*shape(*natural)[0]
    if abs(total-1)>1e-9:raise ValueError('Incomplete disk support')
   a=np.prod(length)*(1 if local<4 else 5)/24;matrix[ids[local]]+=a*row;areas[ids[local]]+=a
 predicted=matrix[np.ix_(labels-1,columns-1)]/areas[labels-1,None]
 error=predicted-native
 print(mode,'relative',np.linalg.norm(error)/np.linalg.norm(native),'maximum',abs(error).max())
 np.savetxt(root/('c3d20rt_disk_'+mode+'_prediction.csv'),predicted,delimiter=',',fmt='%.16e')
