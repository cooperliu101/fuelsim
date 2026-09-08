"""Diagnostic parent-face interpolants; these are not production quadrature rules."""
from pathlib import Path
import numpy as np
import argparse
import netCDF4 as nc
from analyze_h20_28_finite_transfer import extract
from diagnose_c3d20rt_contact_gap import ROOT,FACES,shape,project

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
local_to_row=[0,1,3,2,4,7,5,6]
parser=argparse.ArgumentParser()
parser.add_argument('--grid',type=int,choices=[8,16,32],nargs='+',default=[8,16])
args=parser.parse_args()
for grid in args.grid:
 stem='h20_28_finite_transfer_probe' if grid==8 else 'h20_28_finite_transfer_'+str(grid)+'_probe' 
 reference=extract(stem);points=reference[3];weights=reference[2][local_to_row]
 matrix=np.zeros((88,88));area=np.zeros(88)
 for ids in secondary:
  face=coordinates[ids];face_area=np.linalg.norm(np.cross(face[1]-face[0],face[3]-face[0]))
  samples=[]
  for x,y in points:
   n=shape(2*x-1,2*y-1)[0];p=n@face
   candidates=[(j,v) for j,f in enumerate(primary) if (v:=project(p,coordinates[f])) is not None]
   if not candidates:raise RuntimeError('No primary projection for interpolant node')
   j,(_,_,natural)=min(candidates,key=lambda v:abs(v[1][1]))
   basis=np.zeros(88);basis[ids]+=n;basis[primary[j]]-=shape(*natural)[0];samples.append(basis)
  rows=weights@np.array(samples)
  for local,node in enumerate(ids):
   a=face_area*(1 if local<4 else 5)/24;matrix[node]+=a*rows[local];area[node]+=a
 predicted=matrix[np.ix_(labels-1,columns-1)]/area[labels-1,None]
 difference=predicted-native
 print('diagnostic_grid',grid,'relative_frobenius',np.linalg.norm(difference)/np.linalg.norm(native),
       'maximum_absolute_coefficient_difference',np.max(abs(difference)),flush=True)
 np.savetxt(root/('h20_transfer_interpolant_'+str(grid)+'_prediction.csv'),predicted,delimiter=',',fmt='%.12e')
