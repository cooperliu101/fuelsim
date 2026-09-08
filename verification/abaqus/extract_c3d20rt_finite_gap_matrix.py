"""Extract the prescribed-displacement native gap matrix; no production solve."""
from pathlib import Path
import netCDF4 as nc
import numpy as np
from diagnose_c3d20rt_contact_gap import ROOT,FACES,rules,project,shape

r=ROOT/'verification/abaqus'
steps=np.loadtxt(r/'c3d20rt_finite_gap_matrix_steps.csv',delimiter=',',skiprows=1)
contact=np.loadtxt(r/'c3d20rt_finite_gap_matrix_contact.csv',delimiter=',',skiprows=1).reshape(-1,13,14)
assert len(steps)==len(contact) and np.max(abs(contact[:,:,0]-np.arange(1,len(steps)+1)[:,None]*.025))<1e-6
labels=contact[0,:,1].astype(int)
assert np.all(contact[:,:,1]==labels) and np.max(abs(contact[0,:,2]+1e-5))<1e-12
columns=steps[1::2,1].astype(int)
assert np.array_equal(columns,steps[2::2,1].astype(int)) and np.all(steps[1::2,2]==-steps[2::2,2])
native=((contact[1::2,:,2]-contact[2::2,:,2])/(2*steps[1::2,2,None])).T
with nc.Dataset(ROOT/'verification/meshes/c3d20rt_finite_contact.e') as d:
 coordinates=np.array([d['coord'+a][:] for a in 'xyz']).T
 elements=np.vstack([d['connect1'][:],d['connect2'][:]])-1
 names=nc.chartostring(d['ss_names'][:]).tolist()
 def faces(name):
  k=names.index(name)+1
  return [elements[int(e)-1,FACES[int(s)-1]] for e,s in zip(d['elem_ss'+str(k)][:],d['side_ss'+str(k)][:])]
 primary=faces('primary_contact');secondary=faces('secondary_contact')
mat=np.zeros((88,88));area=np.zeros(88);samples=rules()
for ids in secondary:
 face=coordinates[ids]
 face_area=np.linalg.norm(np.cross(face[1]-face[0],face[3]-face[0]))
 for local in range(8):
  transfer=samples[0 if local<4 else 1].copy()
  for v in transfer:
   x,y=v[:2]
   if local==1:v[0]=1-x
   elif local==2:v[:2]=[1-x,1-y]
   elif local==3:v[1]=1-y
   elif local==5:v[:2]=[1-y,x]
   elif local==6:v[1]=1-y
   elif local==7:v[:2]=[y,x]
  transfer[:,2]/=transfer[:,2].sum()
  weight=face_area*(1 if local<4 else 5)/24;row=ids[local];area[row]+=weight
  for x,y,w in transfer:
   n=shape(2*x-1,2*y-1)[0];p=n@face
   candidates=[(j,v) for j,f in enumerate(primary) if (v:=project(p,coordinates[f])) is not None]
   j,(_,_,natural)=min(candidates,key=lambda v:abs(v[1][1]))
   mat[row,ids]+=weight*w*n
   mat[row,primary[j]]-=weight*w*shape(*natural)[0]
predicted=mat[np.ix_(labels-1,columns-1)]/area[labels-1,None]
np.savetxt(r/'c3d20rt_finite_gap_matrix_native.csv',native,delimiter=',',header=','.join(map(str,columns)),comments='',fmt='%.12e')
np.savetxt(r/'c3d20rt_finite_gap_matrix_predicted.csv',predicted,delimiter=',',header=','.join(map(str,columns)),comments='',fmt='%.12e')
print('diagnostic_only_native_gap_matrix')
print('rows',labels.tolist());print('columns',columns.tolist())
for mask,name in [(columns<57,'primary'),(columns>=57,'secondary'),(np.ones(len(columns),dtype=bool),'all')]:
 diff=native[:,mask]-predicted[:,mask];idx=np.unravel_index(np.argmax(abs(diff)),diff.shape)
 print(name,'maximum_absolute_coefficient_difference',np.max(abs(diff)),'relative_frobenius',np.linalg.norm(diff)/np.linalg.norm(native[:,mask]),'row_node',labels[idx[0]],'column_node',columns[mask][idx[1]])
print('native_rigid_translation_max',np.max(abs(native.sum(axis=1))))
print('predicted_rigid_translation_max',np.max(abs(predicted.sum(axis=1))))

nodes=np.loadtxt(r/'c3d20rt_finite_gap_matrix_nodes.csv',delimiter=',',skiprows=1).reshape(-1,88,10)
assert len(nodes)==len(steps)
reaction=((nodes[1::2,:,6]-nodes[2::2,:,6])/(2*steps[1::2,2,None])).T
primary=columns[columns<57]
weighted=reaction[np.ix_(primary-1,np.flatnonzero(columns>=57))]@np.linalg.inv(native[:,columns>=57])
force_transfer=(weighted/weighted.sum(axis=0)).T
gap_transfer=-native[:,columns<57]
print('force_vs_gap_transfer_relative',np.linalg.norm(force_transfer-gap_transfer)/np.linalg.norm(gap_transfer))
print('force_vs_gap_transfer_maximum_absolute',np.max(abs(force_transfer-gap_transfer)))
