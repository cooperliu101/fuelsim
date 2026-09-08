"""Identify the friction virtual-work geometry at fixed Abaqus states.

Uses native accumulated slip and the unchanged Coulomb law. The local gap and
surface gradient are computed from geometry. No coefficient is fitted to force.
This prescribed-state diagnostic does not replace a production solve.
"""
import sys
sys.path.insert(0,'verification/abaqus')
import diagnose_c3d20rt_finite_friction as d
from diagnose_c3d20rt_contact_gap import project
import numpy as np
import netCDF4 as nc
with nc.Dataset(d.ROOT/'verification/meshes/c3d20rt_finite_contact.e') as ds:
 i=nc.chartostring(ds['ss_names'][:]).tolist().index('primary_contact')+1
 primaryids=[d.elements[int(e)-1,d.FACES[int(s)-1]] for e,s in zip(ds['elem_ss'+str(i)][:],ds['side_ss'+str(i)][:])]
for step in [0,39,159,170,207,279]:
 matrices={k:np.zeros((13,13,2,3)) for k in ['native_convention']}
 for ids in d.faces:
  face=d.X[step,ids];target=[d.lookup[int(i)] for i in ids]
  for i,row in enumerate(target):
   for a,b,w in d.samples(i):
    n,dx,dy=d.shape(2*a-1,2*b-1);tx=dx@face;ty=dy@face
    sn=-np.cross(tx,ty);measure=np.linalg.norm(sn);sn/=measure
    candidates=[]
    for pids in primaryids:
     pface=d.X[step,pids];pr=project(n@face,pface)
     if pr is not None:
      _,pdx,pdy=d.shape(*pr[2]);pn=np.cross(pdx@pface,pdy@pface);pn/=np.linalg.norm(pn)
      if pn@sn<0:pn=-pn
      candidates.append((abs(pr[1]),pn,pr[0]-n@face))
    if not candidates: raise RuntimeError('No primary normal')
    _,pn,rel=min(candidates,key=lambda v:v[0]);gp=-rel@sn
    normal=sn;A=np.array([tx,ty]);grad=np.array([dx,dy]).T@np.linalg.solve(A@A.T,A)
    wt=4*(1 if i<4 else 5)/24*w*measure
    for mode in matrices:
     t1=ty/np.linalg.norm(ty);t2=np.cross(normal,t1);basis=np.array([t1,t2])
     matrices[mode][row,target]+=wt*(n[:,None,None]*basis[None,:,:]+gp*(grad@basis.T)[:,:,None]*normal[None,None,:])
 for key,matrix in matrices.items():
  M=matrix.transpose(1,3,0,2).reshape(39,26);force=-d.contact[step,:,7:10].reshape(39)
  value=d.traction_components[step].ravel();error=M@value-force
  print(step+1,key,'relative',np.linalg.norm(error)/np.linalg.norm(force),'max',abs(error).max(),flush=True)
