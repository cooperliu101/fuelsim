"""Read-only native-state friction replay; no production solve or fitted parameters.

Separate secondary force transfer and Coulomb history from geometric slip
integration. Uses the previously identified Q8 integration samples.
"""
from pathlib import Path
import numpy as np
import netCDF4 as nc
from diagnose_c3d20rt_contact_gap import FACES, LOCATIONS, shape
from check_h20_disk_transfer import samples

ROOT=Path(__file__).resolve().parents[2]
STEM=ROOT/'verification/abaqus/c3d20rt_finite_contact_friction'
nodes=np.loadtxt(str(STEM)+'_nodes.csv',delimiter=',',skiprows=1).reshape(-1,88,10)
contact=np.loadtxt(str(STEM)+'_contact.csv',delimiter=',',skiprows=1).reshape(-1,13,14)
frames=np.loadtxt(str(STEM)+'_frames.csv',delimiter=',',skiprows=1).reshape(-1,13,10)
assert np.array_equal(contact[:,:,:2],frames[:,:,:2])
with nc.Dataset(ROOT/'verification/meshes/c3d20rt_finite_contact.e') as d:
 x0=np.array([d['coord'+a][:] for a in 'xyz']).T
 elements=np.vstack([d['connect1'][:],d['connect2'][:]])-1
 i=nc.chartostring(d['ss_names'][:]).tolist().index('secondary_contact')+1
 faces=[elements[int(e)-1,FACES[int(s)-1]] for e,s in zip(d['elem_ss'+str(i)][:],d['side_ss'+str(i)][:])]
labels=contact[0,:,1].astype(int)-1; lookup={int(n):i for i,n in enumerate(labels)}
X=x0[None,:,:]+nodes[:,:,3:6]; nframes=len(X)
scalar=np.zeros((nframes,13,13)); localnorm=np.zeros((nframes,13,3)); areanorm=np.zeros_like(localnorm)
tangent_transfer=np.zeros((nframes,13,13,2,3))
normal_transfer=np.zeros((nframes,13,13,3))
for ids in faces:
 face=X[:,ids]; target=[lookup[int(i)] for i in ids]
 for i,row in enumerate(target):
  q=samples(i); shapes=np.array([shape(2*a-1,2*b-1) for a,b,w in q])
  tx=np.einsum('qn,tnc->tqc',shapes[:,1],face);ty=np.einsum('qn,tnc->tqc',shapes[:,2],face)
  av=np.cross(tx,ty);measure=np.linalg.norm(av,axis=2);weights=4*(1 if i<4 else 5)/24*q[:,2]
  c=np.einsum('q,tq,qn->tn',weights,measure,shapes[:,0]);scalar[:,row,target]+=c
  pn=-av/measure[:,:,None];pt2=tx/np.linalg.norm(tx,axis=2)[:,:,None];pt1=-np.cross(pn,pt2)
  tb=np.stack([pt1,pt2],axis=2)
  normal_transfer[:,row,target]+=np.einsum('q,tq,qn,tqc->tnc',weights,measure,shapes[:,0],pn)
  tangent_transfer[:,row,target]+=np.einsum('q,tq,qn,tqkc->tnkc',weights,measure,shapes[:,0],tb)
  areanorm[:,row]+=np.einsum('q,tqc->tc',weights,av)
  _,dx,dy=shape(*LOCATIONS[i]);n=np.cross(np.einsum('n,tnc->tc',dx,face),np.einsum('n,tnc->tc',dy,face));n/=np.linalg.norm(n,axis=1)[:,None]
  localnorm[:,row]+=np.sum(weights[None,:]*measure,axis=1)[:,None]*n
localnorm/=np.linalg.norm(localnorm,axis=2)[:,:,None];areanorm/=np.linalg.norm(areanorm,axis=2)[:,:,None]
basis=frames[:,:,2:8].reshape(nframes,13,2,3);normal=np.cross(basis[:,:,0],basis[:,:,1])
for name,n in [('constraint_location_normal',localnorm),('integrated_area_normal',areanorm)]:
 print(name,'maximum_native_difference',np.min(np.stack([np.linalg.norm(n-normal,axis=2),np.linalg.norm(n+normal,axis=2)]),axis=0).max())
# Invert the identified scalar transfer to diagnose whether recovered traction is tangential.
traction=np.linalg.solve(scalar.transpose(0,2,1),-contact[:,:,7:10])
print('scalar_transfer_recovered_normal_traction_maximum_pa',np.max(abs(np.sum(traction*normal,axis=2))))
print('scalar_transfer_recovered_traction_peak_pa',np.max(np.linalg.norm(traction,axis=2)))
length=1e-5*(np.sqrt(1.2)+np.sqrt(.8))/2
old=np.zeros((13,2));last=np.zeros((13,2));prediction=[];sampled_prediction=[];traction_components=[]
for step in range(nframes):
 limit=-.3e11*contact[step,:,2];stiffness=limit/length
 slip=frames[step,:,8:10];trial=(old+slip-last)*stiffness[:,None]
 magnitude=np.linalg.norm(trial,axis=1);scale=np.minimum(1,limit/np.maximum(magnitude,np.finfo(float).tiny))
 value=trial*scale[:,None];old=value/stiffness[:,None];last=slip
 traction_components.append(value.copy())
 global_value=np.einsum('ik,ikc->ic',value,basis[step])
 prediction.append(-scalar[step].T@global_value)
 sampled_prediction.append(-np.einsum('ijkc,ik->jc',tangent_transfer[step],value))
error=np.array(prediction)-contact[:,:,7:10]
print('native_slip_scalar_transfer_force_relative_l2',np.linalg.norm(error)/np.linalg.norm(contact[:,:,7:10]))
print('native_slip_scalar_transfer_force_maximum_absolute_n',np.max(np.linalg.norm(error,axis=2)))
print('native_slip_scalar_transfer_force_maximum_relative',np.max(np.linalg.norm(error,axis=2)/np.linalg.norm(contact[:,:,7:10],axis=2)))

error=np.array(sampled_prediction)-contact[:,:,7:10]
print('native_slip_sampled_tangent_force_relative_l2',np.linalg.norm(error)/np.linalg.norm(contact[:,:,7:10]))
print('native_slip_sampled_tangent_force_maximum_absolute_n',np.max(np.linalg.norm(error,axis=2)))
print('native_slip_sampled_tangent_force_maximum_relative',np.max(np.linalg.norm(error,axis=2)/np.linalg.norm(contact[:,:,7:10],axis=2)))
for step in [0,39,79,159,207,279]:
 matrix=tangent_transfer[step].transpose(1,3,0,2).reshape(39,26)
 value,residual,rank,sv=np.linalg.lstsq(matrix,-contact[step,:,7:10].reshape(39),rcond=None)
 error=(matrix@value).reshape(13,3)+contact[step,:,7:10]
 print('sampled_tangent_recovery',step+1,'relative_l2',np.linalg.norm(error)/np.linalg.norm(contact[step,:,7:10]),'max_abs_n',np.max(abs(error)))

normal_prediction=np.einsum('ti,tijc->tjc',-1e11*contact[:,:,2],normal_transfer)
normal_error=normal_prediction-contact[:,:,4:7]
print('native_normal_force_relative_l2',np.linalg.norm(normal_error)/np.linalg.norm(contact[:,:,4:7]))
print('native_normal_force_maximum_absolute_n',abs(normal_error).max())
