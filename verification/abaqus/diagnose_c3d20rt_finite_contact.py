"""Compare only mutually completed finite-contact increments; this is not acceptance."""
from pathlib import Path
import netCDF4 as nc
import numpy as np
root=Path(__file__).resolve().parent
nodes=np.loadtxt(root/'c3d20rt_finite_contact_partial_nodes.csv',delimiter=',',skiprows=1).reshape(-1,88,10)
contact=np.loadtxt(root/'c3d20rt_finite_contact_partial_contact.csv',delimiter=',',skiprows=1).reshape(-1,13,14)
with nc.Dataset(root.parent/'fuelsim/transient_c3d20rt_finite_contact_equal_height_results.e') as data:
 names=nc.chartostring(data['name_nod_var'][:]).tolist()
 def field(name,frame):return np.asarray(data['vals_nod_var'+str(names.index(name)+1)][frame])
 times=np.asarray(data['time_whole'][:]);count=min(len(times)-1,len(nodes));print('diagnostic_only_common_increments',count)
 print('time,temperature_relative_l2,displacement_relative_l2,pressure_relative_l2,gap_relative_l2,normal_force_relative_l2,tangential_force_relative_l2')
 for frame in range(1,count+1):
  reference=nodes[frame-1];rcontact=contact[frame-1];indices=rcontact[:,1].astype(int)-1
  assert abs(times[frame]-reference[0,0])<1e-6
  pairs=[(field('temperature',frame),reference[:,2]),(np.array([field('displacement_'+axis,frame) for axis in 'xyz']).T,reference[:,3:6]),(field('contact_pressure_interface',frame)[indices],rcontact[:,3]),(field('contact_gap_interface',frame)[indices],rcontact[:,2]),(np.array([field('contact_normal_force_'+axis+'_interface',frame)[indices] for axis in 'xyz']).T,-rcontact[:,4:7]),(np.array([field('contact_tangential_force_'+axis+'_interface',frame)[indices] for axis in 'xyz']).T,-rcontact[:,7:10])]
  errors=[np.linalg.norm(a-b)/np.linalg.norm(b) for a,b in pairs]
  print(','.join([format(times[frame],'.8g')]+[format(v,'.12e') for v in errors]))
