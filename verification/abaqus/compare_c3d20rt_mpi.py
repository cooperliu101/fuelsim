import netCDF4 as nc
import numpy as np
import sys
if len(sys.argv) != 3:
 raise SystemExit('Expected one-rank and two-rank Exodus result paths')
paths=sys.argv[1:]

with nc.Dataset(paths[0]) as a,nc.Dataset(paths[1]) as b:
 names=nc.chartostring(a['name_nod_var'][:]).tolist()
 assert names==nc.chartostring(b['name_nod_var'][:]).tolist()
 assert np.array_equal(a['time_whole'][:],b['time_whole'][:])
 assert all(name in names for name in ['temperature','displacement_x','displacement_y','displacement_z'])
 for i,name in enumerate(names):
  if name not in ['temperature','displacement_x','displacement_y','displacement_z']:continue
  x=np.asarray(a[f'vals_nod_var{i+1}'][:]);y=np.asarray(b[f'vals_nod_var{i+1}'][:]);assert np.array_equal(np.isnan(x),np.isnan(y))
  error=float(np.nanmax(abs(x-y)));print(name+'_all_frame_maximum_absolute_difference',error)
  assert error<(1e-10 if name=='temperature' else 1e-12)
 print('compared_frames',len(a['time_whole'][:]))
 print('compared_nodes',len(a.dimensions['num_nodes']))
