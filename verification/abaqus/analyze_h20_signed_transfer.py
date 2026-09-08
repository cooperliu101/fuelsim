"""Check signed finite-sliding transfer using an independent primary perturbation."""
from pathlib import Path
import csv
import numpy as np
from analyze_h20_28_finite_transfer import extract

root=Path(__file__).resolve().parent
for stem,grid in [('h20_28_hex20_refined_primary_sts_default',8),
                  ('h20_28_finite_transfer_probe',8),('h20_28_finite_transfer_16_probe',16)]:
 result=extract(stem);indices=np.rint(result[3]*2*grid).astype(int)
 midpoint=(indices[:,0]%2+indices[:,1]%2)==1
 for row in [0,4]:
  j=np.flatnonzero(midpoint)[np.argmin(result[2][row,midpoint])]
  print(stem,'constraint_row',row,'minimum_mid_edge_coefficient',result[2][row,j],
        'primary_yz',result[3][j].tolist())

result=extract('h20_28_finite_transfer_probe')
column=np.flatnonzero(np.all(abs(result[3]-[.3125,0])<1e-12,axis=1))
assert len(column)==1
with (root/'h20_28_primary_signed_gap_probe_gap.csv').open() as stream:rows=list(csv.DictReader(stream))
assert len(rows)==24
labels=[int(row['node']) for row in rows if row['step']=='BASE']
for name in ['PRIMARY_PLUS','PRIMARY_MINUS']:
 assert [int(row['node']) for row in rows if row['step']==name]==labels
plus=np.array([float(row['gap']) for row in rows if row['step']=='PRIMARY_PLUS'])
minus=np.array([float(row['gap']) for row in rows if row['step']=='PRIMARY_MINUS'])
direct=-(plus-minus)/2e-7
force=result[2][:,column[0]]
print('primary_node_170_direct_gap_transfer',direct.tolist())
print('primary_node_170_force_transfer',force.tolist())
print('maximum_absolute_difference',np.max(abs(direct-force)))
assert direct[0]<-0.003 and np.max(abs(direct-force))<1e-7
