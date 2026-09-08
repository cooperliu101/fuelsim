"""Independent smoothing and closure controls for the signed-transfer diagnosis."""
from pathlib import Path
import csv
import numpy as np
from analyze_h20_28_finite_transfer import extract

root=Path(__file__).resolve().parent
q=extract('h20_28_finite_transfer_probe')
l=extract('h20_28_finite_linear_transfer_probe')
s=extract('h20_28_hex20_refined_primary_sts_linear')
m=q[0]@np.linalg.inv(l[0])
print('linear_finite_vs_small_primary_relative',np.linalg.norm(l[2]-s[2])/np.linalg.norm(s[2]))
print('linear_to_quadratic_reconstruction_relative',np.linalg.norm(m@l[2]-q[2])/np.linalg.norm(q[2]))
indices=np.rint(l[3]*16).astype(int);midpoint=(indices[:,0]%2+indices[:,1]%2)==1
print('linear_minimum_mid_edge_transfer',np.min(l[2][:,midpoint]))
previous=None
for stem,closure in [('h20_28_primary_signed_gap_probe',1e-4),('h20_28_primary_signed_gap_shallow_probe',1e-5)]:
 with (root/(stem+'_gap.csv')).open() as stream:rows=list(csv.DictReader(stream))
 assert len(rows)==24
 base=np.array([float(row['gap']) for row in rows if row['step']=='BASE'])
 assert np.max(abs(base+closure))<1e-12
 plus=np.array([float(row['gap']) for row in rows if row['step']=='PRIMARY_PLUS'])
 minus=np.array([float(row['gap']) for row in rows if row['step']=='PRIMARY_MINUS'])
 transfer=-(plus-minus)/2e-7
 print(stem,'transfer',transfer.tolist())
 if previous is not None:print('closure_control_maximum_coefficient_difference',np.max(abs(transfer-previous)))
 previous=transfer
