"""Compare uniform tiny-face transfer with independent Quad8 area integration."""
from pathlib import Path
import csv
import numpy as np
from diagnose_c3d20rt_contact_gap import shape

root=Path(__file__).resolve().parent
for stem,y,z in [('h20_28_finite_point_probe',.18,.18),('h20_28_finite_point_edge_probe',.1249,.18)]:
 with (root/(stem+'_operator.csv')).open() as stream:
  rows=[row for row in csv.DictReader(stream) if row['step']=='BASE' and row['side']=='primary']
 assert len(rows)==225
 xy=np.array([[float(row['coord_y_m']),float(row['coord_z_m'])] for row in rows])
 native=np.array([float(row['cnormf_x_n']) for row in rows]);native/=native.sum()
 lookup={tuple(np.rint(16*p).astype(int)):i for i,p in enumerate(xy)}
 predicted=np.zeros(225);g,w=np.polynomial.legendre.leggauss(3)
 for u,wu in zip(g,w):
  for v,wv in zip(g,w):
   a=y+5e-5*u;b=z+5e-5*v;i=int(8*a);j=int(8*b)
   n=shape(2*(8*a-i)-1,2*(8*b-j)-1)[0]
   ids=[lookup[c] for c in [(2*i,2*j),(2*i+2,2*j),(2*i+2,2*j+2),(2*i,2*j+2),(2*i+1,2*j),(2*i+2,2*j+1),(2*i+1,2*j+2),(2*i,2*j+1)]]
   predicted[ids]+=.25*wu*wv*n
 print(stem,'relative_l2',np.linalg.norm(native-predicted)/np.linalg.norm(predicted),
       'maximum_absolute_coefficient_difference',np.max(abs(native-predicted)),
       'native_nonzero_above_1e-8',np.count_nonzero(abs(native)>1e-8))
