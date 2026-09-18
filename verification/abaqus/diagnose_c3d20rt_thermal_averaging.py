"""Validate corner-patch averaged thermal constraints against two independent Abaqus matrices."""
import csv
from pathlib import Path
import re
import numpy as np

source = Path('src/core/cartesian3d_assembly.cpp').read_text() + Path('elements/src/quad8_face.cpp').read_text()
rules = []
for name in ('corner', 'edge'):
    begin = source.index('static const std::array<AbaqusQuad8TransferSample',
                         source.index('abaqus_quad8_' + name + '_transfer_rule()'))
    end = source.index('return value;', begin)
    rule = np.array([float(value) for value in re.findall(
        r'-?\d+\.\d+(?:e[+-]?\d+)?', source[begin:end])]).reshape(-1, 3)
    rule[:, 2] /= rule[:, 2].sum()
    rules.append(rule)
# A sample on a quadrant boundary contributes equally to its adjacent patches.
# The tolerance only identifies coincident samples in the independently recovered rule.
A=np.zeros((4,4));D=np.zeros(4)
for i in range(8):
 for u,v,w in rules[i//4]:
  if i in [1,2]:u=1-u
  if i in [2,3,6]:v=1-v
  if i==5:u,v=1-v,u
  if i==7:u,v=v,u
  n=np.array([(1-u)*(1-v),u*(1-v),(1-u)*v,u*v])
  ids=[j for j in range(4) if ((u<=.5+1e-12) if j%2==0 else (u>=.5-1e-12)) and ((v<=.5+1e-12) if j<2 else (v>=.5-1e-12))]
  for j in ids:
   weight=w*(1/24 if i<4 else 5/24)/len(ids);A[j]+=weight*n;D[j]+=weight
A/=D[:,None]
C=np.loadtxt('verification/abaqus/c3d20rt_contact_matching_matrix.csv',delimiter=',',skiprows=1)[:,1:]
matching_error = np.max(abs(1.2*A.T@np.diag(D)@A-C[:4,:4]))
print('matching_patch_areas', D.tolist())
print('matching_temperature_average', A.tolist())
print('matching_maximum_absolute_matrix_error', matching_error)
assert matching_error < 1e-10
# Merge adjacent-face patch integrals at shared temperature corner nodes
# before forming B.T * inverse(diag(area)) * B.
B=np.zeros((6,16));D=np.zeros(6)
for side,(lo,hi) in enumerate([(0,1.2),(1.2,2)]):
 ss=([10,11,12,13] if side==0 else [11,14,13,15])
 for i in range(8):
  for u,v,w in rules[i//4]:
   if i in [1,2]:u=1-u
   if i in [2,3,6]:v=1-v
   if i==5:u,v=1-v,u
   if i==7:u,v=v,u
   y=lo+(hi-lo)*u;j=min(int(y/.5),3);f=y/.5-j;d=np.zeros(16)
   remap=[0,2,1,3,4,5,6,7,8,9]
   ps=[remap[k] for k in [2*j,2*(j+1),2*j+1,2*(j+1)+1]]
   d[ps]=-np.array([(1-f)*(1-v),f*(1-v),(1-f)*v,f*v]);d[ss]=[(1-u)*(1-v),u*(1-v),(1-u)*v,u*v]
   ids=[q for q in range(4) if ((u<=.5+1e-12) if q%2==0 else (u>=.5-1e-12)) and ((v<=.5+1e-12) if q<2 else (v>=.5-1e-12))]
   for q in ids:
    weight=(hi-lo)*w*(1/24 if i<4 else 5/24)/len(ids);idx=ss[q]-10;B[idx]+=weight*d;D[idx]+=weight
C=np.loadtxt('verification/abaqus/c3d20rt_contact_thermal_matrix.csv',delimiter=',',skiprows=1)[:,1:]
nonmatching_error = np.max(abs(B.T@np.diag(1/D)@B-C))
print('nonmatching_patch_areas', D.tolist())
print('nonmatching_maximum_absolute_matrix_error', nonmatching_error)
assert nonmatching_error < 1e-10
