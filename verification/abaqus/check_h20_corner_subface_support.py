import sys,numpy as np
sys.path.insert(0,'verification/abaqus')
from analyze_h20_28_finite_transfer import extract
from diagnose_c3d20rt_contact_gap import rules
r=extract('h20_28_finite_transfer_32_probe');p=r[3];lookup={tuple(np.rint(64*x).astype(int)):i for i,x in enumerate(p)}
quads=[np.array([[0,0],[1/6,0],[1/8,1/8],[0,1/6]]),np.array([[1/6,0],[.5,0],[.25,.25],[.125,.125]]),np.array([[0,1/6],[.125,.125],[.25,.25],[0,.5]])]
for order in [2,3,4,5,6]:
 g,_=np.polynomial.legendre.leggauss(order);mask=np.zeros(len(p),bool);points=[]
 for q in quads:
  for a in g:
   for b in g:
    x,y=np.array([(1-a)*(1-b),(1+a)*(1-b),(1+a)*(1+b),(1-a)*(1+b)])@q/4;points.append([x,y]);i=int(32*x);j=int(32*y)
    ids=[(2*i,2*j),(2*i+2,2*j),(2*i+2,2*j+2),(2*i,2*j+2),(2*i+1,2*j),(2*i+2,2*j+1),(2*i+1,2*j+2),(2*i,2*j+1)]
    for key in ids:mask[lookup[key]]=True
 target=abs(r[2][0])>1e-4
 print(order,'points',len(points),'covered_nodes',sum(mask),'uncovered_target',sum(target&~mask),'uncovered_max',max(abs(r[2][0,~mask])))
 if order==2:print('nearest_old_sample_distance',max(min(np.linalg.norm(np.array(points)-row[:2],axis=1)) for row in rules()[0]))
