"""Identify averaged primary transfer from prescribed secondary perturbations."""
from pathlib import Path
import csv
import argparse
import numpy as np

ROOT=Path(__file__).resolve().parent

def extract(name):
 with (ROOT/(name+'_operator.csv')).open() as f: rows=list(csv.DictReader(f))
 grouped={}
 for row in rows:grouped.setdefault((row['step'],row['side']),[]).append(row)
 def values(step,side,field):return np.array([float(row[field]) for row in grouped[(step,side)]])
 delta=np.array([float(grouped[('S%d_PLUS'%j,'secondary')][0]['closure_delta_m'])-float(grouped[('S%d_MINUS'%j,'secondary')][0]['closure_delta_m']) for j in range(1,9)])
 if not np.all(np.isfinite(delta)) or np.any(delta<=0):raise ValueError('Invalid native perturbation increments')
 a=np.column_stack([-(values('S%d_PLUS'%j,'secondary','copen_m')-values('S%d_MINUS'%j,'secondary','copen_m'))/delta[j-1] for j in range(1,9)])
 reaction=np.column_stack([(values('S%d_PLUS'%j,'primary','cnormf_x_n')-values('S%d_MINUS'%j,'primary','cnormf_x_n'))/delta[j-1] for j in range(1,9)])
 weighted=reaction@np.linalg.inv(a)
 area=weighted.sum(axis=0)/1e8
 transfer=(weighted/weighted.sum(axis=0)).T
 coordinates=np.array([[float(row['coord_y_m']),float(row['coord_z_m'])] for row in grouped[('BASE','primary')]])
 secondary=np.array([[float(row['coord_y_m']),float(row['coord_z_m'])] for row in grouped[('BASE','secondary')]])
 return a,area,transfer,coordinates,secondary

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--refined-grid',type=int,choices=[16,32],default=16)
 grid=parser.parse_args().refined_grid
 small=extract('h20_28_hex20_refined_primary_sts_default')
 finite=extract('h20_28_finite_transfer_probe')
 print('diagnostic_only_refined_primary_transfer')
 for i,name in enumerate(['secondary_averaging','constraint_area','primary_transfer']):
  a,b=finite[i],small[i]
  print(name,'maximum_absolute_difference',np.max(abs(a-b)),'relative_frobenius',np.linalg.norm(a-b)/np.linalg.norm(b))
 for name,result in [('small',small),('finite',finite)]:
  print(name,'areas',result[1].tolist())
  print(name,'nonzero_primary_coefficients',np.sum(abs(result[2])>1e-8,axis=1).tolist())
  print(name,'effective_centers',(result[2]@result[3]).tolist())
 np.savetxt(ROOT/'h20_28_finite_transfer_matrix.csv',finite[2],delimiter=',',fmt='%.12e')
 np.savetxt(ROOT/'h20_28_finite_transfer_primary_coordinates.csv',finite[3],delimiter=',',header='y,z',comments='',fmt='%.12e')
 np.savetxt(ROOT/'h20_28_finite_transfer_secondary_coordinates.csv',finite[4],delimiter=',',header='y,z',comments='',fmt='%.12e')
 refined=extract('h20_28_finite_transfer_'+str(grid)+'_probe')
 print('refined_grid',grid)
 print('refined_secondary_averaging_max_difference',np.max(abs(refined[0]-finite[0])))
 print('refined_area_max_difference',np.max(abs(refined[1]-finite[1])))
 for tolerance in [1e-8,1e-6,1e-4]:
  print('refined_nonzero',tolerance,np.sum(abs(refined[2])>tolerance,axis=1).tolist())
 print('refined_effective_centers',(refined[2]@refined[3]).tolist())
 # Rows are normalized during extraction; this is arithmetic QA, not an
 # independent physical rigid-translation check.
 print('refined_normalization_roundoff',np.max(abs(refined[2].sum(axis=1)-1)))
 print('refined_transfer_row_l1',np.sum(abs(refined[2]),axis=1).tolist())
 indices=np.rint(refined[3]*2*grid).astype(int)
 midpoint=(indices[:,0]%2+indices[:,1]%2)==1
 print('refined_minimum_mid_edge_coefficients',np.min(refined[2][:,midpoint],axis=1).tolist())
 np.savetxt(ROOT/('h20_28_finite_transfer_'+str(grid)+'_matrix.csv'),refined[2],delimiter=',',fmt='%.12e')
 np.savetxt(ROOT/('h20_28_finite_transfer_'+str(grid)+'_primary_coordinates.csv'),refined[3],delimiter=',',header='y,z',comments='',fmt='%.12e')


if __name__=='__main__':main()
