"""Check uniform-temperature native histories against the exact zero heat flux."""
import csv
from pathlib import Path
import netCDF4
import numpy as np

HERE=Path(__file__).resolve().parent
rows=[]
for temperature,path in [(300,HERE/'b114_uniform_300_nodes.csv'),
        (600,HERE.parent/'production/b114_cax8t_recovery_nodes.csv'),
        (1200,HERE/'b114_uniform_1200_nodes.csv')]:
    with path.open() as stream: data=list(csv.DictReader(stream))
    assert all(float(row['temperature'])==temperature for row in data)
    times=sorted({float(row['time']) for row in data})
    assert len(times)==14
    values=np.array([float(row['reaction_heat']) for row in data])
    rows.append(dict(solver='Abaqus2025',uniform_temperature_k=temperature,
        frames=len(times),samples=len(data),maximum_absolute_reaction_w=max(abs(values)),
        maximum_absolute_frame_sum_w=max(abs(sum(float(row['reaction_heat']) for row in data
            if float(row['time'])==time)) for time in times)))
root=HERE.parents[3]
with netCDF4.Dataset(root/'build/blackbox/b114_cax8t_recovery/verification/fuelsim/transient_b114_cax8t_recovery_results.e') as ds:
    names=netCDF4.chartostring(ds['name_nod_var'][:]).tolist()
    def field(name): return np.array(ds['vals_nod_var%d'%(names.index(name)+1)][1:])
    assert np.max(abs(field('temperature')-600))==0
    q=field('reaction_heat_flux')
    rows.append(dict(solver='Fuelsim',uniform_temperature_k=600,frames=q.shape[0],samples=q.size,
        maximum_absolute_reaction_w=np.max(abs(q)),maximum_absolute_frame_sum_w=np.max(abs(q.sum(axis=1)))))
with (HERE/'b114_uniform_temperature.tsv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
for row in rows: print(row)
