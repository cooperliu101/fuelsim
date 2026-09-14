"""Quantify saved native thermal responses without modifying reference values."""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def heat(path):
    return np.array([float(r['reaction_heat']) for r in rows(path)])

q1 = heat(HERE/'b14_cax4t_bulk1e6_nodes.csv')
q2 = heat(HERE/'b14_cax4t_bulk1e3_nodes.csv')
# Every nodal temperature and displacement is prescribed. Bulk conduction is
# linear in its constant conductivity, so the two runs isolate that contribution.
for a,b in zip(rows(HERE/'b14_cax4t_bulk1e6_nodes.csv'),rows(HERE/'b14_cax4t_bulk1e3_nodes.csv')):
    assert all(a[key]==b[key] for key in ['time','node','temperature','ur','uz'])
bulk = (q2-q1)/(.001-.000001)
contact = q1-1e-6*bulk
index = int(np.argmax(abs(contact)))
assert abs(contact[index])>1
samples = []
variants = [('original',1e-12,1000,HERE.parent/'production/b14_nts_cax4t_nodes.csv')]
for variant,k,h in [('closed1000',1e-12,1000),('closed500',1e-12,500),('bulk1e9',1e-9,1000),
                    ('closed1e4',1e-12,.0001),('bulk1e6',1e-6,1000),('bulk1e3',1e-3,1000),
                    ('bulk1e6_half',1e-6,500),('bulk1e6_gap',1e-6,1000),
                    ('pressure_only',1e-12,1000),('pressure_first',1e-12,1000),
                    ('pressure_last',1e-12,1000),('equilibrium',1e-12,1000)]:
    variants.append((variant,k,h,HERE/('b14_cax4t_'+variant+'_nodes.csv')))
for variant,k,h,path in variants:
    q=heat(path)
    samples.append(dict(variant=variant,bulk_conductivity=k,requested_contact_conductance=h,
        sample_time=rows(path)[index]['time'],sample_node=rows(path)[index]['node'],
        native_heat_reaction_w=q[index],isolated_contact_heat_w=q[index]-k*bulk[index],
        inferred_conductance=1000*(q[index]-k*bulk[index])/contact[index]))
with (HERE/'conductance_probe.tsv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(samples[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader();writer.writerows(samples)
for sample in samples:print(sample['variant'],sample['inferred_conductance'])
# This is a diagnostic response pattern, not a production acceptance threshold.
cap=samples[0]['inferred_conductance']
for sample in samples:
    prediction=min(sample['requested_contact_conductance'],cap*sample['bulk_conductivity']/1e-12)
    assert abs(sample['inferred_conductance']/prediction-1)<1e-7,sample
