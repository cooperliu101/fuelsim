"""Extract complete prescribed-state probe output; no problem construction."""
import csv
import sys
from odbAccess import openOdb

def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data

odb = openOdb(sys.argv[1], readOnly=True)
try:
    with open('nodes.csv', 'w', newline='') as nf, open('points.csv', 'w', newline='') as pf:
        nw = csv.writer(nf); pw = csv.writer(pf)
        nw.writerow(['step', 'node', 'temperature', 'rfl', 'u1', 'u2', 'u3'])
        pw.writerow(['step', 'element', 'point', 'temperature', 'ivol'])
        for name in odb.steps.keys():
            step = odb.steps[name]
            if len(step.frames) != 2 or abs(step.frames[-1].frameValue - 1) > 1e-12:
                raise RuntimeError('Expected exactly one increment: '+name)
            frame = step.frames[-1]
            fields = {key: {v.nodeLabel: data(v) for v in frame.fieldOutputs[key].values}
                      for key in ['NT11', 'RFL11', 'U']}
            if any(len(v) != 960 for v in fields.values()):
                raise RuntimeError('Missing nodal samples')
            for node in sorted(fields['NT11']):
                u=list(fields['U'][node]);u += [0]*(3-len(u))
                nw.writerow([name,node,format(fields['NT11'][node],'.17g'),format(fields['RFL11'][node],'.17g')]+[format(x,'.17g') for x in u])
            fields = {key: {(v.elementLabel,v.integrationPoint): data(v) for v in frame.fieldOutputs[key].values}
                      for key in ['TEMP','IVOL']}
            if fields['TEMP'].keys() != fields['IVOL'].keys():
                raise RuntimeError('Inconsistent integration point coverage')
            for e,q in sorted(fields['TEMP']):
                pw.writerow([name,e,q,format(fields['TEMP'][e,q],'.17g'),format(fields['IVOL'][e,q],'.17g')])
finally:
    odb.close()
