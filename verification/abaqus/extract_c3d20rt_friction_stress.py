"""Read native contact stress components without changing the analysis."""
from __future__ import print_function
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
    with open(sys.argv[2], 'wb') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time', 'node', 'shear1', 'shear2'])
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                fields = []
                for name in ['CSHEAR1', 'CSHEAR2']:
                    keys = [k for k in frame.fieldOutputs.keys() if k.strip() == name or k.strip().startswith(name + ' ')]
                    if len(keys) != 1:
                        raise RuntimeError('Expected unique field ' + name)
                    fields.append(dict((v.nodeLabel, data(v)) for v in frame.fieldOutputs[keys[0]].values))
                for node in [57, 60, 61, 64, 68, 69, 72, 76, 78, 80, 83, 85, 88]:
                    writer.writerow([elapsed + frame.frameValue, node, fields[0][node], fields[1][node]])
            elapsed += step.timePeriod
finally:
    odb.close()
