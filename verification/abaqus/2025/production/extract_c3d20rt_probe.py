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
    with open(sys.argv[2], "w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["time", "node", "temperature", "ux", "uy", "uz", "reaction_heat"])
        elapsed = 0.0
        for step in odb.steps.values():
            frame = step.frames[-1]
            fields = {}
            for name in ["NT11", "U", "RFL11"]:
                fields[name] = dict((v.nodeLabel, data(v)) for v in frame.fieldOutputs[name].values)
            for node in range(1, 9):
                writer.writerow([elapsed + frame.frameValue, node, fields["NT11"][node]] +
                                list(fields["U"][node]) + [fields["RFL11"][node]])
            elapsed += step.timePeriod
finally:
    odb.close()
