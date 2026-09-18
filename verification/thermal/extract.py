"""Run using Abaqus Python; extract native node and integration point values."""
from odbAccess import openOdb
import csv
import sys
odb = openOdb(sys.argv[1], readOnly=True)
with open(sys.argv[2], 'w') as stream:
    writer = csv.writer(stream)
    writer.writerow(['time','node','temperature','reaction'])
    for step in odb.steps.values():
        for frame in step.frames:
            if frame.frameValue <= 0:
                continue
            temperature = {v.nodeLabel: v.data for v in frame.fieldOutputs['NT11'].values}
            reaction = {v.nodeLabel: v.data for v in frame.fieldOutputs['RFL11'].values}
            for node in sorted(temperature):
                writer.writerow([frame.frameValue,node,temperature[node],reaction.get(node,0.0)])
with open(sys.argv[3], 'w') as stream:
    writer = csv.writer(stream)
    writer.writerow(['time','element','point','x','y','z','qx','qy','qz'])
    for step in odb.steps.values():
        for frame in step.frames:
            if frame.frameValue <= 0:
                continue
            coordinates = {(v.elementLabel,v.integrationPoint): list(v.data) for v in frame.fieldOutputs['COORD'].values}
            for v in frame.fieldOutputs['HFL'].values:
                xyz = coordinates[(v.elementLabel,v.integrationPoint)]
                flux = list(v.data)
                writer.writerow([frame.frameValue,v.elementLabel,v.integrationPoint]+xyz+[0.0]*(3-len(xyz))+flux+[0.0]*(3-len(flux)))
odb.close()
