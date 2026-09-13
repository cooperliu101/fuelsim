from __future__ import print_function
import csv
import sys
from odbAccess import openOdb

odb = openOdb(sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2], 'w', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(['time', 'element', 'point', 'heat_r', 'heat_z'])
        elapsed = 0.0
        for step in odb.steps.values():
            for index in range(1, len(step.frames)):
                frame = step.frames[index]
                for value in frame.fieldOutputs['HFL'].values:
                    try:
                        vector = value.dataDouble
                    except Exception:
                        vector = value.data
                    writer.writerow([elapsed + frame.frameValue, value.elementLabel,
                                     value.integrationPoint, vector[0], vector[1]])
            elapsed += step.timePeriod
finally:
    odb.close()
