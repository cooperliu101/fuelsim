from __future__ import print_function
import csv
import sys
from odbAccess import openOdb

def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data

def field(frame, name, nodal=False):
    if name not in list(frame.fieldOutputs.keys()):
        return {}
    return dict(((v.nodeLabel if nodal else (v.elementLabel, v.integrationPoint)), data(v))
                for v in frame.fieldOutputs[name].values)

odb = openOdb(sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2] + "_nodes.csv", "w", newline="") as nf, open(sys.argv[2] + "_points.csv", "w", newline="") as pf:
        nodes, points = csv.writer(nf), csv.writer(pf)
        nodes.writerow(["time", "node", "temperature", "ux", "uy", "uz", "rx", "ry", "rz"])
        components = ["xx", "yy", "zz", "xy", "yz", "xz"]
        points.writerow(["time", "element", "point"] +
                        [prefix + c for prefix in ["s_", "e_", "p_", "c_"] for c in components] +
                        ["equiv_plastic", "equiv_creep"])
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = elapsed + frame.frameValue
                temp, disp, force = [field(frame, name, True) for name in ["NT11", "U", "RF"]]
                for node in sorted(disp):
                    nodes.writerow([time, node, temp[node]] + list(disp[node]) + list(force[node]))
                fields = dict((name, field(frame, name)) for name in ["S", "EE", "PE", "CE", "PEEQ", "CEEQ"])
                for key in sorted(fields["S"]):
                    row = [time, key[0], key[1]]
                    for name in ["S", "EE", "PE", "CE"]:
                        vector = fields[name].get(key, [0.0] * 6)
                        # Abaqus tensor output order is 11,22,33,12,13,23.
                        row += [vector[i] * (0.5 if name != "S" and i >= 3 else 1.0)
                                for i in [0, 1, 2, 3, 5, 4]]
                    row += [fields[name].get(key, 0.0) for name in ["PEEQ", "CEEQ"]]
                    points.writerow(row)
            elapsed += step.timePeriod
finally:
    odb.close()
