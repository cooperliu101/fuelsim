from __future__ import print_function
import csv
import sys
from odbAccess import openOdb


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def values(frame, name, nodal=False):
    if name not in frame.fieldOutputs:
        return {}
    result = {}
    for value in frame.fieldOutputs[name].values:
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate output entry in " + name)
        result[key] = data(value)
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2] + "_nodes.csv", "w", newline="") as nf, open(sys.argv[2] + "_points.csv", "w", newline="") as pf:
        nodes, points = csv.writer(nf), csv.writer(pf)
        nodes.writerow(["time", "node", "temperature", "ur", "uz", "rf_r", "rf_z", "reaction_heat", "temperature_active"])
        components = ["rr", "zz", "hoop", "rz"]
        points.writerow(["time", "element", "point"] +
                        [prefix + c for prefix in ["stress_", "elastic_", "plastic_", "creep_"] for c in components] +
                        ["equiv_plastic", "equiv_creep", "radius", "axial", "volume"])
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = elapsed + frame.frameValue
                fields = dict((n, values(frame, n, True)) for n in ["NT11", "U", "RF", "RFL11"])
                for node in sorted(fields["U"]):
                    active = node in fields["RFL11"]
                    nodes.writerow([time, node, fields["NT11"].get(node, float("nan"))] +
                                   list(fields["U"][node][:2]) + list(fields["RF"][node][:2]) +
                                   [fields["RFL11"].get(node, float("nan")), int(active)])
                fields = dict((n, values(frame, n)) for n in ["S", "EE", "PE", "CE", "PEEQ", "CEEQ", "COORD", "IVOL"])
                for key in sorted(fields["S"]):
                    row = [time, key[0], key[1]]
                    for name in ["S", "EE", "PE", "CE"]:
                        tensor = list(fields[name].get(key, (0., 0., 0., 0.)))
                        if name != "S":
                            tensor[3] *= .5
                        row += tensor
                    row += [fields["PEEQ"].get(key, 0.), fields["CEEQ"].get(key, 0.)]
                    row += list(fields["COORD"].get(key, (float("nan"), float("nan")))[:2])
                    row += [fields["IVOL"].get(key, float("nan"))]
                    points.writerow(row)
            elapsed += step.timePeriod
finally:
    odb.close()
