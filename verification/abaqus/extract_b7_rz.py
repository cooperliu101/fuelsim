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
    if name not in list(frame.fieldOutputs.keys()):
        return {}
    result = {}
    for value in frame.fieldOutputs[name].values:
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate output entry in " + name)
        result[key] = data(value)
    return result


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b7_rz.py job.odb output_prefix")
odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2] + "_nodes.csv", "wb") as node_file, open(sys.argv[2] + "_points.csv", "wb") as point_file:
        nodes = csv.writer(node_file)
        points = csv.writer(point_file)
        contact = sys.argv[2].endswith("_contact")
        nodes.writerow(["time", "node", "temperature", "ur", "uz", "rf_r", "rf_z", "reaction_heat"] +
                       (["contact_pressure", "contact_gap", "contact_force_z", "contact_heat_flux"] if contact else []))
        components = ["rr", "zz", "hoop", "rz"]
        points.writerow(["time", "element", "point"] +
                        [prefix + c for prefix in ["stress_", "elastic_", "plastic_", "creep_"] for c in components] +
                        ["equiv_plastic", "equiv_creep"])
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = elapsed + frame.frameValue
                fields = dict((name, values(frame, name, True)) for name in ["NT11", "U", "RF", "RFL11"])
                if contact:
                    contact_fields = {}
                    for short in ["CPRESS", "COPEN", "CNORMF", "HFL"]:
                        keys = [key for key in frame.fieldOutputs.keys()
                                if key.split()[0] == short and "S_LOWER/S_UPPER" in key]
                        if len(keys) != 1:
                            raise RuntimeError("Missing or ambiguous contact field " + short)
                        contact_fields[short] = values(frame, keys[0], True)
                for node in sorted(fields["U"]):
                    row = [time, node, fields["NT11"][node]] + list(fields["U"][node][:2]) +\
                          list(fields["RF"][node][:2]) + [fields["RFL11"][node]]
                    if contact:
                        # COPEN and interface HFL exist on secondary nodes only.
                        if node in (3, 4):
                            row += [contact_fields["CPRESS"][node], contact_fields["COPEN"][node],
                                    contact_fields["CNORMF"][node][1], contact_fields["HFL"][node]]
                        else:
                            row += [0.0, 0.0, contact_fields["CNORMF"].get(node, (0.0, 0.0))[1], 0.0]
                    nodes.writerow(row)
                fields = dict((name, values(frame, name)) for name in ["S", "EE", "PE", "CE", "PEEQ", "CEEQ"])
                for key in sorted(fields["S"]):
                    row = [time, key[0], key[1]]
                    for name in ["S", "EE", "PE", "CE"]:
                        tensor = list(fields[name].get(key, (0.0, 0.0, 0.0, 0.0)))
                        if name != "S":
                            tensor[3] *= 0.5
                        row += tensor
                    row += [fields["PEEQ"].get(key, 0.0), fields["CEEQ"].get(key, 0.0)]
                    points.writerow(row)
            elapsed += step.timePeriod
finally:
    odb.close()
