from __future__ import print_function

import csv
import sys
from odbAccess import openOdb


def payload(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def vector(value):
    try:
        return list(value)
    except TypeError:
        return [value]


def field(frame, name, nodal=True):
    if name not in frame.fieldOutputs:
        raise RuntimeError("Missing required field " + name)
    result = {}
    for value in frame.fieldOutputs[name].values:
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate field key in " + name)
        result[key] = payload(value)
    return result


def output_file(path):
    if sys.version_info[0] < 3:
        return open(path, "w", newline="")
    return open(path, "w", newline="")


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract.py job.odb output_prefix")

odb = openOdb(sys.argv[1], readOnly=True)
prefix = sys.argv[2]
files = [output_file(prefix + suffix) for suffix in ("_nodes.csv", "_points.csv", "_contact.csv")]
writers = [csv.writer(stream, lineterminator="\n") for stream in files]
nodes, points, contact = writers
nodes.writerow(["time", "node", "r", "z", "temperature", "ur", "uz", "rf_r", "rf_z", "reaction_heat"])
points.writerow(["time", "element", "point", "r", "z", "volume"]
                + [field_prefix + component for field_prefix in ("stress_", "elastic_")
                   for component in ("rr", "zz", "hoop", "rz")]
                + ["heat_flux_r", "heat_flux_z"])
contact.writerow(["time", "field", "node", "component", "value"])
frame_count = 0
try:
    coordinates = {}
    for instance in odb.rootAssembly.instances.values():
        for node in instance.nodes:
            if node.label in coordinates:
                raise RuntimeError("Duplicate global node label")
            coordinates[node.label] = node.coordinates
    with open(prefix + "_fields.txt", "w") as fields_file:
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = elapsed + frame.frameValue
                frame_count += 1
                fields_file.write("time=%.17g\n" % time)
                for name in sorted(frame.fieldOutputs.keys()):
                    fields_file.write(name + "\n")
                nodal = dict((name, field(frame, name)) for name in ("NT11", "U", "RF", "RFL11"))
                if set(nodal["U"]) != set(coordinates):
                    raise RuntimeError("Incomplete nodal displacement coverage")
                for node in sorted(coordinates):
                    nodes.writerow([time, node] + list(coordinates[node][:2])
                                   + [nodal["NT11"][node]] + list(nodal["U"][node][:2])
                                   + list(nodal["RF"][node][:2]) + [nodal["RFL11"][node]])
                material = dict((name, field(frame, name, False)) for name in ("S", "EE", "COORD", "IVOL", "HFL"))
                for key in sorted(material["S"]):
                    strain = list(material["EE"][key])
                    strain[3] *= 0.5
                    points.writerow([time, key[0], key[1]] + list(material["COORD"][key][:2])
                                    + [material["IVOL"][key]] + list(material["S"][key])
                                    + strain + list(material["HFL"][key][:2]))
                for name in sorted(frame.fieldOutputs.keys()):
                    short = name.split()[0]
                    if short not in ("CPRESS", "COPEN", "CSHEAR1", "CSLIP1", "CNORMF", "CSHEARF", "CSTATUS", "HFL"):
                        continue
                    if short == "HFL" and name == "HFL":
                        continue
                    for value in frame.fieldOutputs[name].values:
                        for component, number in enumerate(vector(payload(value))):
                            contact.writerow([time, name, value.nodeLabel, component, number])
            elapsed += step.timePeriod
    if frame_count != 10:
        raise RuntimeError("Expected all ten accepted increments, got %d" % frame_count)
finally:
    for stream in files:
        stream.close()
    odb.close()

print("exported_accepted_frames=%d" % frame_count)
with open(prefix + "_export_complete.txt", "w") as complete:
    complete.write("accepted_frames=%d\n" % frame_count)
