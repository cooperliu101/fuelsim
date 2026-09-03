from __future__ import print_function

import json
import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError("usage: extract_b60_finite.py <job.odb> <mesh.json> <nodal.csv> <final step time>")


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field(frame, name):
    if name not in frame.fieldOutputs:
        raise RuntimeError("missing Abaqus field %s" % name)
    return frame.fieldOutputs[name]


def nodal_values(output, labels):
    result = {}
    for value in output.values:
        if value.nodeLabel in labels:
            if value.nodeLabel in result:
                raise RuntimeError("duplicate nodal value for node %d" % value.nodeLabel)
            result[value.nodeLabel] = data(value)
    if len(result) != len(labels):
        raise RuntimeError("incomplete nodal field: expected %d values, got %d" % (len(labels), len(result)))
    return result


with open(sys.argv[2], "rb") as source:
    mesh = json.load(source)
coordinates = dict((node["label"], node["coordinates"]) for node in mesh["nodes"])
labels = set(coordinates.keys())

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    if len(odb.steps) != 1:
        raise RuntimeError("B6.0 finite-strain expects exactly one Abaqus step")
    step = list(odb.steps.values())[0]
    frame = step.frames[-1]
    final_time = float(sys.argv[4])
    if abs(frame.frameValue - final_time) > 1.0e-10:
        raise RuntimeError("B6.0 finite-strain final Abaqus frame has the wrong step time")
    temperature = nodal_values(field(frame, "NT11"), labels)
    displacement = nodal_values(field(frame, "U"), labels)
    output = open(sys.argv[3], "wb")
    output.write("id,x,y,z,temperature,displacement_x,displacement_y,displacement_z\n")
    for label in sorted(labels):
        output.write(
            "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n"
            % ((label,) + tuple(coordinates[label]) + (temperature[label],) + tuple(displacement[label]))
        )
    output.close()
finally:
    odb.close()
