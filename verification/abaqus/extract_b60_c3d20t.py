from __future__ import print_function

import json
import math
import sys

from odbAccess import openOdb


if len(sys.argv) != 6:
    raise RuntimeError(
        "usage: extract_b60_c3d20t.py <job.odb> <mesh.json> <temperature.csv> <displacement.csv> <material.csv>"
    )


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


def integration_values(output, labels):
    result = {}
    for value in output.values:
        if value.elementLabel in labels and value.integrationPoint is not None:
            result[(value.elementLabel, value.integrationPoint)] = data(value)
    return result


def equivalent_stress(stress):
    mean = sum(stress[:3]) / 3.0
    return math.sqrt(
        1.5
        * (
            sum((component - mean) ** 2 for component in stress[:3])
            + 2.0 * sum(component * component for component in stress[3:])
        )
    )


with open(sys.argv[2], "rb") as source:
    mesh = json.load(source)
coordinates = dict((node["label"], node["coordinates"]) for node in mesh["nodes"])
all_labels = set(coordinates.keys())
corner_labels = set(node for element in mesh["elements"] for node in element["nodes"][:8])
all_elements = set(element["label"] for element in mesh["elements"])

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    if len(odb.steps) != 1:
        raise RuntimeError("B6.0 C3D20T expects exactly one Abaqus step")
    step = list(odb.steps.values())[0]
    frame = step.frames[-1]
    temperature = nodal_values(field(frame, "NT11"), corner_labels)
    displacement = nodal_values(field(frame, "U"), all_labels)
    with open(sys.argv[3], "wb") as output:
        output.write("id,x,y,z,temperature\n")
        for label in sorted(corner_labels):
            output.write(
                "%d,%.17g,%.17g,%.17g,%.17g\n"
                % ((label,) + tuple(coordinates[label]) + (temperature[label],))
            )
    with open(sys.argv[4], "wb") as output:
        output.write("id,x,y,z,displacement_x,displacement_y,displacement_z\n")
        for label in sorted(all_labels):
            output.write(
                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n"
                % ((label,) + tuple(coordinates[label]) + tuple(displacement[label]))
            )

    stress = integration_values(field(frame, "S"), all_elements)
    point_coordinates = integration_values(field(frame, "COORD"), all_elements)
    volume = integration_values(field(frame, "IVOL"), all_elements)
    if len(stress) != 27 * len(all_elements) or len(point_coordinates) != len(stress) or len(volume) != len(stress):
        raise RuntimeError("incomplete C3D20T integration-point output")
    with open(sys.argv[5], "wb") as output:
        output.write(
            "element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress\n"
        )
        for label in sorted(all_elements):
            for point in range(1, 28):
                key = (label, point)
                output.write(
                    "%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n"
                    % (
                        label,
                        point,
                        point_coordinates[key][0],
                        point_coordinates[key][1],
                        point_coordinates[key][2],
                        volume[key],
                        equivalent_stress(stress[key]),
                    )
                )
finally:
    odb.close()
