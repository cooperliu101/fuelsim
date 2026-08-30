from __future__ import print_function

import json
import math
import sys

from odbAccess import openOdb


if len(sys.argv) != 6:
    raise RuntimeError(
        "usage: extract_b546.py <job.odb> <mesh.json> <nodal.csv> <contact.csv> <clad-state.csv>"
    )


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def exact_field(frame, name):
    if name not in frame.fieldOutputs:
        raise RuntimeError(
            "missing Abaqus field %s; available fields are %s"
            % (name, sorted(frame.fieldOutputs.keys()))
        )
    return frame.fieldOutputs[name]


def contact_field(frame, prefix):
    matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix) and "/" in name]
    if len(matches) != 1:
        raise RuntimeError("expected one contact field beginning with %s, got %s" % (prefix, matches))
    return frame.fieldOutputs[matches[0]]


def nodal_values(field, labels):
    result = {}
    for value in field.values:
        if value.nodeLabel in labels:
            if value.nodeLabel in result:
                raise RuntimeError("duplicate nodal value for node %d" % value.nodeLabel)
            result[value.nodeLabel] = data(value)
    if len(result) != len(labels):
        raise RuntimeError("nodal field has %d requested values, expected %d" % (len(result), len(labels)))
    return result


def integration_values(field, labels):
    result = {}
    for value in field.values:
        if value.elementLabel in labels and value.integrationPoint:
            result[(value.elementLabel, value.integrationPoint)] = data(value)
    return result


def determinant(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def reference_weights(points):
    gauss = 1.0 / math.sqrt(3.0)
    signs = (
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    )
    result = []
    for zeta in (-gauss, gauss):
        for eta in (-gauss, gauss):
            for xi in (-gauss, gauss):
                jacobian = [[0.0 for unused_natural in range(3)] for unused_physical in range(3)]
                for node, (sx, sy, sz) in enumerate(signs):
                    derivative = (
                        0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
                        0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
                        0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta),
                    )
                    for physical in range(3):
                        for natural in range(3):
                            jacobian[physical][natural] += points[node][physical] * derivative[natural]
                weight = determinant(jacobian)
                if weight <= 0.0:
                    raise RuntimeError("mesh manifest contains a nonpositive HEX8 reference Jacobian")
                result.append(weight)
    return result


def equivalent_stress(stress):
    mean = (stress[0] + stress[1] + stress[2]) / 3.0
    return math.sqrt(
        1.5
        * (
            (stress[0] - mean) ** 2
            + (stress[1] - mean) ** 2
            + (stress[2] - mean) ** 2
            + 2.0 * (stress[3] ** 2 + stress[4] ** 2 + stress[5] ** 2)
        )
    )


with open(sys.argv[2], "rb") as source:
    mesh = json.load(source)
coordinates = dict((node["label"], node["coordinates"]) for node in mesh["nodes"])
all_labels = set(coordinates.keys())
node_sets = dict((node_set["name"], node_set["nodes"]) for node_set in mesh["node_sets"])
element_records = dict((element["label"], element) for element in mesh["elements"])
clad_labels = set(
    element["label"]
    for block in mesh["blocks"]
    if block["name"] == "CLAD"
    for element in block["elements"]
)
if len(all_labels) != 1617 or len(clad_labels) != 512 or len(node_sets.get("FUEL_OUTER", [])) != 144:
    raise RuntimeError("unexpected B5.46 mesh manifest counts")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    step = odb.steps["PATH"]
    frame = step.frames[-1]
    if abs(frame.frameValue - 1.0) > 1.0e-12:
        raise RuntimeError("B5.46 final Abaqus frame is not at one second")
    temperature = nodal_values(exact_field(frame, "NT11"), all_labels)
    displacement = nodal_values(exact_field(frame, "U"), all_labels)
    fuel_bottom = set(node_sets["FUEL_BOTTOM"])
    clad_bottom = set(node_sets["CLAD_BOTTOM"])
    for label in fuel_bottom | clad_bottom:
        value = list(displacement[label])
        value[0] = 0.0
        value[1] = 0.0
        if label in clad_bottom:
            value[2] = 0.0
        displacement[label] = tuple(value)

    nodal = open(sys.argv[3], "wb")
    nodal.write("T,disp_x,disp_y,disp_z,id,x,y,z\n")
    for label in sorted(all_labels):
        nodal.write(
            "%.16g,%.16g,%.16g,%.16g,%d,%.17g,%.17g,%.17g\n"
            % ((temperature[label],) + tuple(displacement[label]) + (label,) + tuple(coordinates[label]))
        )
    nodal.close()

    contact_labels = set(node_sets["FUEL_OUTER"])
    pressure = nodal_values(contact_field(frame, "CPRESS"), contact_labels)
    contact = open(sys.argv[4], "wb")
    contact.write("contact_pressure,id,x,y,z\n")
    for label in sorted(contact_labels):
        contact.write(
            "%.16g,%d,%.17g,%.17g,%.17g\n" % ((pressure[label], label) + tuple(coordinates[label]))
        )
    contact.close()

    stress = integration_values(exact_field(frame, "S"), clad_labels)
    plastic = integration_values(exact_field(frame, "PEEQ"), clad_labels)
    creep = integration_values(exact_field(frame, "CEEQ"), clad_labels)
    clad = open(sys.argv[5], "wb")
    clad.write("effective_creep_strain,effective_plastic_strain,id,vonmises_stress,x,y,z\n")
    for label in sorted(clad_labels):
        nodes = element_records[label]["nodes"]
        points = [coordinates[node] for node in nodes]
        weights = reference_weights(points)
        total = sum(weights)
        averages = [0.0, 0.0, 0.0]
        for point in range(1, 9):
            key = (label, point)
            if key not in stress or key not in plastic or key not in creep:
                raise RuntimeError("missing clad integration field for element %d point %d" % (label, point))
            weight = weights[point - 1]
            averages[0] += weight * equivalent_stress(stress[key]) / total
            averages[1] += weight * plastic[key] / total
            averages[2] += weight * creep[key] / total
        centroid = [sum(point[component] for point in points) / 8.0 for component in range(3)]
        clad.write(
            "%.16g,%.16g,%d,%.16g,%.17g,%.17g,%.17g\n"
            % (averages[2], averages[1], label, averages[0], centroid[0], centroid[1], centroid[2])
        )
    clad.close()
finally:
    odb.close()
