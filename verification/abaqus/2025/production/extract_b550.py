from __future__ import print_function

import json
import math
import sys

from odbAccess import openOdb


if len(sys.argv) != 7:
    raise RuntimeError(
        "usage: extract_b550.py <job.odb> <mesh.json> <temperature.csv> <displacement.csv> <material.csv> <contact.csv>"
    )


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def clean(value):
    try:
        len(value)
    except TypeError:
        return 0.0 if abs(value) < 1.0e-20 else value
    else:
        return tuple(0.0 if abs(component) < 1.0e-20 else component for component in value)


def exact_field(frame, name):
    if name not in frame.fieldOutputs:
        raise RuntimeError("missing Abaqus field %s" % name)
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
            result[value.nodeLabel] = clean(data(value))
    if len(result) != len(labels):
        raise RuntimeError("nodal field has %d requested values, expected %d" % (len(result), len(labels)))
    return result


def integration_values(field, labels):
    result = {}
    for value in field.values:
        if value.elementLabel in labels and value.integrationPoint:
            result[(value.elementLabel, value.integrationPoint)] = data(value)
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
corner_labels = set(node for element in mesh["elements"] for node in element["nodes"][:8])
all_elements = set(element["label"] for element in mesh["elements"])
element_records = dict((element["label"], element) for element in mesh["elements"])
face_nodes = (
    (0, 1, 5, 4, 8, 13, 16, 12),
    (1, 2, 6, 5, 9, 14, 17, 13),
    (2, 3, 7, 6, 10, 15, 18, 14),
    (3, 0, 4, 7, 11, 12, 19, 15),
    (0, 3, 2, 1, 11, 10, 9, 8),
    (4, 5, 6, 7, 16, 17, 18, 19),
)
secondary_faces = next(
    side_set["faces"] for side_set in mesh["side_sets"] if side_set["name"] == "SECONDARY_CONTACT"
)
secondary_contact_labels = set()
for face in secondary_faces:
    nodes = element_records[face["element"]]["nodes"]
    secondary_contact_labels.update(nodes[index] for index in face_nodes[face["exodus_side"] - 1])
if len(all_labels) != 112 or len(corner_labels) != 40 or len(all_elements) != 8 or len(secondary_contact_labels) != 8:
    raise RuntimeError("unexpected B5.50 manifest topology")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["PATH"].frames[-1]
    print("B5.50 final frame time=%.17g" % frame.frameValue)
    if abs(frame.frameValue - 0.4) > 1.0e-6:
        raise RuntimeError("B5.50 final frame is not at 0.4 seconds")
    temperature = nodal_values(exact_field(frame, "NT11"), corner_labels)
    displacement = nodal_values(exact_field(frame, "U"), all_labels)
    with open(sys.argv[3], "w", newline="") as output:
        output.write("id,x,y,z,temperature\n")
        for label in sorted(corner_labels):
            output.write("%d,%.17g,%.17g,%.17g,%.16g\n" % ((label,) + tuple(coordinates[label]) + (temperature[label],)))
    with open(sys.argv[4], "w", newline="") as output:
        output.write("id,x,y,z,displacement_x,displacement_y,displacement_z\n")
        for label in sorted(all_labels):
            output.write("%d,%.17g,%.17g,%.17g,%.16g,%.16g,%.16g\n" % ((label,) + tuple(coordinates[label]) + tuple(displacement[label])))

    stress = integration_values(exact_field(frame, "S"), all_elements)
    plastic = integration_values(exact_field(frame, "PEEQ"), all_elements)
    creep = integration_values(exact_field(frame, "CEEQ"), all_elements)
    point_coordinates = integration_values(exact_field(frame, "COORD"), all_elements)
    volume = integration_values(exact_field(frame, "IVOL"), all_elements)
    print(
        "B5.50 integration rows S=%d PEEQ=%d CEEQ=%d COORD=%d IVOL=%d sample=%s"
        % (len(stress), len(plastic), len(creep), len(point_coordinates), len(volume), sorted(stress.keys())[:8])
    )
    with open(sys.argv[5], "w", newline="") as output:
        output.write("element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress,effective_plastic_strain,effective_creep_strain\n")
        for label in sorted(all_elements):
            for point in range(1, 28):
                key = (label, point)
                if key not in stress or key not in point_coordinates or key not in volume:
                    raise RuntimeError("missing B5.50 integration value for %s" % (key,))
                output.write("%d,%d,%.17g,%.17g,%.17g,%.17g,%.16g,%.16g,%.16g\n" % (label, point, point_coordinates[key][0], point_coordinates[key][1], point_coordinates[key][2], volume[key], equivalent_stress(stress[key]), plastic.get(key, 0.0), creep.get(key, 0.0)))

    gap = nodal_values(contact_field(frame, "COPEN"), secondary_contact_labels)
    pressure = nodal_values(contact_field(frame, "CPRESS"), secondary_contact_labels)
    normal = nodal_values(contact_field(frame, "CNORMF"), secondary_contact_labels)
    shear = nodal_values(contact_field(frame, "CSHEARF"), secondary_contact_labels)
    with open(sys.argv[6], "w", newline="") as output:
        output.write("id,x,y,z,gap,pressure,normal_x,normal_y,normal_z,shear_x,shear_y,shear_z\n")
        for label in sorted(secondary_contact_labels):
            output.write("%d,%.17g,%.17g,%.17g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n" % ((label,) + tuple(coordinates[label]) + (gap[label], pressure[label]) + tuple(normal[label]) + tuple(shear[label])))
finally:
    odb.close()
