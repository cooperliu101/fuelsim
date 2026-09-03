from __future__ import print_function

import json
import math
import sys

from odbAccess import openOdb


if len(sys.argv) not in (6, 7):
    raise RuntimeError(
        "usage: extract_b548.py <job.odb> <mesh.json> <nodal.csv> <contact.csv> <clad-points.csv> [--friction]"
    )
friction_output = len(sys.argv) == 7 and sys.argv[6] == "--friction"
if len(sys.argv) == 7 and not friction_output:
    raise RuntimeError("unknown extract_b548.py option %s" % sys.argv[6])


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
            key = (value.elementLabel, value.integrationPoint)
            if key in result:
                raise RuntimeError("duplicate integration value for element %d point %d" % key)
            result[key] = data(value)
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
element_records = dict((element["label"], element) for element in mesh["elements"])
face_nodes = (
    (0, 1, 5, 4, 8, 13, 16, 12),
    (1, 2, 6, 5, 9, 14, 17, 13),
    (2, 3, 7, 6, 10, 15, 18, 14),
    (3, 0, 4, 7, 11, 12, 19, 15),
    (0, 3, 2, 1, 11, 10, 9, 8),
    (4, 5, 6, 7, 16, 17, 18, 19),
)
fuel_outer_faces = next(side_set["faces"] for side_set in mesh["side_sets"] if side_set["name"] == "FUEL_OUTER")
fuel_outer_labels = set()
for face in fuel_outer_faces:
    nodes = element_records[face["element"]]["nodes"]
    fuel_outer_labels.update(nodes[index] for index in face_nodes[face["exodus_side"] - 1])
clad_labels = set(
    element["label"]
    for block in mesh["blocks"]
    if block["name"] == "CLAD"
    for element in block["elements"]
)
if len(all_labels) != 5969 or len(clad_labels) != 512 or len(fuel_outer_labels) != 416:
    raise RuntimeError("unexpected B5.48 mesh manifest counts")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["PATH"].frames[-1]
    if abs(frame.frameValue - 1.0) > 1.0e-12:
        raise RuntimeError("B5.48 final Abaqus frame is not at one second")
    print("B5.48 contact fields=%s" % sorted(name for name in frame.fieldOutputs.keys() if "/" in name))

    temperature = nodal_values(exact_field(frame, "NT11"), all_labels)
    displacement = nodal_values(exact_field(frame, "U"), all_labels)
    nodal = open(sys.argv[3], "wb")
    nodal.write("T,disp_x,disp_y,disp_z,id,x,y,z\n")
    for label in sorted(all_labels):
        nodal.write(
            "%.16g,%.16g,%.16g,%.16g,%d,%.17g,%.17g,%.17g\n"
            % ((temperature[label],) + tuple(displacement[label]) + (label,) + tuple(coordinates[label]))
        )
    nodal.close()

    pressure_field = contact_field(frame, "CPRESS")
    opening = nodal_values(contact_field(frame, "COPEN"), fuel_outer_labels)
    normal_force = nodal_values(contact_field(frame, "CNORMF"), fuel_outer_labels)
    shear_force = nodal_values(contact_field(frame, "CSHEARF"), fuel_outer_labels)
    shear_1 = nodal_values(contact_field(frame, "CSHEAR1"), fuel_outer_labels)
    shear_2 = nodal_values(contact_field(frame, "CSHEAR2"), fuel_outer_labels)
    if friction_output:
        slip_1 = nodal_values(contact_field(frame, "CSLIP1"), fuel_outer_labels)
        slip_2 = nodal_values(contact_field(frame, "CSLIP2"), fuel_outer_labels)
        tangent_1 = nodal_values(contact_field(frame, "CTANDIR1"), fuel_outer_labels)
        tangent_2 = nodal_values(contact_field(frame, "CTANDIR2"), fuel_outer_labels)
    heat_flow = nodal_values(contact_field(frame, "HFL"), fuel_outer_labels)
    current_coordinates = nodal_values(exact_field(frame, "COORD"), all_labels)
    pressure = {}
    generated = 0
    print(
        "B5.48 raw CPRESS rows=%d unknown structural labels=%d maximum label=%d"
        % (
            len(pressure_field.values),
            sum(1 for value in pressure_field.values if value.nodeLabel not in coordinates),
            max(value.nodeLabel for value in pressure_field.values),
        )
    )
    for value in pressure_field.values:
        label = value.nodeLabel
        if label not in fuel_outer_labels:
            continue
        if label in pressure:
            raise RuntimeError("duplicate contact pressure for node %d" % label)
        pressure[label] = data(value)
        if label not in coordinates:
            generated += 1
    contact = open(sys.argv[4], "wb")
    contact.write("contact_pressure,contact_opening,shear_1,shear_2,normal_force_x,normal_force_y,normal_force_z,shear_force_x,shear_force_y,shear_force_z,heat_flow,id,x,y,z,current_x,current_y,current_z,generated")
    if friction_output:
        contact.write(",slip_1,slip_2,tangent_1_x,tangent_1_y,tangent_1_z,tangent_2_x,tangent_2_y,tangent_2_z")
    contact.write("\n")
    instance = odb.rootAssembly.instances["PART-1-1"]
    instance_coordinates = dict((node.label, node.coordinates) for node in instance.nodes)
    for label in sorted(pressure.keys()):
        reference = coordinates.get(label, instance_coordinates[label])
        current = current_coordinates.get(label, reference)
        values = (
            (pressure[label], opening[label], shear_1[label], shear_2[label])
            + tuple(normal_force[label])
            + tuple(shear_force[label])
            + (heat_flow[label], label)
            + tuple(reference)
            + tuple(current)
            + (0 if label in coordinates else 1,)
        )
        line = "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d" % values
        if friction_output:
            line += ",%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g" % (
                (slip_1[label], slip_2[label]) + tuple(tangent_1[label]) + tuple(tangent_2[label])
            )
        contact.write(line + "\n")
    contact.close()
    if len(pressure) != len(fuel_outer_labels):
        raise RuntimeError("B5.48 contact output has %d secondary nodes, expected %d" % (len(pressure), len(fuel_outer_labels)))
    print("B5.48 secondary contact output nodes=%d generated=%d" % (len(pressure), generated))

    stress = integration_values(exact_field(frame, "S"), clad_labels)
    plastic = integration_values(exact_field(frame, "PEEQ"), clad_labels)
    creep = integration_values(exact_field(frame, "CEEQ"), clad_labels)
    point_coordinates = integration_values(exact_field(frame, "COORD"), clad_labels)
    volume = integration_values(exact_field(frame, "IVOL"), clad_labels)
    print(
        "B5.48 clad integration rows S=%d PEEQ=%d CEEQ=%d COORD=%d IVOL=%d"
        % (len(stress), len(plastic), len(creep), len(point_coordinates), len(volume))
    )
    clad = open(sys.argv[5], "wb")
    clad.write("element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress,effective_plastic_strain,effective_creep_strain\n")
    for label in sorted(clad_labels):
        for point in range(1, 28):
            key = (label, point)
            if key not in stress or key not in plastic or key not in creep or key not in point_coordinates or key not in volume:
                raise RuntimeError("missing clad integration field for element %d point %d" % key)
            clad.write(
                "%d,%d,%.17g,%.17g,%.17g,%.17g,%.16g,%.16g,%.16g\n"
                % (
                    label,
                    point,
                    point_coordinates[key][0],
                    point_coordinates[key][1],
                    point_coordinates[key][2],
                    volume[key],
                    equivalent_stress(stress[key]),
                    plastic[key],
                    creep[key],
                )
            )
    clad.close()
finally:
    odb.close()
