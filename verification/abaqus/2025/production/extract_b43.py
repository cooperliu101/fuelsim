from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b43.py <job.odb> <nodes.csv> <contact.csv>")


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def field_with_prefix(frame, prefix):
    matches = []
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            matches.append(frame.fieldOutputs[name])
    if len(matches) != 1:
        raise RuntimeError("field %s has %d matches, expected one" % (prefix, len(matches)))
    return matches[0]


def nodal_values(frame, prefix, labels):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(labels):
        raise RuntimeError("%s has %d values, expected %d" % (prefix, len(result), len(labels)))
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
instance = odb.rootAssembly.instances.values()[0]
all_labels = sorted(node.label for node in instance.nodes)
node_by_label = dict((node.label, node) for node in instance.nodes)
contact_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
frames = odb.steps["DEFORM"].frames
if len(frames) != 21:
    raise RuntimeError("DEFORM has %d frames, expected the initial frame and twenty accepted increments" % len(frames))

node_output = open(sys.argv[2], "w", newline="")
contact_output = open(sys.argv[3], "w", newline="")
node_output.write("step,id,x,y,z,ux,uy,uz,rf_x,rf_y,rf_z\n")
contact_output.write(
    "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,"
    "slip_1,slip_2,gap,pressure\n"
)

for step_index in range(1, len(frames)):
    frame = frames[step_index]
    displacement = nodal_values(frame, "U", all_labels)
    reaction = nodal_values(frame, "RF", all_labels)
    for label in all_labels:
        point = node_by_label[label].coordinates
        value = displacement[label]
        force = reaction[label]
        node_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0),
                component(point, 1),
                component(point, 2),
                component(value, 0),
                component(value, 1),
                component(value, 2),
                component(force, 0),
                component(force, 1),
                component(force, 2),
            )
        )
    coordinates = nodal_values(frame, "COORD", contact_labels)
    normal = nodal_values(frame, "CNORMF", contact_labels)
    tangential = nodal_values(frame, "CSHEARF", contact_labels)
    slip_1 = nodal_values(frame, "CSLIP1", contact_labels)
    slip_2 = nodal_values(frame, "CSLIP2", contact_labels)
    gap = nodal_values(frame, "COPEN", contact_labels)
    pressure = nodal_values(frame, "CPRESS", contact_labels)
    for label in contact_labels:
        point = coordinates[label]
        force = normal[label]
        shear = tangential[label]
        contact_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0),
                component(point, 1),
                component(point, 2),
                component(force, 0),
                component(force, 1),
                component(force, 2),
                component(shear, 0),
                component(shear, 1),
                component(shear, 2),
                slip_1[label],
                slip_2[label],
                gap[label],
                pressure[label],
            )
        )

node_output.close()
contact_output.close()
odb.close()
