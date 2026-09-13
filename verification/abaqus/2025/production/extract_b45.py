from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b45.py <job.odb> <nodes.csv> <contact.csv>")


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def clean(value, tolerance=1.0e-14):
    return 0.0 if abs(value) < tolerance else value


def field_with_prefix(frame, prefix):
    matches = [frame.fieldOutputs[name] for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
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
frame = odb.steps["LOAD"].frames[-1]
instance = odb.rootAssembly.instances.values()[0]
all_labels = sorted(node.label for node in instance.nodes)
node_by_label = dict((node.label, node) for node in instance.nodes)
secondary_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
primary_labels = sorted(node.label for node in instance.nodeSets["PRIMARY_CONTACT_NODES"].nodes)
displacement = nodal_values(frame, "U", all_labels)
coordinates = nodal_values(frame, "COORD", all_labels)
reaction = nodal_values(frame, "RF", primary_labels)
normal_force = nodal_values(frame, "CNORMF", secondary_labels)
gap = nodal_values(frame, "COPEN", secondary_labels)
pressure = nodal_values(frame, "CPRESS", secondary_labels)

node_output = open(sys.argv[2], "w", newline="")
node_output.write("id,x,y,z,ux,uy,uz\n")
for label in all_labels:
    point = node_by_label[label].coordinates
    value = displacement[label]
    node_output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            label - 1,
            component(point, 0),
            component(point, 1),
            component(point, 2),
            clean(component(value, 0)),
            clean(component(value, 1)),
            clean(component(value, 2)),
        )
    )
node_output.close()

contact_output = open(sys.argv[3], "w", newline="")
contact_output.write("side,id,x,y,z,force_x,force_y,force_z,gap,pressure\n")
for label in secondary_labels:
    point = coordinates[label]
    force = normal_force[label]
    contact_output.write(
        "secondary,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            label - 1,
            component(point, 0),
            component(point, 1),
            component(point, 2),
            clean(component(force, 0), 1.0e-9),
            clean(component(force, 1), 1.0e-9),
            clean(component(force, 2), 1.0e-9),
            clean(gap[label]),
            clean(pressure[label], 1.0e-9),
        )
    )
for label in primary_labels:
    point = coordinates[label]
    force = reaction[label]
    contact_output.write(
        "primary,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,0,0\n"
        % (
            label - 1,
            component(point, 0),
            component(point, 1),
            component(point, 2),
            clean(component(force, 0), 1.0e-9),
            clean(component(force, 1), 1.0e-9),
            clean(component(force, 2), 1.0e-9),
        )
    )
contact_output.close()
odb.close()
