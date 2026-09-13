from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b46.py <job.odb> <nodes.csv> <contact.csv>")


step_names = ["CLOSE", "SLIDE", "OPEN", "OPEN_CROSS", "RECONTACT", "SLIDE_OUT"]


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def clean(value, tolerance=1.0e-10):
    return 0.0 if abs(value) < tolerance else value


def field_with_prefix(frame, prefix):
    matches = [
        frame.fieldOutputs[name]
        for name in frame.fieldOutputs.keys()
        if name.strip().startswith(prefix)
    ]
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
secondary_labels = sorted(
    node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes
)
node_by_label = dict((node.label, node) for node in instance.nodes)
node_output = open(sys.argv[2], "w", newline="")
contact_output = open(sys.argv[3], "w", newline="")
node_output.write("step,id,x,y,z,ux,uy,uz,rfx,rfy,rfz\n")
contact_output.write("step,id,x,y,z,normal_x,normal_y,normal_z,gap,pressure\n")

for step_index, step_name in enumerate(step_names, 1):
    frame = odb.steps[step_name].frames[-1]
    displacement = nodal_values(frame, "U", all_labels)
    reaction = nodal_values(frame, "RF", all_labels)
    normal_force = nodal_values(frame, "CNORMF", secondary_labels)
    gap = nodal_values(frame, "COPEN", secondary_labels)
    pressure = nodal_values(frame, "CPRESS", secondary_labels)
    for label in all_labels:
        point = node_by_label[label].coordinates
        u = displacement[label]
        force = reaction[label]
        node_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0),
                component(point, 1),
                component(point, 2),
                clean(component(u, 0)),
                clean(component(u, 1)),
                clean(component(u, 2)),
                clean(component(force, 0), 1.0e-7),
                clean(component(force, 1), 1.0e-7),
                clean(component(force, 2), 1.0e-7),
            )
        )
    for label in secondary_labels:
        point = node_by_label[label].coordinates
        u = displacement[label]
        force = normal_force[label]
        contact_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0) + component(u, 0),
                component(point, 1) + component(u, 1),
                component(point, 2) + component(u, 2),
                clean(component(force, 0)),
                clean(component(force, 1)),
                clean(component(force, 2)),
                clean(gap[label]),
                clean(pressure[label]),
            )
        )

node_output.close()
contact_output.close()
odb.close()
