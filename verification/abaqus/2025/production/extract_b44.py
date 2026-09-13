from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5 or sys.argv[4] not in ("base", "swapped"):
    raise RuntimeError("usage: extract_b44.py <job.odb> <nodes.csv> <contact.csv> <base|swapped>")


if sys.argv[4] == "base":
    step_names = ["FULL_CLOSE", "PARTIAL_ORIGIN", "PARTIAL_LOW", "CROSS_PARTIAL", "CROSS_VERTEX", "FULL_RETURN"]
    contact_set_name = "SECONDARY_CONTACT_NODES"
else:
    step_names = ["PARTIAL_ORIGIN"]
    contact_set_name = "PRIMARY_PATCH_NODES"


def clean(value, tolerance=1.0e-14):
    return 0.0 if abs(value) < tolerance else value


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
contact_labels = sorted(node.label for node in instance.nodeSets[contact_set_name].nodes)

node_output = open(sys.argv[2], "w", newline="")
contact_output = open(sys.argv[3], "w", newline="")
node_output.write("step,id,x,y,z,ux,uy,uz\n")
contact_output.write(
    "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,"
    "slip_1,slip_2,gap,pressure\n"
)

for step_index, step_name in enumerate(step_names, 1):
    frames = odb.steps[step_name].frames
    if len(frames) != 2:
        raise RuntimeError("step %s has %d frames, expected the initial and one accepted increment" % (step_name, len(frames)))
    frame = frames[-1]
    displacement = nodal_values(frame, "U", all_labels)
    for label in all_labels:
        point = node_by_label[label].coordinates
        value = displacement[label]
        node_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0),
                component(point, 1),
                component(point, 2),
                clean(component(value, 0)),
                clean(component(value, 1)),
                clean(component(value, 2)),
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
                clean(component(force, 0), 1.0e-9),
                clean(component(force, 1), 1.0e-9),
                clean(component(force, 2), 1.0e-9),
                clean(component(shear, 0), 1.0e-9),
                clean(component(shear, 1), 1.0e-9),
                clean(component(shear, 2), 1.0e-9),
                clean(slip_1[label]),
                clean(slip_2[label]),
                clean(gap[label]),
                clean(pressure[label], 1.0e-9),
            )
        )

node_output.close()
contact_output.close()
odb.close()
