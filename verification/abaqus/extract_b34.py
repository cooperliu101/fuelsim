from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: extract_b34.py <job.odb> <nodes.csv> <contact.csv> <reaction.csv>"
    )


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
        raise RuntimeError(
            "field %s has %d matches, expected one" % (prefix, len(matches))
        )
    return matches[0]


def nodal_values(frame, prefix, labels):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(labels):
        raise RuntimeError(
            "%s has %d values, expected %d" % (prefix, len(result), len(labels))
        )
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
instance = odb.rootAssembly.instances.values()[0]
frame = odb.steps["LOAD"].frames[-1]
node_by_label = dict((node.label, node) for node in instance.nodes)
all_labels = sorted(node.label for node in instance.nodeSets["ALL_NODES"].nodes)
secondary_labels = sorted(
    node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes
)
primary_left = set(node.label for node in instance.nodeSets["PRIMARY_LEFT"].nodes)
secondary_right = set(
    node.label for node in instance.nodeSets["SECONDARY_RIGHT"].nodes
)

displacement = nodal_values(frame, "U", all_labels)
nodes_output = open(sys.argv[2], "wb")
nodes_output.write("id,x,y,z,u1,u2,u3\n")
for label in all_labels:
    point = node_by_label[label].coordinates
    value = [component(displacement[label], index) for index in range(3)]
    if label in primary_left:
        value = [0.0, 0.0, 0.0]
    elif label in secondary_right:
        value[2] = 0.0
    nodes_output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            label - 1,
            component(point, 0),
            component(point, 1),
            component(point, 2),
            value[0],
            value[1],
            value[2],
        )
    )
nodes_output.close()

cnormf = nodal_values(frame, "CNORMF", secondary_labels)
cshearf = nodal_values(frame, "CSHEARF", secondary_labels)
cslip1 = nodal_values(frame, "CSLIP1", secondary_labels)
cslip2 = nodal_values(frame, "CSLIP2", secondary_labels)
copen = nodal_values(frame, "COPEN", secondary_labels)
cpress = nodal_values(frame, "CPRESS", secondary_labels)
contact_output = open(sys.argv[3], "wb")
contact_output.write(
    "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,"
    "cslip1,cslip2,copen,cpress\n"
)
for label in secondary_labels:
    point = node_by_label[label].coordinates
    normal = cnormf[label]
    shear = cshearf[label]
    contact_output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            label - 1,
            component(point, 0),
            component(point, 1),
            component(point, 2),
            component(normal, 0),
            component(normal, 1),
            component(normal, 2),
            component(shear, 0),
            component(shear, 1),
            component(shear, 2),
            cslip1[label],
            cslip2[label],
            copen[label],
            cpress[label],
        )
    )
contact_output.close()

reaction = [0.0, 0.0, 0.0]
for value in frame.fieldOutputs["RF"].values:
    if value.nodeLabel in secondary_right:
        for index in range(3):
            reaction[index] += component(value.dataDouble, index)
reaction_output = open(sys.argv[4], "wb")
reaction_output.write("reaction_x,reaction_y,reaction_z\n")
reaction_output.write("%.16g,%.16g,%.16g\n" % tuple(reaction))
reaction_output.close()
odb.close()
