from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b40.py <job.odb> <contact.csv> <reaction.csv>")


step_names = ["STICK", "SLIDE_A", "SLIDE_B", "SLIDE_C"]


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
node_by_label = dict((node.label, node) for node in instance.nodes)
secondary_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
secondary_all = set(node.label for node in instance.nodeSets["SECONDARY_ALL"].nodes)

contact_output = open(sys.argv[2], "w", newline="")
reaction_output = open(sys.argv[3], "w", newline="")
contact_output.write(
    "step,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,"
    "cslip1,cslip2,copen,cpress\n"
)
reaction_output.write("step,reaction_x,reaction_y,reaction_z\n")

for step_index, step_name in enumerate(step_names, 1):
    frame = odb.steps[step_name].frames[-1]
    cnormf = nodal_values(frame, "CNORMF", secondary_labels)
    cshearf = nodal_values(frame, "CSHEARF", secondary_labels)
    cslip1 = nodal_values(frame, "CSLIP1", secondary_labels)
    cslip2 = nodal_values(frame, "CSLIP2", secondary_labels)
    copen = nodal_values(frame, "COPEN", secondary_labels)
    cpress = nodal_values(frame, "CPRESS", secondary_labels)
    for label in secondary_labels:
        point = node_by_label[label].coordinates
        normal = cnormf[label]
        shear = cshearf[label]
        contact_output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
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
    reaction = [0.0, 0.0, 0.0]
    for value in frame.fieldOutputs["RF"].values:
        if value.nodeLabel in secondary_all:
            for index in range(3):
                reaction[index] += component(value.dataDouble, index)
    reaction_output.write(
        "%d,%.16g,%.16g,%.16g\n" % (step_index, reaction[0], reaction[1], reaction[2])
    )

contact_output.close()
reaction_output.close()
odb.close()
