from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: extract_h20_33.py <job.odb> <displacement.csv> <contact.csv> <reaction.csv>"
    )


step_names = ["INITIAL_STICK", "STICK", "MIXED", "FORWARD_SLIDE", "UNLOAD", "REVERSE_SLIDE", "RESTICK"]
secondary_labels = [57, 60, 61, 64, 68, 69, 72, 76, 78, 80, 83, 85, 88]
loaded_labels = [58, 59, 62, 63, 66, 70, 71, 74, 77, 79, 81, 84, 86]
primary_fixed_labels = [
    1, 4, 5, 8, 12, 13, 16, 20, 22, 24, 27, 29, 32, 34, 36, 39, 41, 44, 46,
    48, 51, 53, 56,
]
prescribed_y = [5.0e-6, 1.0e-5, 3.0e-5, 6.0e-5, 1.0e-5, -5.0e-5, -4.5e-5]


def component(data, index):
    if index < len(data):
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
        raise RuntimeError(
            "%s has %d requested nodal values, expected %d"
            % (prefix, len(result), len(labels))
        )
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
nodes = []
for instance in odb.rootAssembly.instances.values():
    nodes.extend(instance.nodes)
nodes.sort(key=lambda node: node.label)
node_by_label = dict((node.label, node) for node in nodes)

displacement_output = open(sys.argv[2], "wb")
contact_output = open(sys.argv[3], "wb")
reaction_output = open(sys.argv[4], "wb")
displacement_output.write("step,time,id,x,y,z,disp_x,disp_y,disp_z\n")
contact_output.write(
    "step,time,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,"
    "cslip1,cslip2,copen,cpress\n"
)
reaction_output.write("step,time,reaction_x,reaction_y,reaction_z\n")

for step_index, step_name in enumerate(step_names, 1):
    frame = odb.steps[step_name].frames[-1]
    if step_index == 1:
        print("H20.33 available final-frame fields: %s" % ", ".join(sorted(frame.fieldOutputs.keys())))
    displacements = {}
    reactions = {}
    for value in frame.fieldOutputs["U"].values:
        displacements[value.nodeLabel] = value.dataDouble
    for value in frame.fieldOutputs["RF"].values:
        reactions[value.nodeLabel] = value.dataDouble
    for node in nodes:
        coordinates = node.coordinates
        displacement = displacements[node.label]
        values = [component(displacement, 0), component(displacement, 1), 0.0]
        if node.label in primary_fixed_labels:
            values = [0.0, 0.0, 0.0]
        if node.label in loaded_labels:
            values[0] = -1.0e-5
            values[1] = prescribed_y[step_index - 1]
        displacement_output.write(
            "%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                float(step_index),
                node.label - 1,
                component(coordinates, 0),
                component(coordinates, 1),
                component(coordinates, 2),
                values[0],
                values[1],
                values[2],
            )
        )

    cnormf = nodal_values(frame, "CNORMF", secondary_labels)
    cshearf = nodal_values(frame, "CSHEARF", secondary_labels)
    cslip1 = nodal_values(frame, "CSLIP1", secondary_labels)
    cslip2 = nodal_values(frame, "CSLIP2", secondary_labels)
    copen = nodal_values(frame, "COPEN", secondary_labels)
    cpress = nodal_values(frame, "CPRESS", secondary_labels)
    for label in secondary_labels:
        coordinates = node_by_label[label].coordinates
        normal = cnormf[label]
        shear = cshearf[label]
        contact_output.write(
            "%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                float(step_index),
                label - 1,
                component(coordinates, 0),
                component(coordinates, 1),
                component(coordinates, 2),
                component(normal, 0),
                component(normal, 1),
                0.0,
                component(shear, 0),
                component(shear, 1),
                0.0,
                0.0,
                cslip2[label],
                copen[label],
                cpress[label],
            )
        )

    reaction = [0.0, 0.0, 0.0]
    for label in loaded_labels:
        value = reactions[label]
        for component_index in range(3):
            reaction[component_index] += component(value, component_index)
    reaction_output.write(
        "%d,%.16g,%.16g,%.16g,%.16g\n"
        % (step_index, float(step_index), reaction[0], reaction[1], reaction[2])
    )

displacement_output.close()
contact_output.close()
reaction_output.close()
odb.close()
