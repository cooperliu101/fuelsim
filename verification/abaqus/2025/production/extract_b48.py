from __future__ import print_function

import math
import sys

from odbAccess import openOdb


if len(sys.argv) != 6:
    raise RuntimeError(
        "usage: extract_b48.py <job.odb> <pair-a.csv> <pair-b.csv> <reaction.csv> <energy.csv>"
    )


step_names = ["STICK", "CROSS", "FORWARD", "FAR_SLIDE"]
pair_names = ["A", "B"]
friction = [0.3, 0.5]


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def contact_values(frame, prefix, labels):
    result = {}
    for name in frame.fieldOutputs.keys():
        if not name.strip().startswith(prefix):
            continue
        for value in frame.fieldOutputs[name].values:
            if value.nodeLabel in labels:
                if value.nodeLabel in result:
                    raise RuntimeError(
                        "%s has duplicate values for node %d" % (prefix, value.nodeLabel)
                    )
                result[value.nodeLabel] = value.dataDouble
    if len(result) != len(labels):
        raise RuntimeError(
            "%s has %d requested values, expected %d"
            % (prefix, len(result), len(labels))
        )
    return result


def allfd_value(step):
    matches = []
    for region in step.historyRegions.values():
        for name in region.historyOutputs.keys():
            if name.strip().startswith("ALLFD"):
                matches.append(region.historyOutputs[name])
    if len(matches) != 1:
        raise RuntimeError("ALLFD has %d matches, expected one" % len(matches))
    return matches[0].data[-1][1]


odb = openOdb(path=sys.argv[1], readOnly=True)
instance = odb.rootAssembly.instances.values()[0]
node_by_label = dict((node.label, node) for node in instance.nodes)
secondary_contact = []
secondary_all = []
for name in pair_names:
    secondary_contact.append(
        sorted(
            node.label
            for node in instance.nodeSets[
                "SECONDARY_%s_CONTACT_NODES" % name
            ].nodes
        )
    )
    secondary_all.append(
        set(node.label for node in instance.nodeSets["SECONDARY_%s_ALL" % name].nodes)
    )

contact_outputs = [open(sys.argv[2], "w", newline=""), open(sys.argv[3], "w", newline="")]
reaction_output = open(sys.argv[4], "w", newline="")
energy_output = open(sys.argv[5], "w", newline="")
for output in contact_outputs:
    output.write(
        "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,"
        "slip1,slip2,gap,pressure,state\n"
    )
reaction_output.write(
    "step,pair_a_x,pair_a_y,pair_a_z,pair_b_x,pair_b_y,pair_b_z,global_x,global_y,global_z\n"
)
energy_output.write("step,allfd\n")

for step_index, step_name in enumerate(step_names, 1):
    step = odb.steps[step_name]
    frame = step.frames[-1]
    reactions = []
    for pair in range(2):
        labels = secondary_contact[pair]
        cnormf = contact_values(frame, "CNORMF", labels)
        cshearf = contact_values(frame, "CSHEARF", labels)
        cslip1 = contact_values(frame, "CSLIP1", labels)
        cslip2 = contact_values(frame, "CSLIP2", labels)
        copen = contact_values(frame, "COPEN", labels)
        cpress = contact_values(frame, "CPRESS", labels)
        for label in labels:
            point = node_by_label[label].coordinates
            normal = cnormf[label]
            tangent = cshearf[label]
            normal_magnitude = math.sqrt(
                sum(component(normal, index) ** 2 for index in range(3))
            )
            tangent_magnitude = math.sqrt(
                sum(component(tangent, index) ** 2 for index in range(3))
            )
            if cpress[label] <= 0.0:
                state = 0
            elif tangent_magnitude >= 0.999999 * friction[pair] * normal_magnitude:
                state = 2
            else:
                state = 1
            contact_outputs[pair].write(
                "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%d\n"
                % (
                    step_index,
                    label - 1,
                    component(point, 0),
                    component(point, 1),
                    component(point, 2),
                    component(normal, 0),
                    component(normal, 1),
                    component(normal, 2),
                    component(tangent, 0),
                    component(tangent, 1),
                    component(tangent, 2),
                    cslip1[label],
                    cslip2[label],
                    copen[label],
                    cpress[label],
                    state,
                )
            )
        reaction = [0.0, 0.0, 0.0]
        for value in frame.fieldOutputs["RF"].values:
            if value.nodeLabel in secondary_all[pair]:
                for component_index in range(3):
                    reaction[component_index] += component(
                        value.dataDouble, component_index
                    )
        reactions.append(reaction)
    global_reaction = [
        reactions[0][index] + reactions[1][index] for index in range(3)
    ]
    reaction_output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            step_index,
            reactions[0][0],
            reactions[0][1],
            reactions[0][2],
            reactions[1][0],
            reactions[1][1],
            reactions[1][2],
            global_reaction[0],
            global_reaction[1],
            global_reaction[2],
        )
    )
    energy_output.write("%d,%.16g\n" % (step_index, allfd_value(step)))

for output in contact_outputs:
    output.close()
reaction_output.close()
energy_output.close()
odb.close()
