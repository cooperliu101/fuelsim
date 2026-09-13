from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_h20_40.py <job.odb> <result.csv>")


step_names = ["CLOSE", "CROSS", "STRADDLE", "RETURN"]


def component(data, index):
    if index < len(data):
        return data[index]
    return 0.0


def cleaned(value):
    if abs(value) < 1.0e-10:
        return 0.0
    return value


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
instance = odb.rootAssembly.instances.values()[0]
secondary_labels = sorted(
    node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes
)
node_by_label = dict((node.label, node) for node in instance.nodes)
output = open(sys.argv[2], "w", newline="")
output.write(
    "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,slip_2,gap,pressure\n"
)

for step_index, step_name in enumerate(step_names, 1):
    frame = odb.steps[step_name].frames[-1]
    cnormf = nodal_values(frame, "CNORMF", secondary_labels)
    cshearf = nodal_values(frame, "CSHEARF", secondary_labels)
    cslip1 = nodal_values(frame, "CSLIP1", secondary_labels)
    cslip2 = nodal_values(frame, "CSLIP2", secondary_labels)
    copen = nodal_values(frame, "COPEN", secondary_labels)
    cpress = nodal_values(frame, "CPRESS", secondary_labels)
    for label in secondary_labels:
        coordinates = node_by_label[label].coordinates
        normal = cnormf[label]
        tangent = cshearf[label]
        output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(coordinates, 0),
                component(coordinates, 1),
                component(coordinates, 2),
                cleaned(component(normal, 0)),
                cleaned(component(normal, 1)),
                cleaned(component(normal, 2)),
                cleaned(component(tangent, 0)),
                cleaned(component(tangent, 1)),
                cleaned(component(tangent, 2)),
                cslip1[label],
                cslip2[label],
                copen[label],
                cpress[label],
            )
        )

output.close()
odb.close()
