from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b42.py <job.odb> <result.csv>")


step_names = ["STICK", "SLIDE"] + ["ROTATE_%02d" % index for index in range(1, 21)] + ["RESTICK"]


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
secondary_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)

output = open(sys.argv[2], "wb")
output.write(
    "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,"
    "slip_1,slip_2,gap,pressure\n"
)

for step_index, step_name in enumerate(step_names, 1):
    frames = odb.steps[step_name].frames
    if len(frames) != 2:
        raise RuntimeError("step %s has %d frames, expected the initial and one accepted increment" % (step_name, len(frames)))
    frame = frames[-1]
    coordinates = nodal_values(frame, "COORD", secondary_labels)
    normal = nodal_values(frame, "CNORMF", secondary_labels)
    tangential = nodal_values(frame, "CSHEARF", secondary_labels)
    slip_1 = nodal_values(frame, "CSLIP1", secondary_labels)
    slip_2 = nodal_values(frame, "CSLIP2", secondary_labels)
    gap = nodal_values(frame, "COPEN", secondary_labels)
    pressure = nodal_values(frame, "CPRESS", secondary_labels)
    for label in secondary_labels:
        point = coordinates[label]
        normal_force = normal[label]
        tangential_force = tangential[label]
        output.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_index,
                label - 1,
                component(point, 0),
                component(point, 1),
                component(point, 2),
                component(normal_force, 0),
                component(normal_force, 1),
                component(normal_force, 2),
                component(tangential_force, 0),
                component(tangential_force, 1),
                component(tangential_force, 2),
                slip_1[label],
                slip_2[label],
                gap[label],
                pressure[label],
            )
        )

output.close()
odb.close()
