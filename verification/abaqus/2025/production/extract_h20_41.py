from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_h20_41.py <job.odb> <job.inp> <contact.csv>")


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def component(value, index):
    try:
        return value[index]
    except TypeError:
        return value if index == 0 else 0.0


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
            result[value.nodeLabel] = data(value)
    if len(result) != len(labels):
        raise RuntimeError(
            "%s has %d requested nodal values, expected %d"
            % (prefix, len(result), len(labels))
        )
    return result


def input_coordinates(path):
    result = {}
    reading_nodes = False
    source = open(path, "r")
    for raw_line in source:
        line = raw_line.strip()
        if line.startswith("*"):
            reading_nodes = line.lower() == "*node"
            continue
        if not reading_nodes or not line:
            continue
        values = [value.strip() for value in line.split(",")]
        result[int(values[0])] = tuple(float(value) for value in values[1:4])
    source.close()
    return result


coordinates = input_coordinates(sys.argv[2])
odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["LOAD"].frames[-1]
    instance = odb.rootAssembly.instances.values()[0]
    labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
    gap = nodal_values(frame, "COPEN", labels)
    pressure = nodal_values(frame, "CPRESS", labels)
    normal = nodal_values(frame, "CNORMF", labels)
    output = open(sys.argv[3], "w", newline="")
    output.write("id,x,y,z,gap,pressure,normal_x,normal_y,normal_z\n")
    for label in labels:
        output.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                label - 1,
                coordinates[label][0],
                coordinates[label][1],
                coordinates[label][2],
                component(gap[label], 0),
                component(pressure[label], 0),
                component(normal[label], 0),
                component(normal[label], 1),
                component(normal[label], 2),
            )
        )
    output.close()
finally:
    odb.close()
