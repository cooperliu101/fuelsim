from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: extract_h20_35.py <job.odb> <job.inp> <displacement.csv> <contact.csv>"
    )


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def input_coordinates(path):
    result = {}
    reading_nodes = False
    input_file = open(path, "r")
    for raw_line in input_file:
        line = raw_line.strip()
        if line.startswith("*"):
            reading_nodes = line.lower() == "*node"
            continue
        if not reading_nodes or not line:
            continue
        values = [value.strip() for value in line.split(",")]
        result[int(values[0])] = tuple(float(value) for value in values[1:4])
    input_file.close()
    return result


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
frame = odb.steps["LOAD"].frames[-1]
instance = odb.rootAssembly.instances.values()[0]
nodes = sorted(instance.nodes, key=lambda node: node.label)
coordinates_by_label = input_coordinates(sys.argv[2])
secondary_labels = sorted(
    node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes
)
fixed_x = set(node.label for node in instance.nodeSets["PRIMARY_INNER"].nodes)
fixed_x.update(node.label for node in instance.nodeSets["SECONDARY_THETA_UPPER"].nodes)
fixed_y = set(node.label for node in instance.nodeSets["PRIMARY_INNER"].nodes)
fixed_y.update(node.label for node in instance.nodeSets["SECONDARY_THETA_LOWER"].nodes)
fixed_z = set(node.label for node in instance.nodeSets["PRIMARY_BACK"].nodes)
fixed_z.update(node.label for node in instance.nodeSets["SECONDARY_BACK"].nodes)

displacements = {}
for value in frame.fieldOutputs["U"].values:
    displacements[value.nodeLabel] = value.dataDouble
if len(displacements) != len(nodes):
    raise RuntimeError("H20.35 displacement output does not cover every node")

output = open(sys.argv[3], "wb")
output.write("id,x,y,z,disp_x,disp_y,disp_z\n")
for node in nodes:
    coordinates = coordinates_by_label[node.label]
    displacement = displacements[node.label]
    values = [component(displacement, 0), component(displacement, 1), component(displacement, 2)]
    if node.label in fixed_x:
        values[0] = 0.0
    if node.label in fixed_y:
        values[1] = 0.0
    if node.label in fixed_z:
        values[2] = 0.0
    output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            node.label - 1,
            coordinates[0],
            coordinates[1],
            coordinates[2],
            values[0],
            values[1],
            values[2],
        )
    )
output.close()

cnormf = nodal_values(frame, "CNORMF", secondary_labels)
cshearf = nodal_values(frame, "CSHEARF", secondary_labels)
cslip1 = nodal_values(frame, "CSLIP1", secondary_labels)
cslip2 = nodal_values(frame, "CSLIP2", secondary_labels)
output = open(sys.argv[4], "wb")
output.write(
    "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2\n"
)
for label in secondary_labels:
    coordinates = coordinates_by_label[label]
    normal = cnormf[label]
    shear = cshearf[label]
    output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            label - 1,
            coordinates[0],
            coordinates[1],
            coordinates[2],
            component(normal, 0),
            component(normal, 1),
            component(normal, 2),
            component(shear, 0),
            component(shear, 1),
            component(shear, 2),
            cslip1[label],
            cslip2[label],
        )
    )
output.close()
odb.close()
