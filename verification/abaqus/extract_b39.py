from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 8:
    raise RuntimeError(
        "usage: extract_b39.py <job.odb> <job.inp> <normal_x> <normal_y> "
        "<displacement.csv> <force.csv> <reaction.csv>"
    )


def component(data, index):
    if hasattr(data, "__len__") and index < len(data):
        return data[index]
    return 0.0


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


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


normal = (float(sys.argv[3]), float(sys.argv[4]), 0.0)
odb = openOdb(path=sys.argv[1], readOnly=True)
frame = odb.steps["LOAD"].frames[-1]
instance = odb.rootAssembly.instances.values()[0]
nodes = sorted(instance.nodes, key=lambda node: node.label)
coordinates_by_label = input_coordinates(sys.argv[2])
if sorted(coordinates_by_label.keys()) != [node.label for node in nodes]:
    raise RuntimeError("B3.9 input-deck and ODB node labels differ")
secondary_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
outer_labels = set(node.label for node in instance.nodeSets["SECONDARY_OUTER"].nodes)
primary_outer_labels = set(node.label for node in instance.nodeSets["PRIMARY_OUTER"].nodes)

displacements = {}
for value in frame.fieldOutputs["U"].values:
    displacements[value.nodeLabel] = value.dataDouble
if len(displacements) != len(nodes):
    raise RuntimeError("B3.9 displacement output does not cover every node")

output = open(sys.argv[5], "wb")
output.write("normal_displacement,normal_x,normal_y,normal_z,id,x,y,z\n")
for node in nodes:
    value = displacements[node.label]
    displacement = sum(normal[index] * component(value, index) for index in range(3))
    if node.label in primary_outer_labels:
        displacement = 0.0
    coordinates = coordinates_by_label[node.label]
    output.write(
        "%.16g,%.16g,%.16g,%.16g,%d,%.16g,%.16g,%.16g\n"
        % (
            displacement,
            normal[0],
            normal[1],
            normal[2],
            node.label - 1,
            component(coordinates, 0),
            component(coordinates, 1),
            component(coordinates, 2),
        )
    )
output.close()

normal_forces = {}
for value in field_with_prefix(frame, "CNORMF").values:
    if not hasattr(value, "nodeLabel") or value.nodeLabel not in secondary_labels:
        continue
    projected = sum(normal[index] * component(value.dataDouble, index) for index in range(3))
    normal_forces[value.nodeLabel] = normal_forces.get(value.nodeLabel, 0.0) + projected
if sorted(normal_forces.keys()) != secondary_labels:
    raise RuntimeError("B3.9 CNORMF output does not cover every secondary contact node")

output = open(sys.argv[6], "wb")
output.write("normal_force,id,x,y,z\n")
for label in secondary_labels:
    coordinates = coordinates_by_label[label]
    output.write(
        "%.16g,%d,%.16g,%.16g,%.16g\n"
        % (
            normal_forces[label],
            label - 1,
            component(coordinates, 0),
            component(coordinates, 1),
            component(coordinates, 2),
        )
    )
output.close()

reaction = 0.0
for value in frame.fieldOutputs["RF"].values:
    if value.nodeLabel in outer_labels:
        reaction += sum(normal[index] * component(value.dataDouble, index) for index in range(3))
output = open(sys.argv[7], "wb")
output.write("normal_contact_resultant,normal_outer_reaction\n")
output.write("%.16g,%.16g\n" % (sum(normal_forces.values()), reaction))
output.close()
odb.close()
