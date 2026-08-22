from __future__ import print_function

import sys

from odbAccess import openOdb


def component(data, index):
    if index < len(data):
        return data[index]
    return 0.0


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: extract_h20_25.py <job.odb> <displacement.csv> <reaction.csv> <pressure.csv>"
    )

odb = openOdb(path=sys.argv[1], readOnly=True)
frame = odb.steps["LOAD"].frames[-1]
displacements = {}
for value in frame.fieldOutputs["U"].values:
    displacements[value.nodeLabel] = value.data

reactions = {}
for value in frame.fieldOutputs["RF"].values:
    reactions[value.nodeLabel] = value.data

all_nodes = []
for instance in odb.rootAssembly.instances.values():
    all_nodes.extend(instance.nodes)
all_nodes.sort(key=lambda node: node.label)
node_by_label = dict((node.label, node) for node in all_nodes)
nodes = [node for node in all_nodes if node.label <= 40]

output = open(sys.argv[2], "wb")
output.write("disp_x,disp_y,disp_z,id,x,y,z\n")
for node in nodes:
    displacement = displacements[node.label]
    coordinates = node.coordinates
    secondary_reference_shift = 1.0e-6 if node.label >= 21 else 0.0
    output.write(
        "%.16g,%.16g,%.16g,%d,%.16g,%.16g,%.16g\n"
        % (
            component(displacement, 0) + secondary_reference_shift,
            component(displacement, 1),
            component(displacement, 2),
            node.label - 1,
            component(coordinates, 0) - secondary_reference_shift,
            component(coordinates, 1),
            component(coordinates, 2),
        )
    )
output.close()

secondary_right = set([22, 23, 26, 27, 30, 34, 35, 38])
reaction = [0.0, 0.0, 0.0]
for label in secondary_right:
    value = reactions[label]
    for index in range(3):
        reaction[index] += component(value, index)

output = open(sys.argv[3], "wb")
output.write("time,reaction_x,reaction_y,reaction_z\n")
output.write("1,%.16g,%.16g,%.16g\n" % (reaction[0], reaction[1], reaction[2]))
output.close()

pressure_values = []
for name, field in frame.fieldOutputs.items():
    if not name.startswith("CPRESS"):
        continue
    for value in field.values:
        if not hasattr(value, "nodeLabel") or value.nodeLabel not in node_by_label:
            continue
        if value.nodeLabel < 21:
            continue
        data = value.data
        pressure = data if isinstance(data, float) else component(data, 0)
        pressure_values.append((value.nodeLabel, pressure))
if not pressure_values:
    raise RuntimeError("Abaqus ODB contains no nodal CPRESS contact output")
pressure_values.sort(key=lambda value: value[0])
output = open(sys.argv[4], "wb")
output.write("pressure,id,x,y,z\n")
for label, pressure in pressure_values:
    node = node_by_label[label]
    coordinates = node.coordinates
    output.write(
        "%.16g,%d,%.16g,%.16g,%.16g\n"
        % (
            pressure,
            label - 1,
            component(coordinates, 0) - 1.0e-6,
            component(coordinates, 1),
            component(coordinates, 2),
        )
    )
output.close()
odb.close()
