from __future__ import print_function

import sys

from odbAccess import openOdb


def component(data, index):
    if index < len(data):
        return data[index]
    return 0.0


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_h20_21.py <job.odb> <displacement.csv> <reaction.csv>")

odb = openOdb(path=sys.argv[1], readOnly=True)
step = odb.steps["LOAD"]
frame = step.frames[-1]
displacements = {}
for value in frame.fieldOutputs["U"].values:
    displacements[value.nodeLabel] = value.data

reactions = {}
for value in frame.fieldOutputs["RF"].values:
    reactions[value.nodeLabel] = value.data

nodes = []
for instance in odb.rootAssembly.instances.values():
    nodes.extend(instance.nodes)
nodes.sort(key=lambda node: node.label)

output = open(sys.argv[2], "wb")
output.write("disp_x,disp_y,disp_z,id,x,y,z\n")
for node in nodes:
    displacement = displacements[node.label]
    coordinates = node.coordinates
    secondary_reference_shift = 1.0e-6 if node.label >= 21 else 0.0
    output.write("%.16g,%.16g,%.16g,%d,%.16g,%.16g,%.16g\n" %
                 (component(displacement, 0) + secondary_reference_shift,
                  component(displacement, 1), component(displacement, 2), node.label - 1,
                  component(coordinates, 0) - secondary_reference_shift,
                  component(coordinates, 1), component(coordinates, 2)))
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
odb.close()
