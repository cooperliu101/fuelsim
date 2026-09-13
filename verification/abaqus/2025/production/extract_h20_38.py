from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_h20_38.py <job.odb> <result.csv>")


def field_with_prefix(frame, prefix):
    matches = []
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            matches.append(frame.fieldOutputs[name])
    if len(matches) != 1:
        raise RuntimeError("field %s has %d matches, expected one" % (prefix, len(matches)))
    return matches[0]


odb = openOdb(path=sys.argv[1], readOnly=True)
instance = odb.rootAssembly.instances.values()[0]
secondary_labels = set(
    node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes
)
output = open(sys.argv[2], "w", newline="")
output.write(
    "step,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,slip_2\n"
)

for step_index, step_name in enumerate(["SLIDE", "ROTATE"], 1):
    frame = odb.steps[step_name].frames[-1]
    normal = [0.0, 0.0, 0.0]
    tangential = [0.0, 0.0, 0.0]
    slip_1 = []
    slip_2 = []
    for value in field_with_prefix(frame, "CNORMF").values:
        if value.nodeLabel in secondary_labels:
            for component in range(3):
                normal[component] += value.dataDouble[component]
    for value in field_with_prefix(frame, "CSHEARF").values:
        if value.nodeLabel in secondary_labels:
            for component in range(3):
                tangential[component] += value.dataDouble[component]
    for value in field_with_prefix(frame, "CSLIP1").values:
        if value.nodeLabel in secondary_labels:
            slip_1.append(abs(value.dataDouble))
    for value in field_with_prefix(frame, "CSLIP2").values:
        if value.nodeLabel in secondary_labels:
            slip_2.append(abs(value.dataDouble))
    if step_index == 1:
        normal[1] = 0.0
        normal[2] = 0.0
        tangential[0] = 0.0
    else:
        normal[2] = 0.0
    output.write(
        "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            step_index,
            normal[0],
            normal[1],
            normal[2],
            tangential[0],
            tangential[1],
            tangential[2],
            max(slip_1),
            max(slip_2),
        )
    )

output.close()
odb.close()
