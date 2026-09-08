from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_h20_28.py <job.odb> <operator.csv> <expected-primary-nodes>")


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


def nodal_values(frame, prefix, labels):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in labels:
            data = value.dataDouble
            if value.nodeLabel in result:
                previous = result[value.nodeLabel]
                if hasattr(data, "__len__"):
                    data = tuple(previous[index] + data[index] for index in range(len(data)))
                else:
                    data = previous + data
            result[value.nodeLabel] = data
    if len(result) != len(labels):
        raise RuntimeError(
            "%s has %d requested values, expected %d"
            % (prefix, len(result), len(labels))
        )
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
instance = odb.rootAssembly.instances.values()[0]
primary_labels = sorted(node.label for node in instance.nodeSets["PRIMARY_FACE"].nodes)
secondary_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_FACE"].nodes)
input_labels = [instance.nodeSets["S%d" % local_node].nodes[0].label for local_node in range(1, 9)]
if len(primary_labels) != int(sys.argv[3]) or len(secondary_labels) != 8:
    raise RuntimeError(
        "unexpected H20.28 face-node counts: primary=%d secondary=%d"
        % (len(primary_labels), len(secondary_labels))
    )

step_names = ["BASE"]
for local_node in range(1, 9):
    step_names.append("S%d_PLUS" % local_node)
    step_names.append("S%d_MINUS" % local_node)

output = open(sys.argv[2], "wb")
output.write(
    "step,input_local_node,input_label,side,output_local_node,output_label,"
    "closure_delta_m,copen_m,cnormf_x_n,coord_y_m,coord_z_m\n"
)
coordinates = dict((node.label, node.coordinates) for node in instance.nodes)
base_displacements = nodal_values(odb.steps["BASE"].frames[-1], "U", input_labels)
for step_name in step_names:
    frame = odb.steps[step_name].frames[-1]
    secondary_force = nodal_values(frame, "CNORMF", secondary_labels)
    primary_force = nodal_values(frame, "RF", primary_labels)
    opening = nodal_values(frame, "COPEN", secondary_labels)
    input_local_node = 0
    input_label = 0
    closure_delta = 0.0
    if step_name != "BASE":
        input_local_node = int(step_name[1])
        input_label = input_labels[input_local_node - 1]
        displacement = nodal_values(frame, "U", [input_label])[input_label]
        closure_delta = base_displacements[input_label][0] - displacement[0]
        if closure_delta == 0 or (closure_delta > 0) != step_name.endswith("_PLUS"):
            raise RuntimeError("Input displacement does not match the named perturbation")
    for side, labels in (("secondary", secondary_labels), ("primary", primary_labels)):
        for output_local_node, output_label in enumerate(labels, 1):
            point = coordinates[output_label]
            output.write(
                "%s,%d,%d,%s,%d,%d,%.16g,%s,%.16g,%.16g,%.16g\n"
                % (
                    step_name,
                    input_local_node,
                    input_label,
                    side,
                    output_local_node,
                    output_label,
                    closure_delta,
                    "%.16g" % opening[output_label] if side == "secondary" else "",
                    (secondary_force if side == "secondary" else primary_force)[output_label][0],
                    point[1],
                    point[2],
                )
            )
output.close()
odb.close()
