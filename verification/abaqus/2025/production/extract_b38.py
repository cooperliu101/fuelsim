from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b38.py <job.odb> <operator.csv>")


secondary_labels = [9, 12, 16, 13]
step_names = ["BASE"]
for local_node in range(1, 5):
    step_names.append("S%d_PLUS" % local_node)
    step_names.append("S%d_MINUS" % local_node)


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


def nodal_values(frame, prefix):
    result = {}
    field = field_with_prefix(frame, prefix)
    for value in field.values:
        if value.nodeLabel in secondary_labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(secondary_labels):
        raise RuntimeError(
            "%s has %d secondary values, expected %d"
            % (prefix, len(result), len(secondary_labels))
        )
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
output = open(sys.argv[2], "w", newline="")
output.write(
    "step,input_local_node,input_label,output_local_node,output_label,"
    "closure_delta_m,copen_m,cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n\n"
)

for step_name in step_names:
    if step_name not in odb.steps:
        raise RuntimeError("missing step %s" % step_name)
    frame = odb.steps[step_name].frames[-1]
    copen = nodal_values(frame, "COPEN")
    cpress = nodal_values(frame, "CPRESS")
    cnormf = nodal_values(frame, "CNORMF")
    input_local_node = 0
    input_label = 0
    closure_delta = 0.0
    if step_name != "BASE":
        input_local_node = int(step_name[1])
        input_label = secondary_labels[input_local_node - 1]
        closure_delta = 1.0e-6 if step_name.endswith("_PLUS") else -1.0e-6
    for output_local_node, output_label in enumerate(secondary_labels, 1):
        force = cnormf[output_label]
        output.write(
            "%s,%d,%d,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_name,
                input_local_node,
                input_label,
                output_local_node,
                output_label,
                closure_delta,
                copen[output_label],
                cpress[output_label],
                force[0],
                force[1],
                force[2],
            )
        )

output.close()
odb.close()
