from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_h20_31.py <job.odb> <operator.csv>")


secondary_labels = [90, 93, 97, 94, 101, 105, 109, 102]
step_names = ["BASE"]
for local_node in range(1, 9):
    for component in ("X", "Y", "Z"):
        step_names.append("S%d_%s_PLUS" % (local_node, component))
        step_names.append("S%d_%s_MINUS" % (local_node, component))


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


def nodal_values(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in secondary_labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(secondary_labels):
        raise RuntimeError("%s has incomplete H20.31 secondary output" % prefix)
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
output = open(sys.argv[2], "wb")
output.write("step,input_local_node,input_component,output_local_node,displacement_delta_m,copen_m,cpress_pa,"
             "cnormf_x_n,cnormf_y_n,cnormf_z_n\n")
for step_name in step_names:
    frame = odb.steps[step_name].frames[-1]
    copen = nodal_values(frame, "COPEN")
    cpress = nodal_values(frame, "CPRESS")
    cnormf = nodal_values(frame, "CNORMF")
    input_node = 0 if step_name == "BASE" else int(step_name[1])
    component = "NONE" if step_name == "BASE" else step_name.split("_")[1]
    delta = 0.0 if step_name == "BASE" else (1.0e-6 if step_name.endswith("_PLUS") else -1.0e-6)
    for output_node, label in enumerate(secondary_labels, 1):
        force = cnormf[label]
        output.write(
            "%s,%d,%s,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (step_name, input_node, component, output_node, delta, copen[label], cpress[label],
               force[0], force[1], force[2])
        )
output.close()
odb.close()
