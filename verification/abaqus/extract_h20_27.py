from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError(
        "usage: extract_h20_27.py <job.odb> <operator.csv> <history.csv>"
    )


secondary_labels = [57, 60, 64, 61, 68, 72, 76, 69, 78, 80, 83, 85, 88]
primary_labels = [2, 3, 7, 6, 10, 15, 18, 14, 21, 23, 25, 28, 30, 33, 35, 37, 40, 42, 45, 47, 49, 52, 54]
step_names = ["BASE"]
for local_node in range(1, len(secondary_labels) + 1):
    step_names.append("S%02d_PLUS" % local_node)
    step_names.append("S%02d_MINUS" % local_node)
step_names.extend(["OPEN", "RECLOSE"])


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


def nodal_values(frame, prefix, labels):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in labels:
            result[value.nodeLabel] = value.dataDouble
    return result


def history_value(step, prefix):
    matches = []
    for region in step.historyRegions.values():
        for name, history in region.historyOutputs.items():
            if name.strip().startswith(prefix):
                matches.append(history.data[-1][1])
    if len(matches) != 1:
        raise RuntimeError(
            "history %s has %d matches, expected one" % (prefix, len(matches))
        )
    return matches[0]


odb = openOdb(path=sys.argv[1], readOnly=True)
operator_output = open(sys.argv[2], "wb")
history_output = open(sys.argv[3], "wb")
operator_output.write(
    "step,input_local_node,input_label,side,output_local_node,output_label,"
    "closure_delta_m,copen_m,cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n\n"
)
history_output.write(
    "step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
    "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,"
    "center_x_m,center_y_m,center_z_m\n"
)

all_labels = primary_labels + secondary_labels
for step_name in step_names:
    if step_name not in odb.steps:
        raise RuntimeError("missing step %s" % step_name)
    step = odb.steps[step_name]
    frame = step.frames[-1]
    opening = nodal_values(frame, "COPEN", secondary_labels)
    pressure = nodal_values(frame, "CPRESS", all_labels)
    force = nodal_values(frame, "CNORMF", all_labels)
    if len(opening) != len(secondary_labels):
        raise RuntimeError("COPEN does not contain all secondary constraints in %s" % step_name)
    if len(pressure) != len(all_labels) or len(force) != len(all_labels):
        raise RuntimeError("contact fields do not contain all primary and secondary nodes in %s" % step_name)

    input_local_node = 0
    input_label = 0
    closure_delta = 0.0
    if step_name.startswith("S"):
        input_local_node = int(step_name[1:3])
        input_label = secondary_labels[input_local_node - 1]
        closure_delta = 1.0e-6 if step_name.endswith("_PLUS") else -1.0e-6

    for side, labels in (("secondary", secondary_labels), ("primary", primary_labels)):
        for output_local_node, output_label in enumerate(labels, 1):
            vector = force[output_label]
            copen = opening[output_label] if side == "secondary" else ""
            operator_output.write(
                "%s,%d,%d,%s,%d,%d,%.16g,%s,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    step_name,
                    input_local_node,
                    input_label,
                    side,
                    output_local_node,
                    output_label,
                    closure_delta,
                    "%.16g" % copen if copen != "" else "",
                    pressure[output_label],
                    vector[0],
                    vector[1],
                    vector[2],
                )
            )

    history_output.write(
        "%s,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            step_name,
            history_value(step, "CAREA"),
            history_value(step, "CFN1"),
            history_value(step, "CFN2"),
            history_value(step, "CFN3"),
            history_value(step, "CMN1"),
            history_value(step, "CMN2"),
            history_value(step, "CMN3"),
            history_value(step, "XN1"),
            history_value(step, "XN2"),
            history_value(step, "XN3"),
        )
    )

operator_output.close()
history_output.close()
odb.close()
