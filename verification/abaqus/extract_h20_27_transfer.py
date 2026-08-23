from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError(
        "usage: extract_h20_27_transfer.py <job.odb> <nodal.csv> <history.csv>"
    )


secondary_labels = [57, 60, 64, 61, 68, 72, 76, 69]
primary_labels = [2, 3, 7, 6, 10, 15, 18, 14, 21, 23, 25, 28, 30, 33, 35, 37, 40, 42, 45, 47, 49, 52, 54]
step_names = ["BASE", "SHIFT", "OPEN", "RECLOSE", "RETURN"]


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
nodal_output = open(sys.argv[2], "wb")
history_output = open(sys.argv[3], "wb")
nodal_output.write(
    "step,side,local_node,node_label,copen_m,cpress_pa,"
    "cnormf_x_n,cnormf_y_n,cnormf_z_n\n"
)
history_output.write(
    "step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
    "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,"
    "center_x_m,center_y_m,center_z_m\n"
)

all_labels = primary_labels + secondary_labels
for step_name in step_names:
    step = odb.steps[step_name]
    frame = step.frames[-1]
    opening = nodal_values(frame, "COPEN", secondary_labels)
    pressure = nodal_values(frame, "CPRESS", all_labels)
    force = nodal_values(frame, "CNORMF", all_labels)
    if len(opening) != len(secondary_labels):
        raise RuntimeError("COPEN does not contain all secondary constraints in %s" % step_name)
    if len(pressure) != len(all_labels) or len(force) != len(all_labels):
        raise RuntimeError("contact fields do not contain all nodes in %s" % step_name)
    for side, labels in (("secondary", secondary_labels), ("primary", primary_labels)):
        for local_node, label in enumerate(labels, 1):
            vector = force[label]
            copen = opening[label] if side == "secondary" else ""
            nodal_output.write(
                "%s,%s,%d,%d,%s,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    step_name,
                    side,
                    local_node,
                    label,
                    "%.16g" % copen if copen != "" else "",
                    pressure[label],
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

nodal_output.close()
history_output.close()
odb.close()
