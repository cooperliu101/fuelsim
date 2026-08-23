from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: extract_h20_26_geometry.py <job.odb> <case> <nodal.csv> <history.csv>"
    )


primary_labels = [2, 3, 7, 6, 10, 15, 18, 14]
secondary_labels = [21, 24, 28, 25, 32, 36, 40, 33]


def field_with_prefix(frame, prefix):
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            return frame.fieldOutputs[name]
    raise RuntimeError("missing field output beginning with %s" % prefix)


def nodal_values(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in primary_labels or value.nodeLabel in secondary_labels:
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
case_name = sys.argv[2]
nodal_output = open(sys.argv[3], "wb")
history_output = open(sys.argv[4], "wb")
nodal_output.write(
    "case,step,side,local_node,node_label,cnormf_x_n,cnormf_y_n,cnormf_z_n,"
    "copen_m,cpress_pa,coord_x_m,coord_y_m,coord_z_m\n"
)
history_output.write(
    "case,step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
    "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,"
    "center_x_m,center_y_m,center_z_m\n"
)

for step_name in odb.steps.keys():
    step = odb.steps[step_name]
    frame = step.frames[-1]
    force = nodal_values(frame, "CNORMF")
    pressure = nodal_values(frame, "CPRESS")
    opening = nodal_values(frame, "COPEN")
    coordinates = nodal_values(frame, "COORD")
    for side, labels in (("primary", primary_labels), ("secondary", secondary_labels)):
        for local_node, label in enumerate(labels, 1):
            if label not in force or label not in pressure or label not in coordinates:
                raise RuntimeError("missing contact value at node %d" % label)
            vector = force[label]
            point = coordinates[label]
            copen = opening[label] if label in opening else ""
            nodal_output.write(
                "%s,%s,%s,%d,%d,%.16g,%.16g,%.16g,%s,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    case_name,
                    step_name,
                    side,
                    local_node,
                    label,
                    vector[0],
                    vector[1],
                    vector[2],
                    "%.16g" % copen if copen != "" else "",
                    pressure[label],
                    point[0],
                    point[1],
                    point[2],
                )
            )
    history_output.write(
        "%s,%s,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            case_name,
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
