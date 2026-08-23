from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_h20_32.py <job.odb> <contact.csv>")


secondary_labels = [21, 24, 28, 25, 32, 36, 40, 33]


def field_with_prefix(frame, prefix):
    matches = []
    for name in frame.fieldOutputs.keys():
        if name.strip().startswith(prefix):
            matches.append(frame.fieldOutputs[name])
    if len(matches) != 1:
        raise RuntimeError("field %s has %d matches, expected one" % (prefix, len(matches)))
    return matches[0]


def nodal_values(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in secondary_labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(secondary_labels):
        raise RuntimeError(
            "%s has %d secondary values, expected %d"
            % (prefix, len(result), len(secondary_labels))
        )
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
step = odb.steps["LOAD"]
frame = step.frames[-1]
available = sorted([name.strip() for name in frame.fieldOutputs.keys()])
print("H20.32 available final-frame fields: %s" % ", ".join(available))

copen = nodal_values(frame, "COPEN")
cpress = nodal_values(frame, "CPRESS")
cnormf = nodal_values(frame, "CNORMF")
cshearf = nodal_values(frame, "CSHEARF")
cshear1 = nodal_values(frame, "CSHEAR1")
cshear2 = nodal_values(frame, "CSHEAR2")
cslip1 = nodal_values(frame, "CSLIP1")
cslip2 = nodal_values(frame, "CSLIP2")

output = open(sys.argv[2], "wb")
output.write(
    "local_node,label,copen_m,cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n,"
    "cshearf_x_n,cshearf_y_n,cshearf_z_n,cshear1_pa,cshear2_pa,cslip1_m,cslip2_m\n"
)
for local_node, label in enumerate(secondary_labels, 1):
    normal = cnormf[label]
    shear = cshearf[label]
    output.write(
        "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            local_node,
            label,
            copen[label],
            cpress[label],
            normal[0],
            normal[1],
            normal[2],
            shear[0],
            shear[1],
            shear[2],
            cshear1[label],
            cshear2[label],
            cslip1[label],
            cslip2[label],
        )
    )
output.close()
odb.close()
