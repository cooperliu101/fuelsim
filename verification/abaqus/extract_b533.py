from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b533.py <job.odb> <nodal.csv>")


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["FINITE_HOURGLASS"].frames[-1]
    reaction = dict((value.nodeLabel, double_data(value)) for value in frame.fieldOutputs["RF"].values)
    names = []
    for geometry in ("regular", "warped"):
        for state in (
            "general",
            "hourglass_vector",
            "hourglass_x_small",
            "hourglass_x_large",
            "hourglass_y_large",
            "hourglass_z_large",
        ):
            for control in ("default", "weak"):
                names.append("%s_%s_%s" % (geometry, state, control))
    output = open(sys.argv[2], "wb")
    output.write("case,node,rf_x_n,rf_y_n,rf_z_n\n")
    for copy_index, name in enumerate(names):
        for local in range(8):
            node = 8 * copy_index + local + 1
            value = reaction[node]
            output.write("%s,%d,%.16g,%.16g,%.16g\n" % (name, local + 1, value[0], value[1], value[2]))
    output.close()
finally:
    odb.close()
