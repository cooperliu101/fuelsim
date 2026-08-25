from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b52.py <job.odb> <nodal.csv>")


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def nodal_field(frame, prefix):
    field = frame.fieldOutputs[prefix] if prefix in frame.fieldOutputs else None
    if field is None:
        matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
        if len(matches) != 1:
            raise RuntimeError("expected one field beginning with %s, got %s" % (prefix, matches))
        field = frame.fieldOutputs[matches[0]]
    result = {}
    for value in field.values:
        if value.nodeLabel in range(1, 17):
            result[value.nodeLabel] = double_data(value)
    if len(result) != 16:
        raise RuntimeError("%s has %d nodal values, expected 16" % (prefix, len(result)))
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["HEAT"].frames[-1]
    temperatures = nodal_field(frame, "NT11")
    reaction_fluxes = nodal_field(frame, "RFL")
    output = open(sys.argv[2], "wb")
    output.write("node,temperature_k,reaction_heat_flux_w\n")
    for node in range(1, 17):
        output.write("%d,%.16g,%.16g\n" % (node, temperatures[node], reaction_fluxes[node]))
    output.close()
finally:
    odb.close()
