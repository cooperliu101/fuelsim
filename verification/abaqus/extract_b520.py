from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b520.py <job.odb> <nodal.csv>")


STEPS = (
    "BASE",
    "GAP_PLUS",
    "GAP_MINUS",
    "SECONDARY_PLUS",
    "SECONDARY_MINUS",
    "PRIMARY_PLUS",
    "PRIMARY_MINUS",
)


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
    output = open(sys.argv[2], "wb")
    output.write("step,node,temperature_k,reaction_heat_flux_w,u1_m,reaction_force_x_n\n")
    for step_name in STEPS:
        frame = odb.steps[step_name].frames[-1]
        temperatures = nodal_field(frame, "NT11")
        reaction_fluxes = nodal_field(frame, "RFL")
        displacements = nodal_field(frame, "U")
        reaction_forces = nodal_field(frame, "RF")
        for node in range(1, 17):
            output.write(
                "%s,%d,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    step_name,
                    node,
                    temperatures[node],
                    reaction_fluxes[node],
                    displacements[node][0],
                    reaction_forces[node][0],
                )
            )
    output.close()
finally:
    odb.close()
