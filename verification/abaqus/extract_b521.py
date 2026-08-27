from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b521.py <job.odb> <nodal.csv>")


STEPS = ("CLOSE", "SLIDE", "OPEN", "OPEN_CROSS", "RECONTACT")


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
        if value.nodeLabel in range(1, 21):
            result[value.nodeLabel] = double_data(value)
    if len(result) != 20:
        raise RuntimeError("%s has %d nodal values, expected 20" % (prefix, len(result)))
    return result


def canonical_zero(value):
    return 0.0 if abs(value) < 1.0e-10 else value


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    output = open(sys.argv[2], "wb")
    output.write("step,node,temperature_k,reaction_heat_flux_w,u1_m,u2_m,u3_m,rf1_n,rf2_n,rf3_n\n")
    for step_index, step_name in enumerate(STEPS, 1):
        frame = odb.steps[step_name].frames[-1]
        temperature = nodal_field(frame, "NT11")
        reaction_flux = nodal_field(frame, "RFL")
        displacement = nodal_field(frame, "U")
        reaction_force = nodal_field(frame, "RF")
        for node in range(1, 21):
            output.write(
                "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    step_index,
                    node,
                    temperature[node],
                    canonical_zero(reaction_flux[node]),
                    displacement[node][0],
                    displacement[node][1],
                    displacement[node][2],
                    canonical_zero(reaction_force[node][0]),
                    canonical_zero(reaction_force[node][1]),
                    canonical_zero(reaction_force[node][2]),
                )
            )
    output.close()
finally:
    odb.close()
