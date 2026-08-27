from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b57.py <job.odb> <capacity.csv>")


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field_with_prefix(frame, prefix):
    if prefix in frame.fieldOutputs:
        return frame.fieldOutputs[prefix]
    matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(
            "expected one field beginning with %s, got %s; available fields are %s"
            % (prefix, matches, sorted(frame.fieldOutputs.keys()))
        )
    return frame.fieldOutputs[matches[0]]


def nodal_field(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in range(1, 273):
            result[value.nodeLabel] = double_data(value)
    if len(result) != 272:
        raise RuntimeError("%s has %d nodal values, expected 272" % (prefix, len(result)))
    return result


def state_name(state):
    if state == 0:
        return "BASE"
    node = (state - 1) // 2
    sign = "PLUS" if (state - 1) % 2 == 0 else "MINUS"
    return "D%02d_%s" % (node, sign)


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["TRANSIENT"].frames[-1]
    temperatures = nodal_field(frame, "NT")
    reaction_fluxes = nodal_field(frame, "RFL")
    output = open(sys.argv[2], "wb")
    output.write("state,copy,local_node,node,temperature_k,reaction_heat_flux_w\n")
    for state in range(17):
        for transient_copy in (False, True):
            element = 2 * state + (1 if transient_copy else 0)
            copy_name = "transient" if transient_copy else "steady"
            for local_node in range(8):
                label = 8 * element + local_node + 1
                output.write(
                    "%s,%s,%d,%d,%.16g,%.16g\n"
                    % (
                        state_name(state),
                        copy_name,
                        local_node + 1,
                        label,
                        temperatures[label],
                        reaction_fluxes[label],
                    )
                )
    output.close()
finally:
    odb.close()
