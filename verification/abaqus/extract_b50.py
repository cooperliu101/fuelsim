from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b50.py <job.odb> <capacity.csv>")


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


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["CAPACITY"].frames[-1]
    temperatures = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "NT").values}
    reaction_fluxes = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RFL").values}
    if len(temperatures) != 64 or len(reaction_fluxes) != 64:
        raise RuntimeError("expected 64 temperatures and reaction heat fluxes")
    output = open(sys.argv[2], "wb")
    output.write("input_local_node,output_local_node,node,temperature_k,reaction_heat_flux_w\n")
    for case in range(8):
        for local_node in range(8):
            label = 8 * case + local_node + 1
            output.write(
                "%d,%d,%d,%.16g,%.16g\n"
                % (case + 1, local_node + 1, label, temperatures[label], reaction_fluxes[label])
            )
    output.close()
finally:
    odb.close()
