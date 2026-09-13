from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b54.py <job.odb> <capacity.csv> <loads.csv>")


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


def nodal_field(frame, name):
    result = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, name).values}
    if len(result) != 88:
        raise RuntimeError("%s has %d nodal values, expected 88" % (name, len(result)))
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["CAPACITY"].frames[-1]
    temperatures = nodal_field(frame, "NT")
    reaction_fluxes = nodal_field(frame, "RFL")
    capacity = open(sys.argv[2], "w", newline="")
    capacity.write("input_local_node,output_local_node,node,temperature_k,reaction_heat_flux_w\n")
    for case in range(8):
        for local_node in range(8):
            label = 8 * case + local_node + 1
            capacity.write(
                "%d,%d,%d,%.16g,%.16g\n"
                % (case + 1, local_node + 1, label, temperatures[label], reaction_fluxes[label])
            )
    capacity.close()

    load_output = open(sys.argv[3], "w", newline="")
    load_output.write("step,element,local_node,reaction_heat_flux_w\n")
    for step_name, element in (("BODY", 9), ("SURFACE", 10), ("FILM_BASE", 11), ("FILM_ACTIVE", 11)):
        reaction_flux = nodal_field(odb.steps[step_name].frames[-1], "RFL")
        for local_node in range(8):
            label = 8 * (element - 1) + local_node + 1
            load_output.write("%s,%d,%d,%.16g\n" % (step_name, element, local_node + 1, reaction_flux[label]))
    load_output.close()
finally:
    odb.close()
