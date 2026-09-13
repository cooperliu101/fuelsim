from __future__ import print_function

import sys

from odbAccess import openOdb


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def nodal_field(frame, name):
    result = {}
    for value in frame.fieldOutputs[name].values:
        if value.nodeLabel in result:
            raise RuntimeError("Duplicate %s node %d" % (name, value.nodeLabel))
        result[value.nodeLabel] = double_data(value)
    if sorted(result) != list(range(1, 129)):
        raise RuntimeError("Incomplete nodal field %s" % name)
    return result


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_nonaffine_capacity.py job.odb nodal.csv integration.csv")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    step = odb.steps["PROBE"]
    if len(step.frames) != 2 or abs(step.frames[-1].frameValue - 1.0) > 1.0e-12:
        raise RuntimeError("Expected exactly one complete time increment")
    frame = step.frames[-1]
    temperature = nodal_field(frame, "NT11")
    reaction = nodal_field(frame, "RFL11")
    displacement = nodal_field(frame, "U")
    with open(sys.argv[2], "wb") as stream:
        stream.write("node,temperature_k,u1_m,u2_m,u3_m,reaction_heat_flux_w\n")
        for node in range(1, 129):
            values = [temperature[node]] + list(displacement[node]) + [reaction[node]]
            stream.write("%d,%s\n" % (node, ",".join("%.17g" % value for value in values)))
    ivol, material_temperature = {}, {}
    for name, target in [("IVOL", ivol), ("TEMP", material_temperature)]:
        for value in frame.fieldOutputs[name].values:
            if str(value.position) != "INTEGRATION_POINT":
                raise RuntimeError("Expected integration-point output for %s" % name)
            key = (value.elementLabel, value.integrationPoint)
            if key in target:
                raise RuntimeError("Duplicate integration point")
            target[key] = double_data(value)
        expected = [(element, point) for element in range(1, 17) for point in range(1, 9)]
        if sorted(target) != expected:
            raise RuntimeError("Incomplete integration-point field %s" % name)
    with open(sys.argv[3], "wb") as stream:
        stream.write("element,integration_point,temperature_k,ivol_m3\n")
        for key in sorted(ivol):
            stream.write("%d,%d,%.17g,%.17g\n" % (key[0], key[1], material_temperature[key], ivol[key]))
finally:
    odb.close()
