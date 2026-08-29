from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError("usage: extract_b530.py <job.odb> <nodal.csv> <integration.csv> <surface_geometry.csv>")


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
    result = {}
    for value in field_with_prefix(frame, name).values:
        if value.nodeLabel in range(1, 65):
            result[value.nodeLabel] = double_data(value)
    if len(result) != 64:
        raise RuntimeError("%s has %d nodal values, expected 64" % (name, len(result)))
    return result


def integration_field(frame, name, element):
    result = {}
    for value in field_with_prefix(frame, name).values:
        if value.elementLabel == element and getattr(value, "integrationPoint", None) is not None:
            result[value.integrationPoint] = double_data(value)
    if len(result) != 1:
        raise RuntimeError("%s has %d values for element %d, expected 1" % (name, len(result), element))
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    output = open(sys.argv[2], "wb")
    output.write("step,element,local_node,temperature_k,reaction_heat_flux_w\n")
    for step_name in ("BODY", "SURFACE", "FILM_BASE", "FILM_ACTIVE"):
        frame = odb.steps[step_name].frames[-1]
        temperature = nodal_field(frame, "NT11")
        reaction_flux = nodal_field(frame, "RFL")
        for element in range(1, 4):
            for local_node in range(1, 9):
                node = 8 * (element - 1) + local_node
                output.write(
                    "%s,%d,%d,%.16g,%.16g\n"
                    % (step_name, element, local_node, temperature[node], reaction_flux[node])
                )
    output.close()
    points = open(sys.argv[3], "wb")
    points.write("step,element,integration_point,temperature_k,integration_volume_m3\n")
    for step_name in ("BODY", "SURFACE", "FILM_BASE", "FILM_ACTIVE"):
        frame = odb.steps[step_name].frames[-1]
        for element in range(1, 4):
            temperature = integration_field(frame, "TEMP", element)
            volume = integration_field(frame, "IVOL", element)
            for integration_point in range(1, 2):
                points.write(
                    "%s,%d,%d,%.16g,%.16g\n"
                    % (step_name, element, integration_point, temperature[integration_point], volume[integration_point])
                )
    points.close()
    surface = open(sys.argv[4], "wb")
    surface.write("element,local_node,reaction_heat_flux_w\n")
    reaction_flux = nodal_field(odb.steps["SURFACE_GEOMETRY"].frames[-1], "RFL")
    for element in (2, 4, 5, 6, 7, 8):
        for local_node in range(1, 9):
            node = 8 * (element - 1) + local_node
            surface.write("%d,%d,%.16g\n" % (element, local_node, reaction_flux[node]))
    surface.close()
finally:
    odb.close()
