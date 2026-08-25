from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b51.py <job.odb> <nodal.csv> <integration_points.csv>")


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
    frame = odb.steps["HEAT"].frames[-1]
    coordinates = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "COORD").values if value.nodeLabel}
    temperatures = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "NT").values}
    reaction_fluxes = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RFL").values}
    displacements = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "U").values}
    nodal = open(sys.argv[2], "wb")
    nodal.write("node,x_current_m,y_current_m,z_current_m,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w\n")
    for node in range(1, 9):
        nodal.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % ((node,) + tuple(coordinates[node]) + (temperatures[node],) + tuple(displacements[node]) + (reaction_fluxes[node],))
        )
    nodal.close()

    ip_coordinates = [value for value in field_with_prefix(frame, "COORD").values if value.elementLabel]
    heat_fluxes = field_with_prefix(frame, "HFL").values
    if len(ip_coordinates) != 8 or len(heat_fluxes) != 8:
        raise RuntimeError("expected eight integration-point coordinates and heat fluxes")
    integration = open(sys.argv[3], "wb")
    integration.write("element,integration_point,x_current_m,y_current_m,z_current_m,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2\n")
    by_point = {(value.elementLabel, value.integrationPoint): double_data(value) for value in ip_coordinates}
    for value in sorted(heat_fluxes, key=lambda item: (item.elementLabel, item.integrationPoint)):
        key = (value.elementLabel, value.integrationPoint)
        integration.write(
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (key + tuple(by_point[key]) + tuple(double_data(value)))
        )
    integration.close()
finally:
    odb.close()
