from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b528.py <job.odb> <nodal.csv> <integration_points.csv>")


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
        if value.nodeLabel in range(1, 9):
            result[value.nodeLabel] = double_data(value)
    if len(result) != 8:
        raise RuntimeError("%s has %d nodal values, expected 8" % (prefix, len(result)))
    return result


def integration_field(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.elementLabel == 1 and getattr(value, "integrationPoint", None) is not None:
            result[value.integrationPoint] = double_data(value)
    if len(result) != 1:
        raise RuntimeError("%s has %d integration-point values, expected 1" % (prefix, len(result)))
    return result


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    nodal = open(sys.argv[2], "wb")
    nodal.write("step,node,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n\n")
    step_names = ["BASE"]
    for column in range(32):
        step_names.append("D%02d_PLUS" % column)
        step_names.append("D%02d_MINUS" % column)
    for step_name in step_names:
        if step_name not in odb.steps:
            raise RuntimeError("missing step %s" % step_name)
        frame = odb.steps[step_name].frames[-1]
        temperature = nodal_field(frame, "NT11")
        displacement = nodal_field(frame, "U")
        reaction_flux = nodal_field(frame, "RFL")
        reaction_force = nodal_field(frame, "RF")
        for node in range(1, 9):
            u = displacement[node]
            rf = reaction_force[node]
            nodal.write(
                "%s,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    step_name,
                    node,
                    temperature[node],
                    u[0],
                    u[1],
                    u[2],
                    reaction_flux[node],
                    rf[0],
                    rf[1],
                    rf[2],
                )
            )
    nodal.close()

    frame = odb.steps["BASE"].frames[-1]
    coordinates = integration_field(frame, "COORD")
    temperature = integration_field(frame, "TEMP")
    heat_flux = integration_field(frame, "HFL")
    strain = integration_field(frame, "E")
    stress = integration_field(frame, "S")
    points = open(sys.argv[3], "wb")
    points.write(
        "element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
        "e11,e22,e33,e12,e13,e23,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa\n"
    )
    for integration_point in range(1, 2):
        x = coordinates[integration_point]
        hfl = heat_flux[integration_point]
        e = strain[integration_point]
        s = stress[integration_point]
        points.write(
            "1,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
            "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                integration_point,
                x[0],
                x[1],
                x[2],
                temperature[integration_point],
                hfl[0],
                hfl[1],
                hfl[2],
                e[0],
                e[1],
                e[2],
                e[3],
                e[4],
                e[5],
                s[0],
                s[1],
                s[2],
                s[3],
                s[4],
                s[5],
            )
        )
    points.close()
finally:
    odb.close()
