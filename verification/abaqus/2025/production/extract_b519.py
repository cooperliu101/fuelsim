from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b519.py <job.odb> <nodal.csv> <integration.csv>")


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


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["PROBE"].frames[-1]
    coordinates = nodal_field(frame, "COORD")
    temperature = nodal_field(frame, "NT11")
    displacement = nodal_field(frame, "U")
    reaction_flux = nodal_field(frame, "RFL")
    reaction_force = nodal_field(frame, "RF")
    nodal = open(sys.argv[2], "w", newline="")
    nodal.write("node,x_m,y_m,z_m,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n\n")
    for node in range(1, 9):
        x = coordinates[node]
        u = displacement[node]
        rf = reaction_force[node]
        nodal.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                node,
                x[0],
                x[1],
                x[2],
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
    fields = {}
    for name in ("COORD", "IVOL", "LE", "S", "TEMP"):
        if name not in frame.fieldOutputs:
            raise RuntimeError("Abaqus B5.19 result does not contain exact %s output" % name)
        fields[name] = {
            (value.elementLabel, value.integrationPoint): double_data(value)
            for value in frame.fieldOutputs[name].values
            if value.elementLabel and value.integrationPoint
        }
    output = open(sys.argv[3], "w", newline="")
    output.write(
        "integration_point,x_m,y_m,z_m,temperature_k,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
        "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,ivol_m3\n"
    )
    for point in range(1, 9):
        key = (1, point)
        output.write(
            "%d,%.16g,%.16g,%.16g,%.16g,"
            "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
            "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                (point,)
                + tuple(fields["COORD"][key])
                + (fields["TEMP"][key],)
                + tuple(fields["S"][key])
                + tuple(fields["LE"][key])
                + (fields["IVOL"][key],)
            )
        )
    output.close()
finally:
    odb.close()
