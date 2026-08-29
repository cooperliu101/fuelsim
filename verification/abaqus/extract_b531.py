from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b531.py <job.odb> <nodal.csv> <integration.csv>")


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
    frame = odb.steps["TRANSIENT"].frames[-1]
    temperature = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "NT").values}
    displacement = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "U").values}
    reaction_flux = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RFL").values}
    reaction_force = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RF").values}
    nodal = open(sys.argv[2], "wb")
    nodal.write("node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n\n")
    for node in range(1, 13):
        nodal.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                (node, temperature[node])
                + tuple(displacement[node])
                + (reaction_flux[node],)
                + tuple(reaction_force[node])
            )
        )
    nodal.close()

    fields = {}
    for name in ("COORD", "E", "HFL", "IVOL", "S", "TEMP"):
        fields[name] = {
            (value.elementLabel, value.integrationPoint): double_data(value)
            for value in field_with_prefix(frame, name).values
            if value.elementLabel and value.integrationPoint
        }
    integration = open(sys.argv[3], "wb")
    integration.write(
        "element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,ivol_m3\n"
    )
    for element in range(1, 3):
        for point in range(1, 2):
            key = (element, point)
            integration.write(
                "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    key
                    + tuple(fields["COORD"][key])
                    + (fields["TEMP"][key],)
                    + tuple(fields["HFL"][key])
                    + tuple(fields["S"][key])
                    + tuple(fields["E"][key])
                    + (fields["IVOL"][key],)
                )
            )
    integration.close()
finally:
    odb.close()
