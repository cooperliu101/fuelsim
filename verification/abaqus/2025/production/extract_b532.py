from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b532.py <job.odb> <nodal.csv> <integration.csv>")


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
    return dict((value.nodeLabel, double_data(value)) for value in field_with_prefix(frame, prefix).values)


def integration_field(frame, prefix):
    values = [
        double_data(value)
        for value in field_with_prefix(frame, prefix).values
        if value.elementLabel == 1 and getattr(value, "integrationPoint", None) == 1
    ]
    if len(values) != 1:
        raise RuntimeError("expected one base integration value for %s, got %d" % (prefix, len(values)))
    return values[0]


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["FINITE_OPERATOR"].frames[-1]
    temperature = nodal_field(frame, "NT11")
    displacement = nodal_field(frame, "U")
    reaction_flux = nodal_field(frame, "RFL")
    reaction_force = nodal_field(frame, "RF")
    output = open(sys.argv[2], "w", newline="")
    output.write("case,node,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n\n")
    case_names = ["BASE"]
    for column in range(32):
        case_names.append("D%02d_PLUS" % column)
        case_names.append("D%02d_MINUS" % column)
    for case_index, case_name in enumerate(case_names):
        for local in range(8):
            node = 8 * case_index + local + 1
            u = displacement[node]
            rf = reaction_force[node]
            output.write(
                "%s,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (case_name, local + 1, temperature[node], u[0], u[1], u[2], reaction_flux[node], rf[0], rf[1], rf[2])
            )
    output.close()
    coordinates = integration_field(frame, "COORD")
    integration_temperature = integration_field(frame, "TEMP")
    heat_flux = integration_field(frame, "HFL")
    strain = integration_field(frame, "LE")
    stress = integration_field(frame, "S")
    volume = integration_field(frame, "IVOL")
    integration = open(sys.argv[3], "w", newline="")
    integration.write(
        "x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
        "e11,e22,e33,e12,e13,e23,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,ivol_m3\n"
    )
    integration.write(
        "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
        "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
        % (
            coordinates[0],
            coordinates[1],
            coordinates[2],
            integration_temperature,
            heat_flux[0],
            heat_flux[1],
            heat_flux[2],
            strain[0],
            strain[1],
            strain[2],
            strain[3],
            strain[4],
            strain[5],
            stress[0],
            stress[1],
            stress[2],
            stress[3],
            stress[4],
            stress[5],
            volume,
        )
    )
    integration.close()
finally:
    odb.close()
