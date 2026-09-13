from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 5:
    raise RuntimeError("usage: extract_b510.py <job.odb> <nodal.csv> <integration.csv> <energy.csv>")


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


def history_value(step, name):
    matches = []
    for region in step.historyRegions.values():
        if name in region.historyOutputs:
            matches.append(region.historyOutputs[name].data[-1][1])
    if len(matches) != 1:
        raise RuntimeError("expected one %s history, got %s" % (name, matches))
    return matches[0]


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    first_frame = odb.steps["STEP1"].frames[-1]
    strain_name = "E" if "E" in first_frame.fieldOutputs else "LE" if "LE" in first_frame.fieldOutputs else None
    if strain_name is None:
        raise RuntimeError("Abaqus result contains neither exact E nor exact LE strain output")
    strain_prefix = strain_name.lower()
    nodal = open(sys.argv[2], "w", newline="")
    nodal.write("stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n\n")
    integration = open(sys.argv[3], "w", newline="")
    integration.write(
        "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
        "%s11,%s22,%s33,%s12_engineering,%s13_engineering,%s23_engineering,"
        "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
        "pe11,pe22,pe33,pe12_engineering,pe13_engineering,pe23_engineering,peeq,ivol_m3\n"
        % ((strain_prefix,) * 6)
    )
    energy = open(sys.argv[4], "w", newline="")
    energy.write("stage,time_s,internal_energy_j,plastic_dissipation_j,external_work_j\n")
    for stage in range(1, 11):
        step = odb.steps["STEP%d" % stage]
        frame = step.frames[-1]
        time = 0.1 * stage
        temperature = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "NT").values}
        displacement = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "U").values}
        reaction_flux = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RFL").values}
        reaction_force = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RF").values}
        for node in range(1, 9):
            nodal.write(
                "%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (stage, time, node, temperature[node])
                    + tuple(displacement[node])
                    + (reaction_flux[node],)
                    + tuple(reaction_force[node])
                )
            )
        fields = {}
        for name in ("COORD", "EE", "IVOL", "PE", "PEEQ", "S", "TEMP"):
            fields[name] = {
                (value.elementLabel, value.integrationPoint): double_data(value)
                for value in field_with_prefix(frame, name).values
                if value.elementLabel and value.integrationPoint
            }
        fields["STRAIN"] = {
            (value.elementLabel, value.integrationPoint): double_data(value)
            for value in frame.fieldOutputs[strain_name].values
            if value.elementLabel and value.integrationPoint
        }
        for key in sorted(key for key in fields["S"] if key[0] == 1):
            point = key[1]
            integration.write(
                "%d,%.16g,1,%d,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (stage, time, point)
                    + tuple(fields["COORD"][key])
                    + (fields["TEMP"][key],)
                    + tuple(fields["S"][key])
                    + tuple(fields["STRAIN"][key])
                    + tuple(fields["EE"][key])
                    + tuple(fields["PE"][key])
                    + (fields["PEEQ"][key], fields["IVOL"][key])
                )
            )
        energy.write(
            "%d,%.16g,%.16g,%.16g,%.16g\n"
            % (stage, time, history_value(step, "ALLIE"), history_value(step, "ALLPD"), history_value(step, "ALLWK"))
        )
    nodal.close()
    integration.close()
    energy.close()
finally:
    odb.close()
