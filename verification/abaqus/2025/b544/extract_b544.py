from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 6:
    raise RuntimeError("usage: extract_b544.py <job.odb> <nodal.csv> <integration.csv> <contact.csv> <energy.csv>")


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field(frame, name):
    if name in frame.fieldOutputs:
        return frame.fieldOutputs[name]
    matches = [candidate for candidate in frame.fieldOutputs.keys() if candidate.strip().startswith(name)]
    if len(matches) != 1:
        raise RuntimeError("expected one %s field, got %s" % (name, sorted(frame.fieldOutputs.keys())))
    return frame.fieldOutputs[matches[0]]


def nodal_values(frame, name):
    return dict((value.nodeLabel, data(value)) for value in field(frame, name).values if value.nodeLabel)


def integration_values(frame, name):
    return dict(
        ((value.elementLabel, value.integrationPoint), data(value))
        for value in field(frame, name).values
        if value.elementLabel and value.integrationPoint
    )


def component(value, index):
    if hasattr(value, "__len__") and index < len(value):
        return value[index]
    return 0.0


def canonical_zero(value):
    return 0.0 if abs(value) < 1.0e-18 else value


def history_at(step, name, frame_index):
    matches = []
    for region in step.historyRegions.values():
        if name in region.historyOutputs:
            matches.append(region.historyOutputs[name].data)
    if len(matches) != 1:
        raise RuntimeError("expected one %s history, got %d" % (name, len(matches)))
    return matches[0][frame_index][1]


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    nodal = open(sys.argv[2], "w", newline="")
    integration = open(sys.argv[3], "w", newline="")
    contact = open(sys.argv[4], "w", newline="")
    energy = open(sys.argv[5], "w", newline="")
    nodal.write("increment,time_s,node,temperature_k,u1_m,u2_m,u3_m,reaction_heat_flux_w,rf1_n,rf2_n,rf3_n\n")
    integration.write(
        "increment,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
        "hfl1_w_m2,hfl2_w_m2,hfl3_w_m2,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
        "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,"
        "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
        "pe11,pe22,pe33,pe12_engineering,pe13_engineering,pe23_engineering,peeq,"
        "ce11,ce22,ce33,ce12_engineering,ce13_engineering,ce23_engineering,ceeq,ivol_m3\n"
    )
    contact.write(
        "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,normal_force1_n,"
        "normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w,"
        "shear_traction1_pa,shear_traction2_pa,tangent1_x,tangent1_y,tangent1_z,tangent2_x,tangent2_y,"
        "tangent2_z,state\n"
    )
    energy.write(
        "increment,time_s,allie_j,allse_j,allpd_j,allcd_j,allfd_j,allwk_j,boundary_heat_rate_w,allae_j\n"
    )
    step = odb.steps["PATH"]
    if len(step.frames) != 11:
        raise RuntimeError("B5.44 expected ten fixed increments, got %d" % (len(step.frames) - 1))
    zero6 = (0.0,) * 6
    for increment in range(1, len(step.frames)):
        frame = step.frames[increment]
        temperature = nodal_values(frame, "NT11")
        displacement = nodal_values(frame, "U")
        reaction_flux = nodal_values(frame, "RFL11")
        reaction_force = nodal_values(frame, "RF")
        for node in range(1, 31):
            nodal.write(
                "%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (increment, frame.frameValue, node, temperature[node])
                    + tuple(displacement[node])
                    + (canonical_zero(reaction_flux[node]),)
                    + tuple(canonical_zero(value) for value in reaction_force[node])
                )
            )
        fields = {}
        for name in ("COORD", "EE", "HFL", "IVOL", "LE", "S", "TEMP"):
            fields[name] = integration_values(frame, name)
        for name in ("CE", "CEEQ", "PE", "PEEQ"):
            fields[name] = {}
        for element in range(1, 9):
            key = (element, 1)
            pe = fields["PE"].get(key, zero6)
            ce = fields["CE"].get(key, zero6)
            integration.write(
                ("%d,%.16g,%d,1," + ",".join(["%.16g"] * 40) + "\n")
                % (
                    (increment, frame.frameValue, element)
                    + tuple(fields["COORD"][key])
                    + (fields["TEMP"][key],)
                    + tuple(fields["HFL"][key])
                    + tuple(fields["S"][key])
                    + tuple(fields["LE"][key])
                    + tuple(fields["EE"][key])
                    + tuple(pe)
                    + (fields["PEEQ"].get(key, 0.0),)
                    + tuple(ce)
                    + (fields["CEEQ"].get(key, 0.0), fields["IVOL"][key])
                )
            )
        energy.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                increment,
                frame.frameValue,
                history_at(step, "ALLIE", increment),
                history_at(step, "ALLSE", increment),
                history_at(step, "ALLPD", increment),
                history_at(step, "ALLCD", increment),
                history_at(step, "ALLFD", increment),
                history_at(step, "ALLWK", increment),
                sum(reaction_flux.values()),
                history_at(step, "ALLAE", increment),
            )
        )
    nodal.close()
    integration.close()
    contact.close()
    energy.close()
    open(sys.argv[2] + ".ok", "w", newline="").close()
finally:
    odb.close()
