from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) not in (7, 8):
    raise RuntimeError(
        "usage: extract_b524_b525.py <job.odb> <nodal.csv> <integration.csv> "
        "<contact.csv> <energy.csv> <expected_frames> [reduced]"
    )


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def exact_field(frame, name):
    if name not in frame.fieldOutputs:
        raise RuntimeError("missing exact Abaqus field %s" % name)
    return frame.fieldOutputs[name]


def optional_field(frame, name):
    try:
        return frame.fieldOutputs[name]
    except Exception:
        return None


def contact_field(frame, prefix):
    matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix) and "/" in name]
    if len(matches) != 1:
        raise RuntimeError("expected one contact field beginning with %s, got %s" % (prefix, matches))
    return frame.fieldOutputs[matches[0]]


def nodal_values(field, labels):
    result = {}
    for value in field.values:
        if value.nodeLabel in labels:
            if value.nodeLabel in result:
                raise RuntimeError("duplicate nodal field value for node %d" % value.nodeLabel)
            result[value.nodeLabel] = double_data(value)
    if len(result) != len(labels):
        raise RuntimeError("nodal field has %d requested values, expected %d" % (len(result), len(labels)))
    return result


def integration_values(field):
    result = {}
    if field is None:
        return result
    for value in field.values:
        if value.elementLabel and value.integrationPoint:
            result[(value.elementLabel, value.integrationPoint)] = double_data(value)
    return result


def component(value, index):
    if hasattr(value, "__len__") and index < len(value):
        return value[index]
    return 0.0


def canonical_zero(value):
    return 0.0 if abs(value) < 1.0e-18 else value


def history_at(step, name, frame_index):
    matches = []
    for region in step.historyRegions.values():
        try:
            output = region.historyOutputs[name]
        except KeyError:
            continue
        matches.append(output.data)
    if not matches:
        return 0.0
    if len(matches) != 1:
        raise RuntimeError("expected at most one %s history, got %d" % (name, len(matches)))
    data = matches[0]
    if len(data) != len(step.frames):
        raise RuntimeError("%s history/frame count mismatch" % name)
    return data[frame_index][1]


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    instance = list(odb.rootAssembly.instances.values())[0]
    labels = sorted(node.label for node in instance.nodes)
    element_labels = sorted(element.label for element in instance.elements)
    contact_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
    reference_coordinates = dict((node.label, node.coordinates) for node in instance.nodes)
    expected_frames = int(sys.argv[6])

    nodal = open(sys.argv[2], "wb")
    integration = open(sys.argv[3], "wb")
    contact = open(sys.argv[4], "wb")
    energy = open(sys.argv[5], "wb")
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
        "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,"
        "normal_force1_n,normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,"
        "contact_heat_flux_w,shear_traction1_pa,shear_traction2_pa,"
        "tangent1_x,tangent1_y,tangent1_z,tangent2_x,tangent2_y,tangent2_z,state\n"
    )
    include_artificial_energy = len(sys.argv) == 8 and sys.argv[7] == "reduced"
    energy.write(
        "increment,time_s,allie_j,allse_j,allpd_j,allcd_j,allfd_j,allwk_j,boundary_heat_rate_w"
        + (",allae_j" if include_artificial_energy else "")
        + "\n"
    )

    step = odb.steps["PATH"]
    increment = 0
    zero6 = (0.0,) * 6
    for frame_index in range(1, len(step.frames)):
        increment += 1
        frame = step.frames[frame_index]
        time = frame.frameValue
        temperature = nodal_values(exact_field(frame, "NT11"), labels)
        displacement = nodal_values(exact_field(frame, "U"), labels)
        reaction_flux = nodal_values(exact_field(frame, "RFL11"), labels)
        reaction_force = nodal_values(exact_field(frame, "RF"), labels)
        for node in labels:
            nodal.write(
                "%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (increment, time, node, temperature[node])
                    + tuple(displacement[node])
                    + (canonical_zero(reaction_flux[node]),)
                    + tuple(canonical_zero(value) for value in reaction_force[node])
                )
            )

        fields = {}
        for name in ("COORD", "EE", "HFL", "IVOL", "LE", "S", "TEMP"):
            fields[name] = integration_values(exact_field(frame, name))
        for name in ("CE", "CEEQ", "PE", "PEEQ"):
            fields[name] = integration_values(optional_field(frame, name))
        for element in element_labels:
            for key in sorted(key for key in fields["S"] if key[0] == element):
                point = key[1]
                pe = fields["PE"].get(key, zero6)
                ce = fields["CE"].get(key, zero6)
                integration.write(
                    ("%d,%.16g,%d,%d," + ",".join(["%.16g"] * 40) + "\n")
                    % (
                        (increment, time, element, point)
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

        opening = nodal_values(contact_field(frame, "COPEN"), contact_labels)
        pressure = nodal_values(contact_field(frame, "CPRESS"), contact_labels)
        slip1 = nodal_values(contact_field(frame, "CSLIP1"), contact_labels)
        slip2 = nodal_values(contact_field(frame, "CSLIP2"), contact_labels)
        normal_force = nodal_values(contact_field(frame, "CNORMF"), contact_labels)
        shear_force = nodal_values(contact_field(frame, "CSHEARF"), contact_labels)
        shear_traction1 = nodal_values(contact_field(frame, "CSHEAR1"), contact_labels)
        shear_traction2 = nodal_values(contact_field(frame, "CSHEAR2"), contact_labels)
        contact_heat = nodal_values(contact_field(frame, "HFL"), contact_labels)
        tangent_first = nodal_values(exact_field(frame, "CTANDIR1"), contact_labels)
        tangent_second = nodal_values(exact_field(frame, "CTANDIR2"), contact_labels)
        status = nodal_values(contact_field(frame, "CSTATUS"), contact_labels)
        for node in contact_labels:
            current = tuple(
                reference_coordinates[node][index] + displacement[node][index] for index in range(3)
            )
            raw_status = status[node]
            if raw_status < -1.0e30:
                state = 0
            elif abs(raw_status - 1.0) < 1.0e-12:
                state = 2
            elif abs(raw_status - 2.0) < 1.0e-12:
                state = 1
            else:
                raise RuntimeError("unexpected Abaqus CSTATUS value %r at node %d" % (raw_status, node))
            contact.write(
                ("%d,%.16g,%d," + ",".join(["%.16g"] * 22) + ",%d\n")
                % (
                    (increment, time, node)
                    + current
                    + (opening[node], pressure[node], slip1[node], slip2[node])
                    + tuple(component(normal_force[node], index) for index in range(3))
                    + tuple(component(shear_force[node], index) for index in range(3))
                    + (contact_heat[node], shear_traction1[node], shear_traction2[node])
                    + tuple(component(tangent_first[node], index) for index in range(3))
                    + tuple(component(tangent_second[node], index) for index in range(3))
                    + (state,)
                )
            )

        energy_values = (
                increment,
                time,
                history_at(step, "ALLIE", frame_index),
                history_at(step, "ALLSE", frame_index),
                history_at(step, "ALLPD", frame_index),
                history_at(step, "ALLCD", frame_index),
                history_at(step, "ALLFD", frame_index),
                history_at(step, "ALLWK", frame_index),
                sum(reaction_flux.values()),
            )
        if include_artificial_energy:
            energy_values += (history_at(step, "ALLAE", frame_index),)
        energy.write(
            ("%d," + ",".join(["%.16g"] * (len(energy_values) - 1)) + "\n") % energy_values
        )
    if increment != expected_frames:
        raise RuntimeError("expected %d fixed increments, got %d" % (expected_frames, increment))
    nodal.close()
    integration.close()
    contact.close()
    energy.close()
    success = open(sys.argv[2] + ".ok", "wb")
    success.write("ok\n")
    success.close()
finally:
    odb.close()
