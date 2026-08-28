from __future__ import print_function

import math
import sys

from odbAccess import openOdb


if len(sys.argv) != 6:
    raise RuntimeError(
        "usage: extract_b523.py <job.odb> <nodal.csv> <integration.csv> <contact.csv> <energy.csv>"
    )


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def exact_field(frame, name):
    if name not in frame.fieldOutputs:
        raise RuntimeError("missing exact Abaqus field %s; available fields are %s" % (name, sorted(frame.fieldOutputs.keys())))
    return frame.fieldOutputs[name]


def prefixed_field(frame, prefix):
    matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError("expected one field beginning with %s, got %s" % (prefix, matches))
    return frame.fieldOutputs[matches[0]]


def contact_field(frame, prefix):
    matches = [
        name
        for name in frame.fieldOutputs.keys()
        if name.strip().startswith(prefix) and "/" in name
    ]
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
    for value in field.values:
        if value.elementLabel and value.integrationPoint:
            result[(value.elementLabel, value.integrationPoint)] = double_data(value)
    return result


def contact_values(frame, prefix, labels):
    return nodal_values(contact_field(frame, prefix), labels)


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
    data = matches[0]
    if len(data) != len(step.frames):
        raise RuntimeError("%s history/frame count mismatch: %d versus %d" % (name, len(data), len(step.frames)))
    return data[frame_index][1]


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    instance = list(odb.rootAssembly.instances.values())[0]
    labels = list(range(1, 25))
    contact_labels = sorted(node.label for node in instance.nodeSets["SECONDARY_CONTACT_NODES"].nodes)
    reference_coordinates = dict((node.label, node.coordinates) for node in instance.nodes)
    if contact_labels != [13, 16, 19, 22]:
        raise RuntimeError("unexpected B5.23 secondary contact labels: %s" % contact_labels)

    nodal = open(sys.argv[2], "wb")
    integration = open(sys.argv[3], "wb")
    contact = open(sys.argv[4], "wb")
    energy = open(sys.argv[5], "wb")
    nodal.write(
        "increment,time_s,node,temperature_k,u1_m,u2_m,u3_m,reaction_heat_flux_w,rf1_n,rf2_n,rf3_n\n"
    )
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
        "contact_heat_flux_w,state\n"
    )
    energy.write(
        "increment,time_s,allie_j,allse_j,allpd_j,allcd_j,allfd_j,allwk_j,boundary_heat_rate_w\n"
    )

    increment = 0
    friction = 0.05
    step = odb.steps["PATH"]
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
        for name in ("CE", "CEEQ", "COORD", "EE", "HFL", "IVOL", "LE", "PE", "PEEQ", "S", "TEMP"):
            fields[name] = integration_values(exact_field(frame, name))
        zero6 = (0.0,) * 6
        for element in range(1, 5):
            for point in range(1, 9):
                key = (element, point)
                pe = fields["PE"].get(key, zero6)
                ce = fields["CE"].get(key, zero6)
                peeq = fields["PEEQ"].get(key, 0.0)
                ceeq = fields["CEEQ"].get(key, 0.0)
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
                        + (peeq,)
                        + tuple(ce)
                        + (ceeq, fields["IVOL"][key])
                    )
                )

        opening = contact_values(frame, "COPEN", contact_labels)
        pressure = contact_values(frame, "CPRESS", contact_labels)
        slip1 = contact_values(frame, "CSLIP1", contact_labels)
        slip2 = contact_values(frame, "CSLIP2", contact_labels)
        normal_force = contact_values(frame, "CNORMF", contact_labels)
        shear_force = contact_values(frame, "CSHEARF", contact_labels)
        contact_heat = contact_values(frame, "HFL", contact_labels)
        for node in contact_labels:
            current = tuple(
                reference_coordinates[node][component_index] + displacement[node][component_index]
                for component_index in range(3)
            )
            normal_magnitude = math.sqrt(sum(component(normal_force[node], i) ** 2 for i in range(3)))
            shear_magnitude = math.sqrt(sum(component(shear_force[node], i) ** 2 for i in range(3)))
            state = 0 if pressure[node] <= 0.0 else 2 if shear_magnitude >= 0.999999 * friction * normal_magnitude else 1
            contact.write(
                ("%d,%.16g,%d," + ",".join(["%.16g"] * 14) + ",%d\n")
                % (
                    (increment, time, node)
                    + current
                    + (opening[node], pressure[node], slip1[node], slip2[node])
                    + tuple(component(normal_force[node], i) for i in range(3))
                    + tuple(component(shear_force[node], i) for i in range(3))
                    + (contact_heat[node], state)
                )
            )

        energy.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
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
        )
    if increment != 20:
        raise RuntimeError("B5.23 expected twenty fixed increments, got %d" % increment)
    nodal.close()
    integration.close()
    contact.close()
    energy.close()
finally:
    odb.close()
