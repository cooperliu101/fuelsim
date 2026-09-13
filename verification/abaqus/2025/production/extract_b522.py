from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b522.py <job.odb> <nodal.csv> <contact.csv>")


def double_data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def nodal_field(frame, prefix, labels):
    field = frame.fieldOutputs[prefix] if prefix in frame.fieldOutputs else None
    if field is None:
        matches = [name for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
        if len(matches) != 1:
            raise RuntimeError("expected one field beginning with %s, got %s" % (prefix, matches))
        field = frame.fieldOutputs[matches[0]]
    result = {}
    for value in field.values:
        if value.nodeLabel in labels:
            result[value.nodeLabel] = double_data(value)
    if len(result) != len(labels):
        raise RuntimeError("%s has %d nodal values, expected %d" % (prefix, len(result), len(labels)))
    return result


def contact_field(frame, prefix):
    matches = [
        name
        for name in frame.fieldOutputs.keys()
        if name.strip().startswith(prefix) and "/" in name
    ]
    if len(matches) != 1:
        raise RuntimeError("expected one contact field beginning with %s, got %s" % (prefix, matches))
    return frame.fieldOutputs[matches[0]]


def contact_values(frame, prefix, labels):
    result = {}
    for value in contact_field(frame, prefix).values:
        if value.nodeLabel in labels:
            result[value.nodeLabel] = double_data(value)
    if len(result) != len(labels):
        raise RuntimeError("%s has %d contact values, expected %d" % (prefix, len(result), len(labels)))
    return result


def component(value, index):
    if hasattr(value, "__len__") and index < len(value):
        return value[index]
    return 0.0


def canonical_zero(value):
    return 0.0 if abs(value) < 1.0e-10 else value


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    instance = list(odb.rootAssembly.instances.values())[0]
    labels = tuple(sorted(node.label for node in instance.nodes))
    frame = odb.steps["LOAD"].frames[-1]
    temperature = nodal_field(frame, "NT11", labels)
    reaction_flux = nodal_field(frame, "RFL", labels)
    displacement = nodal_field(frame, "U", labels)
    reaction_force = nodal_field(frame, "RF", labels)
    output = open(sys.argv[2], "w", newline="")
    output.write("node,temperature_k,reaction_heat_flux_w,u1_m,u2_m,u3_m,rf1_n,rf2_n,rf3_n\n")
    for node in labels:
        output.write(
            "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                node,
                temperature[node],
                canonical_zero(reaction_flux[node]),
                displacement[node][0],
                displacement[node][1],
                displacement[node][2],
                canonical_zero(reaction_force[node][0]),
                canonical_zero(reaction_force[node][1]),
                canonical_zero(reaction_force[node][2]),
            )
        )
    output.close()

    opening_field = contact_field(frame, "COPEN")
    contact_labels = sorted(set(value.nodeLabel for value in opening_field.values))
    coordinates = nodal_field(frame, "COORD", labels)
    opening = contact_values(frame, "COPEN", contact_labels)
    pressure = contact_values(frame, "CPRESS", contact_labels)
    shear_stress1 = contact_values(frame, "CSHEAR1", contact_labels)
    shear_stress2 = contact_values(frame, "CSHEAR2", contact_labels)
    slip1 = contact_values(frame, "CSLIP1", contact_labels)
    slip2 = contact_values(frame, "CSLIP2", contact_labels)
    normal_force = contact_values(frame, "CNORMF", contact_labels)
    shear_force = contact_values(frame, "CSHEARF", contact_labels)
    contact_heat = contact_values(frame, "HFL", contact_labels)
    contact_output = open(sys.argv[3], "w", newline="")
    contact_output.write(
        "node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,"
        "shear_stress1_pa,shear_stress2_pa,"
        "normal_force1_n,normal_force2_n,normal_force3_n,"
        "shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w\n"
    )
    for node in contact_labels:
        contact_output.write(
            ("%d," + ",".join(["%.16g"] * 16) + "\n")
            % (
                (node,) + tuple(coordinates[node]) +
                (opening[node], pressure[node], slip1[node], slip2[node], shear_stress1[node], shear_stress2[node]) +
                tuple(component(normal_force[node], index) for index in range(3)) +
                tuple(component(shear_force[node], index) for index in range(3)) +
                (contact_heat[node],)
            )
        )
    contact_output.close()
finally:
    odb.close()
