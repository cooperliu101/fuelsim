from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_b59.py <job.odb> <nodal.csv> <integration.csv>")


CASES = (
    ("coarse", 0, 18, 0, 4),
    ("refined", 1000, 75, 1000, 32),
    ("distorted", 2000, 75, 2000, 32),
)
STEP_TIME = 1.0e7


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


def case_node(label):
    for name, node_base, node_count, _, _ in CASES:
        if node_base < label <= node_base + node_count:
            return name, label - node_base
    raise RuntimeError("unexpected B5.9 node label %d" % label)


def case_element(label):
    for name, _, _, element_base, element_count in CASES:
        if element_base < label <= element_base + element_count:
            return name, label - element_base
    raise RuntimeError("unexpected B5.9 element label %d" % label)


odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    nodal = open(sys.argv[2], "wb")
    nodal.write(
        "case,stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n\n"
    )
    integration = open(sys.argv[3], "wb")
    integration.write(
        "case,stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,ivol_m3\n"
    )
    for stage in range(1, 5):
        frame = odb.steps["STEP%d" % stage].frames[-1]
        temperature = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "NT").values}
        displacement = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "U").values}
        reaction_flux = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RFL").values}
        reaction_force = {value.nodeLabel: double_data(value) for value in field_with_prefix(frame, "RF").values}
        for label in sorted(temperature.keys()):
            case, local = case_node(label)
            nodal.write(
                "%s,%d,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (case, stage, stage * STEP_TIME, local, temperature[label])
                    + tuple(displacement[label])
                    + (reaction_flux[label],)
                    + tuple(reaction_force[label])
                )
            )
        fields = {}
        for name in ("COORD", "E", "HFL", "IVOL", "S", "TEMP"):
            fields[name] = {
                (value.elementLabel, value.integrationPoint): double_data(value)
                for value in field_with_prefix(frame, name).values
                if value.elementLabel and value.integrationPoint
            }
        for key in sorted(fields["COORD"].keys()):
            label, point = key
            case, local = case_element(label)
            integration.write(
                "%s,%d,%.16g,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,"
                "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % (
                    (case, stage, stage * STEP_TIME, local, point)
                    + tuple(fields["COORD"][key])
                    + (fields["TEMP"][key],)
                    + tuple(fields["HFL"][key])
                    + tuple(fields["S"][key])
                    + tuple(fields["E"][key])
                    + (fields["IVOL"][key],)
                )
            )
    nodal.close()
    integration.close()
finally:
    odb.close()
