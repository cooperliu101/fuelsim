from __future__ import print_function

import json
import sys

from odbAccess import openOdb


if len(sys.argv) not in (5, 6, 7):
    raise RuntimeError(
        "usage: extract_b61_c3d20t.py <job.odb> <mesh.json> <nodal.csv> <integration.csv> "
        "[frame_number [diagnostics]]"
    )
diagnostics = len(sys.argv) == 7
if diagnostics and sys.argv[6] != "diagnostics":
    raise RuntimeError("B6.1 C3D20T optional extraction mode must be diagnostics")


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


def integration_mises(frame):
    return dict(
        ((value.elementLabel, value.integrationPoint), value.mises)
        for value in field(frame, "S").values
        if value.elementLabel and value.integrationPoint
    )


with open(sys.argv[2], "rb") as source:
    mesh = json.load(source)
element_labels = sorted(element["label"] for block in mesh["blocks"] for element in block["elements"])
if element_labels != list(range(1, 1201)):
    raise RuntimeError("B6.1 C3D20T mesh manifest does not contain element labels 1 through 1200")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    step = odb.steps["FINITE_INELASTIC_BENDING"]
    if len(step.frames) != 6:
        raise RuntimeError("B6.1 C3D20T expected five fixed increments, got %d" % (len(step.frames) - 1))
    frame_number = len(step.frames) - 1 if len(sys.argv) == 5 else int(sys.argv[5])
    if frame_number < 1 or frame_number >= len(step.frames):
        raise RuntimeError("B6.1 C3D20T frame number must be between 1 and 5")
    frame = step.frames[frame_number]
    expected_time = 2.0 * frame_number
    if abs(frame.frameValue - expected_time) > 1.0e-10:
        raise RuntimeError(
            "B6.1 C3D20T frame %d has time %.16g instead of %.16g"
            % (frame_number, frame.frameValue, expected_time)
        )

    temperature = nodal_values(frame, "NT11")
    displacement = nodal_values(frame, "U")
    nodal_labels = sorted(displacement)
    with open(sys.argv[3], "wb") as output:
        output.write("id,x,y,z,temperature,displacement_x,displacement_y,displacement_z\n")
        for node in nodal_labels:
            coordinate = odb.rootAssembly.instances.values()[0].nodes[node - 1].coordinates
            output.write(
                "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n"
                % ((node,) + tuple(coordinate) + (temperature.get(node, 0.0),) + tuple(displacement[node]))
            )

    coordinates = integration_values(frame, "COORD")
    volumes = integration_values(frame, "IVOL")
    mises = integration_mises(frame)
    peeq = integration_values(frame, "PEEQ")
    ceeq = integration_values(frame, "CEEQ")
    expected = set((element, point) for element in element_labels for point in range(1, 28))
    for name, values in (("COORD", coordinates), ("IVOL", volumes), ("S", mises), ("PEEQ", peeq), ("CEEQ", ceeq)):
        if set(values) != expected:
            raise RuntimeError("B6.1 C3D20T %s does not contain all 27 integration points" % name)
    if diagnostics:
        heat_fluxes = integration_values(frame, "HFL")
        stresses = integration_values(frame, "S")
        elastic_strains = integration_values(frame, "EE")
        plastic_strains = integration_values(frame, "PE")
        creep_strains = integration_values(frame, "CE")
        for name, values in (("HFL", heat_fluxes), ("S", stresses), ("EE", elastic_strains),
                             ("PE", plastic_strains), ("CE", creep_strains)):
            if set(values) != expected:
                raise RuntimeError("B6.1 C3D20T %s does not contain all 27 integration points" % name)
    with open(sys.argv[4], "wb") as output:
        if diagnostics:
            output.write(
                "element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress,peeq,ceeq,"
                "hfl1_w_m2,hfl2_w_m2,hfl3_w_m2,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
                "ee11,ee22,ee33,ee12,ee13,ee23,"
                "pe11,pe22,pe33,pe12,pe13,pe23,"
                "ce11,ce22,ce33,ce12,ce13,ce23\n"
            )
        else:
            output.write(
                "element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress,peeq,ceeq\n"
            )
        for element in element_labels:
            for point in range(1, 28):
                key = (element, point)
                common = (element, point) + tuple(coordinates[key]) + (volumes[key], mises[key], peeq[key], ceeq[key])
                if diagnostics:
                    output.write(
                        ("%d,%d," + ",".join(["%.16g"] * 34) + "\n")
                        % (
                            common
                            + tuple(heat_fluxes[key])
                            + tuple(stresses[key])
                            + tuple(elastic_strains[key])
                            + tuple(plastic_strains[key])
                            + tuple(creep_strains[key])
                        )
                    )
                else:
                    output.write("%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n" % common)
finally:
    odb.close()
