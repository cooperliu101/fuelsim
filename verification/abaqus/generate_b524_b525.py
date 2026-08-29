#!/usr/bin/env python3
"""Generate every B5.24, B5.25, and dedicated B5.26 Abaqus C3D8T full-field case."""

import math
from pathlib import Path


BASE = {
    "through": 2,
    "tangential": 1,
    "step": 0.02,
    "distortion": 0.0,
    "penalty": 1.0e9,
    "friction": 0.05,
    "slip_tolerance": 0.005,
    "conductance": 50.0,
    "pressure_conductance": 0.001,
    "primary_poisson": 0.28,
    "secondary_poisson": 0.30,
    "bending_traction": 0.0,
    "traction_controlled": False,
    "elastic_only": False,
    "anchor_bending": False,
    "initial_gap": 0.0,
    "path": "monotonic",
    "end_time": 0.2,
}


def configured(**changes):
    result = dict(BASE)
    result.update(changes)
    return result


CASES = {
    "b524_mesh_coarse": configured(),
    "b524_mesh_medium": configured(through=3),
    "b524_mesh_fine": configured(through=4),
    "b524_time_coarse": configured(step=0.04),
    "b524_time_fine": configured(step=0.01),
    "b524_penalty_low": configured(penalty=5.0e8),
    "b524_penalty_high": configured(penalty=2.0e9),
    "b524_friction_low": configured(friction=0.01),
    "b524_friction_high": configured(friction=0.10),
    "b524_slip_low": configured(slip_tolerance=0.0025),
    "b524_slip_high": configured(slip_tolerance=0.0100),
    "b524_thermal_low": configured(pressure_conductance=0.0005),
    "b524_thermal_high": configured(pressure_conductance=0.0020),
    "b525_poisson_030": configured(
        through=1,
        tangential=16,
        distortion=0.08,
        friction=0.20,
        primary_poisson=0.30,
        secondary_poisson=0.30,
        bending_traction=2.0e4,
        traction_controlled=True,
        elastic_only=True,
        anchor_bending=True,
    ),
    "b525_poisson_045": configured(
        through=1,
        tangential=16,
        distortion=0.08,
        friction=0.20,
        primary_poisson=0.45,
        secondary_poisson=0.45,
        bending_traction=2.0e4,
        traction_controlled=True,
        elastic_only=True,
        anchor_bending=True,
    ),
    "b525_poisson_049": configured(
        through=1,
        tangential=16,
        distortion=0.08,
        friction=0.20,
        primary_poisson=0.49,
        secondary_poisson=0.49,
        bending_traction=2.0e4,
        traction_controlled=True,
        elastic_only=True,
        anchor_bending=True,
    ),
    "b525_poisson_0499": configured(
        through=1,
        tangential=16,
        distortion=0.08,
        friction=0.20,
        primary_poisson=0.499,
        secondary_poisson=0.499,
        bending_traction=2.0e4,
        traction_controlled=True,
        elastic_only=True,
        anchor_bending=True,
    ),
    "b525_poisson_0499_refined": configured(
        through=2,
        tangential=16,
        distortion=0.08,
        friction=0.20,
        primary_poisson=0.499,
        secondary_poisson=0.499,
        bending_traction=2.0e4,
        traction_controlled=True,
        elastic_only=True,
        anchor_bending=True,
    ),
    "b526_contact_cycle": configured(elastic_only=True, initial_gap=1.0e-4, path="contact_cycle", end_time=0.3),
    "b526_friction_reversal": configured(elastic_only=True, path="friction_reversal", end_time=0.4),
}


def append_set(lines, keyword, name, labels):
    lines.append("*%s, %s=%s" % (keyword, keyword.lower(), name))
    for begin in range(0, len(labels), 16):
        lines.append(", ".join(str(value) for value in labels[begin : begin + 16]))


def build_mesh(parameters):
    nx = parameters["through"]
    ny = parameters["tangential"]
    distortion = parameters["distortion"]
    nodes = []
    elements = [[], []]
    node_maps = []
    for block in range(2):
        node_map = {}
        for z in range(2):
            for y in range(ny + 1):
                for x in range(nx + 1):
                    label = len(nodes) + 1
                    node_map[(z, y, x)] = label
                    nodes.append(
                        (
                            block
                            + (parameters["initial_gap"] if block == 1 else 0.0)
                            + float(x) / nx
                            + distortion * math.sin(math.pi * float(y) / ny),
                            float(y) / ny,
                            float(z),
                        )
                    )
        node_maps.append(node_map)
        for y in range(ny):
            for x in range(nx):
                elements[block].append(
                    (
                        node_map[(0, y, x)],
                        node_map[(0, y, x + 1)],
                        node_map[(0, y + 1, x + 1)],
                        node_map[(0, y + 1, x)],
                        node_map[(1, y, x)],
                        node_map[(1, y, x + 1)],
                        node_map[(1, y + 1, x + 1)],
                        node_map[(1, y + 1, x)],
                    )
                )
    return nodes, elements, node_maps


def material_lines(parameters):
    result = [
        "*Material, name=PRIMARY_MATERIAL",
        "*Elastic",
        "1.2000000000000000e8, %.16e, 3.0000000000000000e2" % parameters["primary_poisson"],
        "1.0000000000000000e8, %.16e, 5.0000000000000000e2" % parameters["primary_poisson"],
        "*Expansion, zero=3.0000000000000000e2",
        "8.0000000000000007e-6",
        "*Conductivity",
        "1.5000000000000000e1, 3.0000000000000000e2",
        "1.8000000000000000e1, 5.0000000000000000e2",
        "*Density",
        "1.0000000000000000e2",
        "*Specific Heat",
        "1.0000000000000000e0, 3.0000000000000000e2",
        "1.2000000000000000e0, 5.0000000000000000e2",
        "*Material, name=SECONDARY_MATERIAL",
        "*Elastic",
        "1.0000000000000000e8, %.16e, 3.0000000000000000e2" % parameters["secondary_poisson"],
        "9.0000000000000000e7, %.16e, 5.0000000000000000e2" % parameters["secondary_poisson"],
        "*Expansion, zero=3.0000000000000000e2",
        "1.0000000000000001e-5",
    ]
    if not parameters["elastic_only"]:
        result.extend(
            [
                "*Plastic",
                "2.0000000000000000e5, 0.0, 3.0000000000000000e2",
                "1.0200000000000000e7, 1.0, 3.0000000000000000e2",
                "1.8000000000000000e5, 0.0, 5.0000000000000000e2",
                "9.1800000000000000e6, 1.0, 5.0000000000000000e2",
                "*Creep, law=TIME",
                "1.2500000000000001e-20, 3.0, 0.0, 3.0000000000000000e2",
                "1.7500000000000000e-20, 3.0, 0.0, 5.0000000000000000e2",
            ]
        )
    result.extend(
        [
            "*Conductivity",
            "1.0000000000000000e1, 3.0000000000000000e2",
            "1.2000000000000000e1, 5.0000000000000000e2",
            "*Density",
            "1.0000000000000000e2",
            "*Specific Heat",
            "1.0000000000000000e0, 3.0000000000000000e2",
            "1.2000000000000000e0, 5.0000000000000000e2",
        ]
    )
    return result


def append_amplitude(lines, name, times, values):
    lines.append("*Amplitude, name=%s, time=TOTAL TIME" % name)
    pairs = ["%.16e, %.16e" % pair for pair in zip(times, values)]
    for begin in range(0, len(pairs), 4):
        lines.append(", ".join(pairs[begin : begin + 4]))


def deck(name, parameters):
    nodes, elements, node_maps = build_mesh(parameters)
    nx = parameters["through"]
    ny = parameters["tangential"]
    primary_count = len(elements[0])
    primary_labels = list(range(1, primary_count + 1))
    secondary_labels = list(range(primary_count + 1, primary_count + len(elements[1]) + 1))
    primary_contact_elements = [primary_labels[y * nx + nx - 1] for y in range(ny)]
    secondary_contact_elements = [secondary_labels[y * nx] for y in range(ny)]
    secondary_outer_elements = [secondary_labels[y * nx + nx - 1] for y in range(ny)]
    secondary_lower_elements = [secondary_outer_elements[y] for y in range(ny) if 2 * y < ny]
    secondary_upper_elements = [secondary_outer_elements[y] for y in range(ny) if 2 * y >= ny]
    primary_outer_nodes = [node_maps[0][(z, y, 0)] for z in range(2) for y in range(ny + 1)]
    secondary_outer_nodes = [node_maps[1][(z, y, nx)] for z in range(2) for y in range(ny + 1)]
    secondary_contact_nodes = [node_maps[1][(z, y, 0)] for z in range(2) for y in range(ny + 1)]
    secondary_y0_nodes = [node_maps[1][(z, 0, x)] for z in range(2) for x in range(nx + 1)]

    lines = [
        "*Heading",
        "** %s: generated full-field Abaqus reference." % name,
        "** through=%d tangential=%d step=%.16e distortion=%.16e initial_gap=%.16e path=%s"
        % (nx, ny, parameters["step"], parameters["distortion"], parameters["initial_gap"], parameters["path"]),
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.append("*Element, type=C3D8T, elset=PRIMARY")
    for label, connectivity in zip(primary_labels, elements[0]):
        lines.append("%d, %s" % (label, ", ".join(str(value) for value in connectivity)))
    lines.append("*Element, type=C3D8T, elset=SECONDARY")
    for label, connectivity in zip(secondary_labels, elements[1]):
        lines.append("%d, %s" % (label, ", ".join(str(value) for value in connectivity)))
    append_set(lines, "Elset", "PRIMARY_CONTACT_ELEMENTS", primary_contact_elements)
    append_set(lines, "Elset", "SECONDARY_CONTACT_ELEMENTS", secondary_contact_elements)
    append_set(lines, "Elset", "SECONDARY_OUTER_ELEMENTS", secondary_outer_elements)
    if parameters["traction_controlled"]:
        append_set(lines, "Elset", "SECONDARY_OUTER_LOWER_ELEMENTS", secondary_lower_elements)
        append_set(lines, "Elset", "SECONDARY_OUTER_UPPER_ELEMENTS", secondary_upper_elements)
    append_set(lines, "Nset", "ALL_NODES", list(range(1, len(nodes) + 1)))
    append_set(lines, "Nset", "PRIMARY_OUTER_NODES", primary_outer_nodes)
    append_set(lines, "Nset", "SECONDARY_OUTER_NODES", secondary_outer_nodes)
    append_set(lines, "Nset", "SECONDARY_CONTACT_NODES", secondary_contact_nodes)
    if parameters["anchor_bending"]:
        append_set(lines, "Nset", "SECONDARY_Y0", secondary_y0_nodes)
    lines.extend(
        [
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "PRIMARY_CONTACT_ELEMENTS, S4",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "SECONDARY_CONTACT_ELEMENTS, S6",
            "*Surface, type=ELEMENT, name=SECONDARY_OUTER",
            "SECONDARY_OUTER_ELEMENTS, S4",
        ]
    )
    if parameters["traction_controlled"]:
        lines.extend(
            [
                "*Surface, type=ELEMENT, name=SECONDARY_OUTER_LOWER",
                "SECONDARY_OUTER_LOWER_ELEMENTS, S4",
                "*Surface, type=ELEMENT, name=SECONDARY_OUTER_UPPER",
                "SECONDARY_OUTER_UPPER_ELEMENTS, S4",
            ]
        )
    lines.extend(material_lines(parameters))
    lines.extend(
        [
            "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL",
            ",",
            "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL",
            ",",
            "*Surface Interaction, name=COUPLED_CONTACT",
            "*Surface Behavior, pressure-overclosure=LINEAR",
            "%.16e" % parameters["penalty"],
            "*Friction, slip tolerance=%.16e" % parameters["slip_tolerance"],
            "%.16e" % parameters["friction"],
            "*Gap Conductance, pressure",
            "%.16e, 0.0" % parameters["conductance"],
            "%.16e, 5.0e5"
            % (parameters["conductance"] + 5.0e5 * parameters["pressure_conductance"]),
            "*Contact Pair, interaction=COUPLED_CONTACT, type=SURFACE TO SURFACE, adjust=0.",
            "SECONDARY_CONTACT, PRIMARY_CONTACT",
        ]
    )
    times = [0.0, 0.1, 0.2, 0.3, 0.4]
    if parameters["path"] == "contact_cycle":
        temperature = [300.0, 301.0, 301.0, 301.0, 301.0]
        pressure = [0.0] * 5
        normal = [0.0, -2.0e-4, 0.0, -2.0e-4, -2.0e-4]
        tangential_y = [0.0] * 5
        tangential_z = [0.0] * 5
    elif parameters["path"] == "friction_reversal":
        temperature = [300.0, 400.0, 400.0, 400.0, 400.0]
        pressure = [0.0] * 5
        normal = [0.0, -1.0e-3, -1.0e-3, -1.0e-3, -1.0e-3]
        tangential_y = [0.0, 2.0e-5, 1.2e-2, -4.0e-3, -3.98e-3]
        tangential_z = [0.0] * 5
    else:
        temperature = [300.0, 350.0, 400.0, 450.0, 500.0]
        pressure = [0.0, 8.75e4, 1.75e5, 2.625e5, 3.5e5]
        normal = [0.0] * 5
        tangential_y = [0.0, 5.0e-3, 1.0e-2, 1.5e-2, 2.0e-2]
        tangential_z = [0.0, 2.5e-3, 5.0e-3, 7.5e-3, 1.0e-2]
    append_amplitude(lines, "TEMPERATURE_PATH", times, temperature)
    append_amplitude(lines, "PRESSURE_PATH", times, pressure)
    append_amplitude(lines, "NORMAL_PATH", times, normal)
    append_amplitude(lines, "TANGENTIAL_Y_PATH", times, tangential_y)
    append_amplitude(lines, "TANGENTIAL_Z_PATH", times, tangential_z)
    append_amplitude(lines, "BENDING_PATH", times, [0.0, 0.25, 0.5, 0.75, 1.0])
    lines.extend(
        [
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Diagnostics, nonhybrid=WARNING",
            "*Step, name=PATH, nlgeom=YES, inc=160",
            "*Coupled Temperature-Displacement",
            "%.16e, %.16e, %.16e, %.16e"
            % (parameters["step"], parameters["end_time"], parameters["step"], parameters["step"]),
            "*Controls, parameters=FIELD, field=DISPLACEMENT",
            "1.0000000000000000e-10, 1.0000000000000000e-10, , , , 1.0000000000000000e-10",
            "*Boundary",
            "PRIMARY_OUTER_NODES, 1, 3, 0.0",
            "PRIMARY_OUTER_NODES, 11, 11, 3.0000000000000000e2",
            "*Boundary, amplitude=TEMPERATURE_PATH",
            "SECONDARY_OUTER_NODES, 11, 11, 1.0",
        ]
    )
    if parameters["anchor_bending"]:
        lines.extend(["*Boundary", "SECONDARY_Y0, 2, 3, 0.0"])
    if parameters["path"] != "monotonic":
        lines.extend(
            [
                "*Boundary, amplitude=NORMAL_PATH",
                "SECONDARY_OUTER_NODES, 1, 1, 1.0",
                "*Boundary, amplitude=TANGENTIAL_Y_PATH",
                "SECONDARY_OUTER_NODES, 2, 2, 1.0",
                "*Boundary, amplitude=TANGENTIAL_Z_PATH",
                "SECONDARY_OUTER_NODES, 3, 3, 1.0",
            ]
        )
    elif parameters["traction_controlled"]:
        lines.extend(
            [
                "*Dsload, follower=NO, amplitude=BENDING_PATH",
                "SECONDARY_OUTER_LOWER, TRVEC, %.16e, 0.0, 0.0, 1.0" % parameters["bending_traction"],
                "SECONDARY_OUTER_UPPER, TRVEC, %.16e, 0.0, 0.0, 1.0" % (-parameters["bending_traction"]),
            ]
        )
    else:
        lines.extend(
            [
                "*Boundary, amplitude=TANGENTIAL_Y_PATH",
                "SECONDARY_OUTER_NODES, 2, 2, 1.0",
                "*Boundary, amplitude=TANGENTIAL_Z_PATH",
                "SECONDARY_OUTER_NODES, 3, 3, 1.0",
            ]
        )
    if parameters["path"] == "monotonic":
        lines.extend(["*Dsload, follower=YES, amplitude=PRESSURE_PATH", "SECONDARY_OUTER, P, 1.0"])
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "CE, CEEQ, COORD, EE, HFL, IVOL, LE, PE, PEEQ, S, TEMP",
            "*Contact Output",
            "CDISP, CFORCE, CSTRESS, CTANDIR, CSTATUS, HFL",
            "*Output, history, frequency=1",
            "*Energy Output",
            "ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    manifest_header = (
        "case\tfamily\tthrough_thickness_elements\ttangential_elements\ttime_step_s\tdistortion_m\t"
        "penalty_pa_m\tfriction_coefficient\tslip_tolerance\tpressure_conductance_w_m2_k_pa\t"
        "primary_poisson\tsecondary_poisson\ttraction_controlled\telastic_only\tinitial_gap_m\tpath\tend_time_s\t"
        "expected_frames"
    )
    manifest_lines = [manifest_header]
    b526_manifest_lines = [manifest_header]
    for case_name in sorted(CASES):
        parameters = CASES[case_name]
        path = directory / (case_name + ".inp")
        path.write_text(deck(case_name, parameters), encoding="utf-8")
        print("wrote %s" % path)
        manifest_row = (
            "%s\t%s\t%d\t%d\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\t%d\t%d\t%.16g\t%s\t%.16g\t%d"
            % (
                case_name,
                case_name[:4],
                parameters["through"],
                parameters["tangential"],
                parameters["step"],
                parameters["distortion"],
                parameters["penalty"],
                parameters["friction"],
                parameters["slip_tolerance"],
                parameters["pressure_conductance"],
                parameters["primary_poisson"],
                parameters["secondary_poisson"],
                1 if parameters["traction_controlled"] else 0,
                1 if parameters["elastic_only"] else 0,
                parameters["initial_gap"],
                parameters["path"],
                parameters["end_time"],
                int(round(parameters["end_time"] / parameters["step"])),
            )
        )
        (b526_manifest_lines if case_name.startswith("b526_") else manifest_lines).append(manifest_row)
    manifest = directory / "b524_b525_cases.tsv"
    manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print("wrote %s" % manifest)
    b526_manifest = directory / "b526_cases.tsv"
    b526_manifest.write_text("\n".join(b526_manifest_lines) + "\n", encoding="utf-8")
    print("wrote %s" % b526_manifest)
