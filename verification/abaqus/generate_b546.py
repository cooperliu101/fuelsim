#!/usr/bin/env python3
"""Generate the full-size M5.8-equivalent Abaqus C3D8T or C3D8RT cases."""

from pathlib import Path

from exodus_to_abaqus import convert


JOB = "b546_m58_c3d8t_integrated"


def amplitude(lines, name, pairs):
    lines.append("*Amplitude, name=%s, time=TOTAL TIME" % name)
    values = []
    for time, value in pairs:
        values.extend(("%.16g" % time, "%.16g" % value))
    for begin in range(0, len(values), 8):
        lines.append(", ".join(values[begin : begin + 8]))


def analytic_fuel_conductivity(lines):
    lines.append("*Conductivity")
    temperatures = [300.0, 590.0]
    temperatures.extend(590.0 + 0.1 * index for index in range(1, 301))
    temperatures.append(1200.0)
    for temperature in temperatures:
        lines.append("%.16g, %.16g" % (3824.0 / temperature + 0.61, temperature))


def analytic_gap_conductance(lines):
    lines.append("*Gap Conductance")
    lines.append("4000, 0")
    lines.append("4000, 1e-6")
    for index in range(1, 901):
        gap = 1.0e-6 + index * 1.0e-8
        lines.append("%.16g, %.16g" % (0.004 / gap, gap))
    for gap in (2.0e-5, 5.0e-5, 1.0e-4, 5.0e-4, 1.0e-3):
        lines.append("%.16g, %.16g" % (0.004 / gap, gap))


def deck(job=JOB, case_name="B5.46", element_type="C3D8T"):
    lines = [
        "*Heading",
        "** %s full-size M5.8-equivalent finite-strain %s benchmark." % (case_name, element_type),
        "** The mesh include is converted directly from m58_integrated_hex8_mesh.e.",
        "** SI units: metre, second, kelvin, pascal, watt.",
        "** Temperature and gap tables sample the stated analytic functions; no coefficient is fitted.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Include, input=%s_mesh.inc" % job,
        "*Material, name=FUEL_MATERIAL",
        "*Elastic",
        "2e11, 0.316",
        "*Expansion, zero=600",
        "1e-5",
    ]
    analytic_fuel_conductivity(lines)
    lines.extend(
        [
            "*Density",
            "10970",
            "*Specific Heat",
            "300",
            "*Material, name=CLAD_MATERIAL",
            "*Elastic",
            "7.5e10, 0.3",
            "*Expansion, zero=600",
            "5e-6",
            "*Plastic",
            "1e6, 0",
            "2.0001e10, 1",
            "*Creep, law=TIME",
            "8e-28, 3, 0",
            "*Conductivity",
            "16",
            "*Density",
            "6500",
            "*Specific Heat",
            "330",
            "*Solid Section, elset=FUEL, material=FUEL_MATERIAL",
            ",",
            "*Solid Section, elset=CLAD, material=CLAD_MATERIAL",
            ",",
            "*Surface Interaction, name=FUEL_CLADDING_CONTACT",
            "*Surface Behavior, pressure-overclosure=LINEAR",
            "1e10",
            "*Friction",
            "0.002",
        ]
    )
    analytic_gap_conductance(lines)
    lines.extend(
        [
            "*Contact Pair, interaction=FUEL_CLADDING_CONTACT, type=NODE TO SURFACE, adjust=0.",
            "FUEL_OUTER, CLAD_RMIN",
        ]
    )
    amplitude(lines, "POWER_PATH", ((0, 0), (0.2, 0.5), (0.5, 1), (0.8, 1.2), (1, 1)))
    amplitude(lines, "INTERNAL_PRESSURE_PATH", ((0, 0), (0.2, 0.3), (0.5, 1), (0.8, 0.8), (1, 1)))
    amplitude(lines, "EXTERNAL_PRESSURE_PATH", ((0, 0), (0.2, 0.8), (0.5, 1), (0.8, 1.2), (1, 1)))
    amplitude(lines, "FUEL_AXIAL_PATH", ((0, 0), (0.2, 6e-5), (0.5, 3e-4), (0.8, 4.8e-4), (1, 6e-4)))
    amplitude(lines, "CLAD_AXIAL_PATH", ((0, 0), (0.2, 2e-6), (0.5, 1e-5), (0.8, 1.6e-5), (1, 2e-5)))
    lines.extend(
        [
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 600",
            "*Step, name=PATH, nlgeom=YES, inc=20",
            "*Coupled Temperature-Displacement",
            "0.05, 1, 0.05, 0.05",
            "*Controls, parameters=FIELD, field=DISPLACEMENT",
            "1e-10, 1e-10, , , , 1e-10",
            "*Boundary",
            "FUEL_BOTTOM, 1, 2, 0",
            "CLAD_BOTTOM, 1, 3, 0",
            "CLAD_RMAX, 11, 11, 600",
            "*Boundary, amplitude=FUEL_AXIAL_PATH",
            "FUEL_TOP, 3, 3, 1",
            "*Boundary, amplitude=CLAD_AXIAL_PATH",
            "CLAD_TOP, 3, 3, 1",
            "*Dflux, amplitude=POWER_PATH",
            "FUEL, BF, 2e7",
            "*Dsload, follower=YES, amplitude=INTERNAL_PRESSURE_PATH",
            "CLAD_RMIN, P, 5e5",
            "*Dsload, follower=YES, amplitude=EXTERNAL_PRESSURE_PATH",
            "CLAD_RMAX, P, 2e6",
            "*Output, field, frequency=20",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "CEEQ, IVOL, PEEQ, S",
            "*Contact Output",
            "CDISP, CFORCE, CSTRESS, HFL",
            "*Output, history, frequency=%d" % (1 if element_type == "C3D8RT" else 20),
            "*Energy Output",
            (
                "ALLAE, ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK"
                if element_type == "C3D8RT"
                else "ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK"
            ),
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


def generate(job=JOB, case_name="B5.46", element_type="C3D8T"):
    directory = Path(__file__).resolve().parent
    exodus = directory.parent / "moose" / "m58_integrated_hex8_mesh.e"
    mesh_include = directory / (job + "_mesh.inc")
    manifest = directory / (job + "_mesh.json")
    mesh = convert(exodus, mesh_include, manifest, element_type)
    input_path = directory / (job + ".inp")
    input_path.write_text(deck(job, case_name, element_type), encoding="ascii")
    if len(mesh["nodes"]) != 1617 or len(mesh["elements"]) != 1152:
        raise RuntimeError("%s expected 1617 nodes and 1152 elements" % case_name)
    expected_blocks = ["FUEL", "CLAD"]
    expected_sets = {"FUEL_OUTER", "FUEL_BOTTOM", "FUEL_TOP", "CLAD_RMIN", "CLAD_RMAX", "CLAD_BOTTOM", "CLAD_TOP"}
    if [block["name"] for block in mesh["blocks"]] != expected_blocks:
        raise RuntimeError("%s element block names changed" % case_name)
    if {side_set["name"] for side_set in mesh["side_sets"]} != expected_sets:
        raise RuntimeError("%s side-set names changed" % case_name)
    print("wrote %s, %s, and %s" % (input_path, mesh_include, manifest))


def main():
    generate()


if __name__ == "__main__":
    main()
