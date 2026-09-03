#!/usr/bin/env python3
"""Generate the 376-degree-of-freedom B5.55 Abaqus C3D20T friction case."""

from pathlib import Path

from exodus_to_abaqus import convert
from generate_b549 import deck


JOB = "b555_small_c3d20t_friction"
MESH = "b549_small_c3d20t_mesh"


def contact_deck(mesh):
    text = deck(mesh)
    text = text.replace(
        "B5.49 small C3D20T finite-strain thermo-inelastic volume-element isolation case.",
        "B5.55 small C3D20T finite-strain thermo-inelastic friction case.",
    )
    text = text.replace("contact=NO", "contact=YES")
    text = text.replace(
        "*Amplitude, name=SECONDARY_TEMPERATURE, time=TOTAL TIME",
        "\n".join(
            [
                "*Surface Interaction, name=COUPLED_CONTACT",
                "*Surface Behavior, pressure-overclosure=LINEAR",
                "1e9",
                "*Friction, slip tolerance=1e-4",
                "0.2",
                "*Gap Heat Generation",
                "0, 0.5",
                "*Gap Conductance, pressure",
                "50, 0",
                "550, 5e5",
                "*Contact Pair, interaction=COUPLED_CONTACT, type=SURFACE TO SURFACE, adjust=0.",
                "SECONDARY_CONTACT, PRIMARY_CONTACT",
                "*Amplitude, name=SECONDARY_TEMPERATURE, time=TOTAL TIME",
            ]
        ),
    )
    text = text.replace(
        "SECONDARY_CONTACT_NODES, 1, 1, 0\nSECONDARY_CONTACT_TEMPERATURE, 11, 11, 300\n",
        "",
    )
    text = text.replace(
        "PRIMARY_CONTACT_TEMPERATURE, 11, 11, 1\nSECONDARY_OUTER_TEMPERATURE, 11, 11, 1",
        "SECONDARY_OUTER_TEMPERATURE, 11, 11, 1",
    )
    text = text.replace(
        "*Output, history, frequency=20",
        "*Contact Output\nCDISP, CFORCE, CSTRESS, CTANDIR, HFL\n*Output, history, frequency=20",
    )
    text = text.replace("0, 0, 0.4, 0\n*Amplitude, name=TANGENTIAL_Z", "0, 0, 0.4, 0.0003\n*Amplitude, name=TANGENTIAL_Z")
    text = text.replace("0, 0, 0.4, 0\n*Initial Conditions", "0, 0, 0.4, 0.0004\n*Initial Conditions")
    return text


def main():
    directory = Path(__file__).resolve().parent
    mesh = convert(
        directory / (MESH + ".e"),
        directory / (MESH + ".inc"),
        directory / (MESH + ".json"),
        "C3D20T",
    )
    input_path = directory / (JOB + ".inp")
    input_path.write_text(contact_deck(mesh), encoding="ascii")
    print("wrote %s and refreshed the shared B5.49 mesh artifacts" % input_path)


if __name__ == "__main__":
    main()
