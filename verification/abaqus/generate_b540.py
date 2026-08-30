#!/usr/bin/env python3
"""Generate the deformable nonmatching C3D8RT thermal-mechanical contact path."""

from pathlib import Path

import generate_b524_b525


CASE = "b540_nonmatching_contact_cycle"


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    parameters = generate_b524_b525.configured(
        through=1,
        tangential=1,
        primary_tangential=2,
        secondary_tangential=1,
        y_bounds=((0.0, 2.0), (0.1, 0.9)),
        z_bounds=((-1.0, 2.0), (0.1, 0.9)),
        step=0.02,
        initial_gap=5.0e-4,
        elastic_only=True,
        path="nonmatching_contact_cycle",
        end_time=0.9,
        element_type="C3D8RT",
    )
    output = directory / (CASE + ".inp")
    output.write_text(generate_b524_b525.deck(CASE, parameters), encoding="utf-8")
    manifest = directory / "b540_case.tsv"
    manifest.write_text("case\texpected_frames\n%s\t45\n" % CASE, encoding="utf-8")
    print("wrote %s" % output)
    print("wrote %s" % manifest)
