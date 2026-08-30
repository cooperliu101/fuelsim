#!/usr/bin/env python3
"""Generate C3D8RT thermal-mechanical contact transition and friction-reversal probes."""

from pathlib import Path

import generate_b524_b525


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    cases = {
        "b538_contact_cycle": dict(generate_b524_b525.CASES["b526_contact_cycle"]),
        "b539_friction_reversal": dict(generate_b524_b525.CASES["b526_friction_reversal"]),
    }
    manifest = ["case\texpected_frames"]
    for name, parameters in cases.items():
        parameters["element_type"] = "C3D8RT"
        if name == "b539_friction_reversal":
            parameters["normal_displacement_scale"] = 0.125
        output = directory / (name + ".inp")
        output.write_text(generate_b524_b525.deck(name, parameters), encoding="utf-8")
        expected_frames = int(round(parameters["end_time"] / parameters["step"]))
        manifest.append("%s\t%d" % (name, expected_frames))
        print("wrote %s" % output)
    output_manifest = directory / "b538_b539_cases.tsv"
    output_manifest.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print("wrote %s" % output_manifest)
