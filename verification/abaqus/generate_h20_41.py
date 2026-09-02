#!/usr/bin/env python3
"""Generate the H20.41 finite-sliding partial-contact Abaqus input."""

from pathlib import Path


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError("expected exactly one H20.41 source fragment: %s" % old)
    return text.replace(old, new)


def main():
    directory = Path(__file__).resolve().parent
    source = directory / "h20_29_hex20_sts_graded.inp"
    target = directory / "h20_41_hex20_finite_sliding_partial_contact.inp"
    text = source.read_text(encoding="ascii")
    text = replace_once(
        text,
        "** H20.29: graded C3D20 small-sliding surface-to-surface contact.",
        "** H20.41: C3D20 finite-sliding surface-to-surface partial contact.",
    )
    text = replace_once(
        text,
        "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.",
        "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, adjust=0.",
    )
    text = replace_once(
        text,
        "SECONDARY_TOP, TRVEC, -1000, 1., 0., 0.",
        "SECONDARY_TOP, TRVEC, 20000, 1., 0., 0.",
    )
    target.write_text(text, encoding="ascii")
    print("wrote %s" % target)


if __name__ == "__main__":
    main()
