from pathlib import Path

root = Path(__file__).resolve().parents[1] / "production"
for field in ("nodes", "contact", "points"):
    source = root / ("c3d20rt_finite_contact_frictionless_" + field + ".csv")
    lines = source.read_bytes().splitlines(keepends=True)
    selected = [line for line in lines[1:] if float(line.split(b",", 1)[0]) <= 2.0000001]
    times = {float(line.split(b",", 1)[0]) for line in selected}
    if len(times) != 80:
        raise RuntimeError("Expected the existing 80 native output times")
    target = root / ("c3d20rt_finite_contact_frictionless_short_" + field + ".csv")
    target.write_bytes(b"".join([lines[0]] + selected))
