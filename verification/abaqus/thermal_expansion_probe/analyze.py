"""Read native outputs and identify expansion-temperature weights; no FE solve."""
import csv
import hashlib
import math
from pathlib import Path

root = Path(__file__).resolve().parent
coordinates = [(1.0, 0.0), (2.1, 0.1), (1.9, 1.2), (0.9, 1.0)]
signs = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
gauss = [(-math.sqrt(3 / 5), 5 / 9), (0.0, 8 / 9), (math.sqrt(3 / 5), 5 / 9)]
weights = [0.0] * 4
for xi, wx in gauss:
    for eta, wy in gauss:
        shape = [(1 + sx * xi) * (1 + sy * eta) / 4 for sx, sy in signs]
        dx = [sx * (1 + sy * eta) / 4 for sx, sy in signs]
        dy = [sy * (1 + sx * xi) / 4 for sx, sy in signs]
        radius = sum(n * x[0] for n, x in zip(shape, coordinates))
        a = sum(n * x[0] for n, x in zip(dx, coordinates))
        b = sum(n * x[0] for n, x in zip(dy, coordinates))
        c = sum(n * x[1] for n, x in zip(dx, coordinates))
        d = sum(n * x[1] for n, x in zip(dy, coordinates))
        measure = 2 * math.pi * radius * (a * d - b * c) * wx * wy
        for i in range(4):
            weights[i] += measure * shape[i]
volume = sum(weights)
fractions = [w / volume for w in weights]
print("volume_m3=" + str(volume))
print("volume_temperature_weights=" + str(fractions))
rows = []
for mode in ["small", "finite"]:
    case = "expansion_" + mode
    provenance = dict(line.split("=", 1) for line in (root / (case + "_run.txt")).read_text().splitlines())
    for name, key in [(root / (case + ".inp"), "input_sha256"),
                      (root.parent / "extract_b7_rz.py", "extractor_sha256")]:
        assert hashlib.sha256(name.read_bytes()).hexdigest().upper() == provenance[key]
    assert b"THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in (root / (case + ".sta")).read_bytes()
    with (root / (case + "_nodes.csv")).open() as f:
        for row in csv.DictReader(f):
            hot = round(float(row["time"]) / .1) - 1
            expected = 700 if (int(row["node"]) - 1) % 4 == hot else 600
            assert abs(float(row["temperature"]) - expected) < 1e-9
            assert abs(float(row["ur"])) < 1e-15 and abs(float(row["uz"])) < 1e-15
    with (root / (case + "_points.csv")).open() as f:
        for row in csv.DictReader(f):
            hot = round(float(row["time"]) / .1) - 1
            stresses = [float(row["stress_" + c]) for c in ["rr", "zz", "hoop"]]
            assert max(stresses) - min(stresses) < 1e-8
            assert abs(float(row["stress_rz"])) < 1e-8
            sigma = sum(stresses) / 3
            # sigma = -E*alpha*(T_effective-T_initial)/(1-2*nu).
            effective = 600 - sigma * (1 - 2 * .25) / (1e6 * 1e-5)
            elastic_effective = 600 - float(row["elastic_rr"]) / 1e-5
            assert abs(effective - elastic_effective) < 1e-9
            arithmetic = 625.0
            volumetric = 600 + 100 * fractions[hot]
            rows.append([mode, "CAX4T" if int(row["element"]) == 1 else "CAX4RT",
                         hot + 1, int(row["point"]), sigma, effective, arithmetic, volumetric,
                         effective - arithmetic, effective - volumetric])
with (root / "comparison.tsv").open("w") as f:
    out = csv.writer(f, delimiter="\t")
    out.writerow(["mode", "element", "hot_node", "point", "stress_pa", "inferred_temperature_k",
                  "arithmetic_temperature_k", "volume_temperature_k", "arithmetic_difference_k", "volume_difference_k"])
    out.writerows(rows)
for mode in ["small", "finite"]:
    for element in ["CAX4T", "CAX4RT"]:
        selected = [r for r in rows if r[0] == mode and r[1] == element]
        assert len(selected) == (16 if element == "CAX4T" else 4)
        print(mode, element, "inferred weights", [next((r[5] - 600) / 100 for r in selected if r[2] == n) for n in range(1, 5)])
        print("maximum temperature difference K: arithmetic", max(abs(r[8]) for r in selected),
              "volume", max(abs(r[9]) for r in selected))
