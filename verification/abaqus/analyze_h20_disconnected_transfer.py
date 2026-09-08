"""Compare assembled and separate-face native finite contact operators."""
import csv
from pathlib import Path
import numpy as np
from analyze_h20_28_finite_transfer import extract

root = Path(__file__).resolve().parent
connected = extract("h20_28_finite_transfer_shallow_probe")
separate = extract("h20_28_finite_disconnected_probe")
lookup = {tuple(point): index for index, point in enumerate(connected[3])}
merged = np.zeros_like(connected[2])
for column, point in enumerate(separate[3]):
    merged[:, lookup[tuple(point)]] += separate[2][:, column]
for label, actual, reference in [
        ("secondary_average", separate[0], connected[0]),
        ("constraint_area", separate[1], connected[1]),
        ("merged_primary_transfer", merged, connected[2])]:
    print(label + "_relative_difference=" + str(np.linalg.norm(actual-reference)/np.linalg.norm(reference)))
    print(label + "_maximum_absolute_difference=" + str(abs(actual-reference).max()))

with (root / "h20_28_finite_disconnected_probe_operator.csv").open() as stream:
    rows = [row for row in csv.DictReader(stream)
            if row["step"] == "BASE" and row["side"] == "primary"]
with (root / "h20_28_finite_disconnected_node_map.csv").open() as stream:
    metadata = {int(row["new_node"]): row for row in csv.DictReader(stream)}
faces = np.array([int(metadata[int(row["output_label"])]["primary_element"]) for row in rows])
assert len(rows) == 512 and all(np.count_nonzero(faces == face) == 8 for face in range(1, 65))
print("corner_constraint_separate_face_moments")
for face in range(1, 65):
    mask = faces == face
    weights, xy = separate[2][0, mask], separate[3][mask]
    if abs(weights).max() < 1e-7:
        continue
    mass = weights.sum()
    center = (weights @ xy) / mass if abs(mass) > 1e-14 else np.array([np.nan, np.nan])
    print("face", face, "weight", mass, "center", center.tolist(),
          "bounds", xy.min(axis=0).tolist(), xy.max(axis=0).tolist(),
          "minimum_coefficient", weights.min())
