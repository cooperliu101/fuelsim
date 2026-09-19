"""Compare archived production outputs; never rerun or modify a physical case."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def physical_csv(path):
    with path.open() as stream:
        rows = list(csv.reader(stream))
    if rows and rows[0] == ["metric", "value"]:
        paths = {"input_file", "mesh_file", "results_file", "history_file", "checkpoint_file", "restart_file"}
        rows = [row for row in rows if row and not (
            row[0].startswith(("aggregate_memory.", "aggregate_timing.", "memory."))
            or "seconds" in row[0] or row[0] in paths)]
    return rows


def compare(baseline, current):
    report = dict(exodus_files=0, numeric_arrays=0, numeric_values=0,
                  checkpoint_files=0, csv_files=0, differences=[], missing=[])
    aliases = {}
    # Restart files may have different segment numbers in a reused output directory.
    # Compare the files actually named by the corresponding production summaries.
    for summary in baseline.rglob("*summary.csv"):
        relative = summary.relative_to(baseline)
        updated = current / relative
        if not updated.exists():
            continue
        with summary.open() as stream:
            old = dict(list(csv.reader(stream))[1:])
        with updated.open() as stream:
            new = dict(list(csv.reader(stream))[1:])
        for field in ("history_file", "results_file"):
            if field in old and field in new:
                aliases[relative.parent / Path(old[field]).name] = current / relative.parent / Path(new[field]).name
    for old in sorted(baseline.rglob("*")):
        if not old.is_file():
            continue
        relative = old.relative_to(baseline)
        new = aliases.get(relative, current / relative)
        if not new.exists():
            report["missing"].append(str(relative))
            continue
        if old.suffix == ".e":
            report["exodus_files"] += 1
            with Dataset(old) as first, Dataset(new) as second:
                first.set_auto_maskandscale(False)
                second.set_auto_maskandscale(False)
                if set(first.variables) != set(second.variables):
                    report["differences"].append(dict(file=str(relative), reason="variable names differ"))
                    continue
                for name in first.variables:
                    if first[name].dtype.kind not in "fiu":
                        continue
                    a, b = np.asarray(first[name][:]), np.asarray(second[name][:])
                    report["numeric_arrays"] += 1
                    report["numeric_values"] += int(a.size)
                    if a.shape != b.shape or not np.array_equal(a, b, equal_nan=True):
                        difference = dict(file=str(relative), variable=name, shape_old=a.shape, shape_new=b.shape)
                        if a.shape == b.shape:
                            finite = np.isfinite(a) & np.isfinite(b)
                            difference["maximum_absolute_difference"] = float(np.max(np.abs(a[finite]-b[finite]))) if finite.any() else None
                        report["differences"].append(difference)
        elif old.suffix == ".checkpoint":
            report["checkpoint_files"] += 1
            if old.read_bytes() != new.read_bytes():
                report["differences"].append(dict(file=str(relative), reason="checkpoint bytes differ"))
        elif old.suffix == ".csv":
            report["csv_files"] += 1
            if physical_csv(old) != physical_csv(new):
                report["differences"].append(dict(file=str(relative), reason="physical CSV contents differ"))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = compare(args.baseline, args.current)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["differences"] or result["missing"]:
        raise SystemExit(1)
