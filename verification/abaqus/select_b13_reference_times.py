#!/usr/bin/env python3
"""Select complete native CSV frames without rewriting any retained row."""

import argparse
import collections
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument("--destination-directory", type=Path, required=True)
    args = parser.parse_args()
    selection = Path(__file__).with_name("b13_integrated_reference_times.txt")
    times = [float(value) for value in selection.read_text().split()]
    if times != sorted(set(times)) or times[0] != 0.0625 or times[-1] != 6.0:
        raise ValueError("Invalid reference time selection")
    expected = {step / 16 for step in range(1, 97)}
    if not set(times) <= expected:
        raise ValueError("Selected times are outside the native time grid")
    prepared = []
    for element in ("cax4t", "cax4rt", "cax8t", "cax8rt"):
        for field in ("nodes", "points", "contact"):
            name = f"b13_integrated_{element}_{field}.csv"
            source = args.source_directory / name
            destination = args.destination_directory / name
            if source.resolve() == destination.resolve():
                raise ValueError("Source must be a separate full-history archive")
            original = source.read_bytes()
            lines = original.splitlines(keepends=True)
            if not lines or not lines[0].startswith(b"time,"):
                raise ValueError(f"Missing time column: {source}")
            counts = collections.Counter()
            retained = [lines[0]]
            previous = 0.0
            for row in lines[1:]:
                time = float(row.split(b",", 1)[0])
                if time < previous:
                    raise ValueError(f"Unsorted source times: {source}")
                previous = time
                counts[time] += 1
                if time in times:
                    retained.append(row)
            if set(counts) != expected or len(set(counts.values())) != 1:
                raise ValueError(f"Expected 96 complete, equally sized frames: {source}")
            reduced = b"".join(retained)
            prepared.append((name, original, reduced, len(lines) - 1, len(retained) - 1))

    # Validate all sources before writing any destination. Row bytes, numeric
    # precision, column order and spatial sampling remain exactly as extracted.
    args.destination_directory.mkdir(parents=True, exist_ok=True)
    report = ["file\toriginal_sha256\tselected_sha256\toriginal_bytes\tselected_bytes"
              "\toriginal_rows\tselected_rows\toriginal_frames\tselected_frames\n"]
    for name, original, reduced, original_rows, selected_rows in prepared:
        (args.destination_directory / name).write_bytes(reduced)
        report.append(f"{name}\t{hashlib.sha256(original).hexdigest()}\t"
                      f"{hashlib.sha256(reduced).hexdigest()}\t{len(original)}\t{len(reduced)}\t"
                      f"{original_rows}\t{selected_rows}\t96\t{len(times)}\n")
    (args.destination_directory / "b13_integrated_reference_selection.tsv").write_text("".join(report))
    before = sum(len(entry[1]) for entry in prepared)
    after = sum(len(entry[2]) for entry in prepared)
    print(f"Selected {len(times)}/96 frames in {len(prepared)} CSV files: "
          f"{before} -> {after} bytes ({100 * (1 - after / before):.3f}% smaller)")


if __name__ == "__main__":
    main()
