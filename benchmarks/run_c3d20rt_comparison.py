#!/usr/bin/env python3
"""Compare two production executables using unchanged, isolated C3D20RT cases."""

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--case", choices=("finite_coupled", "small_coupled", "finite_contact_first_increment"))
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    root = Path(__file__).resolve().parents[1]
    executables = {"before": args.before.resolve(), "after": args.after.resolve()}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cases = {
        "finite_coupled": "c3d20rt_probe.e",
        "small_coupled": "c3d20rt_probe.e",
        "finite_contact_first_increment": "c3d20rt_finite_contact.e",
    }
    if args.case:
        cases = {args.case: cases[args.case]}
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    rows = []
    with (output / "executables.txt").open("w") as manifest:
        for label, executable in executables.items():
            digest = hashlib.sha256(executable.read_bytes()).hexdigest()
            manifest.write(f"{label} {digest} {executable}\n")
    with (output / "inputs.sha256").open("w") as manifest:
        for case, mesh in cases.items():
            for relative in (f"verification/fuelsim/transient_c3d20rt_{case}.fsi",
                             f"verification/meshes/{mesh}"):
                digest = hashlib.sha256((root / relative).read_bytes()).hexdigest()
                manifest.write(f"{digest}  {relative}\n")
    for case, mesh in cases.items():
        filename = f"transient_c3d20rt_{case}.fsi"
        for label in executables:
            directory = output / case / label
            (directory / "fuelsim").mkdir(parents=True)
            (directory / "meshes").mkdir()
            shutil.copyfile(root / "verification/fuelsim" / filename, directory / "fuelsim" / filename)
            shutil.copyfile(root / "verification/meshes" / mesh, directory / "meshes" / mesh)
        for repeat in range(-1, args.repeats):
            labels = ["before", "after"] if repeat % 2 == 0 else ["after", "before"]
            for label in labels:
                directory = output / case / label
                command = ["taskset", "-c", str(args.cpu), str(executables[label]),
                           "-i", str(directory / "fuelsim" / filename)]
                log = directory / f"run-{repeat}.log"
                with log.open("w") as stream:
                    start = time.perf_counter()
                    subprocess.run(command, env=environment, stdout=stream, stderr=subprocess.STDOUT, check=True)
                    elapsed = time.perf_counter() - start
                values = dict(line.split("=", 1) for line in log.read_text().splitlines() if "=" in line)
                row = {"case": case, "version": label, "repeat": repeat, "wall_seconds": elapsed}
                for key in ("total_seconds", "accepted_steps", "rejected_steps", "total_cutbacks",
                            "nonlinear_iterations_total", "linear_iterations_total",
                            "residual_evaluations_total", "jacobian_evaluations_total"):
                    row[key] = values[key]
                rows.append(row)
                with (output / "runs.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row), lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
                print(f"{case} {label} repeat={repeat} wall={elapsed:.6f}s", flush=True)
        medians = {label: statistics.median(row["wall_seconds"] for row in rows
                   if row["case"] == case and row["version"] == label and row["repeat"] >= 0)
                   for label in executables}
        print(f"{case}: median before={medians['before']:.6f}s after={medians['after']:.6f}s "
              f"speedup={medians['before'] / medians['after']:.4f}", flush=True)


if __name__ == "__main__":
    main()
