"""Read-only checks of complete cards executed by fuelsim -i; no problem generation."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import netCDF4
import numpy as np
from eccentric import verify_eccentric

p = argparse.ArgumentParser()
p.add_argument("--fuelsim", required=True)
p.add_argument("--source", type=Path, required=True)
p.add_argument("--work", type=Path, required=True)
p.add_argument("--study", choices=["single", "coupled"], required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
for directory in ["cases", "models"]:
    shutil.copytree(a.source / directory, a.work / directory, dirs_exist_ok=True)
exe = str(Path(a.fuelsim).resolve())
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
volume = 2 * np.sqrt(2) * 0.004**2 * 0.008
power = 2e8 * volume
surface = np.loadtxt(a.source / "cases/surface_nodes.txt")[:, 0].astype(int)


def run(name):
    card = a.work / "cases" / (name + ".fsi")
    if card.read_bytes() != (a.source / "cases" / card.name).read_bytes():
        raise AssertionError("Input card changed")
    result = subprocess.run([exe, "-i", str(card)], text=True, capture_output=True, env=env)
    (a.work / (name + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    status = dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line)
    if status.get("completed") != "true":
        raise AssertionError(status)
    with netCDF4.Dataset(a.work / "cases" / (name + "_results.e")) as data:
        names = netCDF4.chartostring(data["name_nod_var"][:]).tolist()
        values = {name: np.asarray(data[f"vals_nod_var{i+1}"][-1]) for i, name in enumerate(names)}
    retained = np.arange(len(values["temperature"]))
    if not name.endswith("_full"):
        retained = retained[retained != 13]
    for field in ["temperature", "heat_reaction"]:
        if not np.isfinite(values[field][retained]).all():
            raise AssertionError("Nonfinite retained production field: " + field)
    return values, status


prefix = "pellet" if a.study == "single" else "coupled"
full, fs = run(prefix + "_full")
exact, es = run(prefix + "_exact")
surrogate, ss = run(prefix + "_surrogate")
if (int(fs["global_state_dofs"]), int(es["global_state_dofs"]), int(ss["global_state_dofs"])) != (
    (27, 26, 26) if a.study == "single" else (75, 74, 74)
):
    raise AssertionError("Interior temperature DOF was not eliminated")
for fields in [exact, surrogate]:
    if not np.isnan(fields["temperature"][13]):
        raise AssertionError("Eliminated temperature must not masquerade as a solved field")
report = {}
if a.study == "single":
    for field_case in ["uniform", "mixed"]:
        if field_case == "mixed":
            full, fs = run("pellet_mixed_full")
            exact, es = run("pellet_mixed_exact")
            surrogate, ss = run("pellet_mixed_surrogate")
        ref = full["heat_reaction"][surface]
        for name, fields, tolerance in [("exact", exact, 1e-10), ("surrogate", surrogate, 1e-3)]:
            r = fields["heat_reaction"][surface]
            error = float(abs(r - ref).max())
            balance = float(abs(r.sum() + power))
            report[field_case + "_" + name] = {
                "residual_max_error_W": error,
                "energy_error_W": balance,
                "residual_relative_l2": float(np.linalg.norm(r - ref) / np.linalg.norm(ref)),
            }
            if error > tolerance or balance > tolerance:
                raise AssertionError(report)

else:

    def cooling(t):
        heat = 0.0
        for z in range(2):
            for i in range(8):
                j = (i + 1) % 8
                nodes = [
                    27 + 16 * z + 8 + i,
                    27 + 16 * z + 8 + j,
                    27 + 16 * (z + 1) + 8 + i,
                    27 + 16 * (z + 1) + 8 + j,
                ]
                area = 2 * 0.0047 * np.sin(np.pi / 8) * 0.004
                heat += 10000 * area * (t[nodes].mean() - 550)
        return heat

    for name, fields, status, tolerance in [("exact", exact, es, 1e-9), ("surrogate", surrogate, ss, 1e-3)]:
        tp = float(abs(fields["temperature"][surface] - full["temperature"][surface]).max())
        tc = float(abs(fields["temperature"][27:] - full["temperature"][27:]).max())
        gap = float(status["contact.gap.total_heat_rate"])
        gap_ref = float(fs["contact.gap.total_heat_rate"])
        energy = abs(cooling(fields["temperature"]) - power)
        report[name] = {
            "pellet_surface_temperature_error_K": tp,
            "cladding_temperature_error_K": tc,
            "gap_heat_rate_W": gap,
            "gap_heat_error_W": abs(gap - gap_ref),
            "energy_error_W": energy,
            "nonlinear_iterations": int(status["nonlinear_iterations_total"]),
        }
        if tp > tolerance or tc > tolerance or abs(gap - gap_ref) > tolerance or energy > tolerance:
            raise AssertionError(report)
    # Exercise the complete global derivative including both sides of the gap.
    result = subprocess.run(
        [exe, "-i", str(a.work / "cases/coupled_surrogate.fsi"), "--check-jacobian"],
        capture_output=True,
        text=True,
        env=env,
    )
    (a.work / "global_jacobian.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    report["eccentric"] = verify_eccentric(run, cooling, surface, power, a.work, exe, env)
(a.work / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
