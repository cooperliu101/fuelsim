"""Eccentric gap qualification from production output, without retraining the pellet."""

import json
from pathlib import Path
import subprocess

import netCDF4
import numpy as np


def verify_eccentric(run, cooling, surface, power, work, executable, environment):
    centered, centered_status = run("concentric_gas_full")
    full, full_status = run("eccentric_full")
    exact, exact_status = run("eccentric_exact")
    surrogate, surrogate_status = run("eccentric_surrogate")
    # The pellet is identical to the training geometry. Only the cladding moves.
    with netCDF4.Dataset(work / "cases/concentric_gas.e") as mesh:
        centered_coordinates = np.column_stack([mesh["coord" + axis][:] for axis in "xyz"])
    with netCDF4.Dataset(work / "cases/eccentric.e") as mesh:
        eccentric_coordinates = np.column_stack([mesh["coord" + axis][:] for axis in "xyz"])
    if not np.array_equal(centered_coordinates[:27], eccentric_coordinates[:27]):
        raise AssertionError("The eccentric test changed the trained pellet geometry")
    shift = eccentric_coordinates[27:] - centered_coordinates[27:]
    if not np.allclose(shift, [-50e-6, 0.0, 0.0], rtol=0.0, atol=1e-18):
        raise AssertionError("Incorrect cladding translation")

    def sectors(status):
        heat = np.array([float(status[f"contact.sector_{i}.total_heat_rate"]) for i in range(8)])
        gap = np.array([float(status[f"contact.sector_{i}.minimum_gap"]) for i in range(8)])
        if not np.isfinite(heat).all() or not np.isfinite(gap).all():
            raise AssertionError("Nonfinite sector result")
        for i in range(8):
            if int(status[f"contact.sector_{i}.projected_contact_nodes"]) != 18:
                raise AssertionError("Missing sector integration points")
            if int(status[f"contact.sector_{i}.unprojected_contact_nodes"]) != 0:
                raise AssertionError("Unprojected sector integration points")
        if np.min(gap) <= 1e-6:
            raise AssertionError("The minimum-gap regularization must remain inactive")
        return heat, gap

    q_centered, gap_centered = sectors(centered_status)
    q_full, gap_full = sectors(full_status)
    angles = np.pi / 8 + np.arange(8) * np.pi / 4
    expected_gap = 100e-6 * np.cos(np.pi / 8) - 50e-6 * np.cos(angles)
    if not np.allclose(gap_full, expected_gap, rtol=0.0, atol=1e-15):
        raise AssertionError("Eccentric sector gaps do not match the fixed octagon geometry")
    if not np.allclose(gap_centered, 100e-6 * np.cos(np.pi / 8), rtol=0.0, atol=1e-15):
        raise AssertionError("Concentric gas reference has unequal gaps")
    for fields, heat in [(centered, q_centered), (full, q_full)]:
        if abs(np.sum(heat) - power) > 1e-9 or abs(cooling(fields["temperature"]) - power) > 1e-9:
            raise AssertionError("Full FEM energy balance failed")
    if np.ptp(q_centered) > 1e-9:
        raise AssertionError("The concentric gas reference must have uniform sector heat flow")
    if np.ptp(q_full) / np.mean(q_full) < 0.10 or np.mean(q_full[[0, 7]]) <= np.mean(q_full[[3, 4]]):
        raise AssertionError("Eccentricity did not produce the expected nonuniform gap heat flow")
    ring = np.array([5, 8, 7, 6, 3, 0, 1, 2]) + 9
    if np.ptp(full["temperature"][ring]) < 1.0:
        raise AssertionError("The eccentric pellet surface temperature is unexpectedly uniform")
    centered_opposite = centered["temperature"][12] - centered["temperature"][14]
    eccentric_opposite = full["temperature"][12] - full["temperature"][14]
    if abs(centered_opposite) > 1e-9 or eccentric_opposite < 1.0:
        raise AssertionError("Eccentricity did not produce the expected opposite-side temperature difference")

    report = {
        "pellet_relative_eccentricity_m": 50e-6,
        "sector_angles_degrees": np.rad2deg(angles).tolist(),
        "sector_minimum_gap_m": gap_full.tolist(),
        "concentric_sector_heat_W": q_centered.tolist(),
        "full_sector_heat_W": q_full.tolist(),
        "full_sector_heat_spread_relative_to_mean": float(np.ptp(q_full) / np.mean(q_full)),
        "full_midheight_pellet_temperature_K": full["temperature"][ring].tolist(),
        "full_midheight_pellet_temperature_spread_K": float(np.ptp(full["temperature"][ring])),
        "concentric_midheight_far_minus_near_temperature_K": float(centered_opposite),
        "full_midheight_far_minus_near_temperature_K": float(eccentric_opposite),
        "full_energy_error_W": float(abs(cooling(full["temperature"]) - power)),
        "full_nonlinear_iterations": int(full_status["nonlinear_iterations_total"]),
    }
    if int(full_status["global_state_dofs"]) != 75:
        raise AssertionError("Unexpected full FEM degree-of-freedom count")
    for name, fields, status, tolerance in [
        ("exact", exact, exact_status, 1e-9),
        ("surrogate", surrogate, surrogate_status, 1e-3),
    ]:
        if int(status["global_state_dofs"]) != 74 or not np.isnan(fields["temperature"][13]):
            raise AssertionError("The eccentric pellet interior DOF was not eliminated")
        heat, gap = sectors(status)
        if not np.array_equal(gap, gap_full):
            raise AssertionError("Replacement changed the gap geometry")
        temperature_error = float(
            np.max(np.abs(fields["temperature"][surface] - full["temperature"][surface]))
        )
        cladding_error = float(np.max(np.abs(fields["temperature"][27:] - full["temperature"][27:])))
        sector_error = float(np.max(np.abs(heat - q_full)))
        total_error = float(abs(np.sum(heat) - np.sum(q_full)))
        energy = float(abs(cooling(fields["temperature"]) - power))
        if temperature_error > tolerance or cladding_error > tolerance:
            raise AssertionError("Eccentric temperature comparison failed")
        if sector_error > tolerance or total_error > tolerance or energy > tolerance:
            raise AssertionError("Eccentric heat-flow comparison failed")
        report[name] = {
            "pellet_surface_temperature_error_K": temperature_error,
            "cladding_temperature_error_K": cladding_error,
            "sector_heat_W": heat.tolist(),
            "sector_heat_max_error_W": sector_error,
            "sector_heat_max_relative_error": float(np.max(np.abs((heat - q_full) / q_full))),
            "total_gap_heat_error_W": total_error,
            "energy_error_W": energy,
            "nonlinear_iterations": int(status["nonlinear_iterations_total"]),
        }

    checked = subprocess.run(
        [
            executable,
            "-i",
            str(work / "cases/eccentric_surrogate.fsi"),
            "--check-jacobian",
        ],
        capture_output=True,
        text=True,
        env=environment,
    )
    (work / "eccentric_global_jacobian.log").write_text(checked.stdout + checked.stderr)
    if checked.returncode or "jacobian.check_passed=true" not in checked.stdout:
        raise AssertionError(checked.stdout + checked.stderr)
    status = dict(line.split("=", 1) for line in checked.stdout.splitlines() if "=" in line)
    report["global_tangent_relative_l2_error"] = float(status["jacobian.temperature.relative_l2"])
    (Path(work) / "eccentric_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
