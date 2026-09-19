"""Read-only qualification of eight complete full/surrogate production-card pairs."""

import json
import numpy as np


def verify_scan(run, cooling, surface, work):
    training = json.loads((work / "models/pellet_mlp.training.json").read_text())
    volume = 2 * np.sqrt(2) * 0.004**2 * 0.008
    report = {"cases": {}, "temperature_tolerance_K": 1e-3, "heat_rate_tolerance_W": 1e-3}
    failed = []
    for eccentricity in [0, 75]:
        for source_name, source in [("low", 5e7), ("high", 3.5e8)]:
            for gas_name, gas in [("low", 0.25), ("high", 1.0)]:
                name = f"scan_{eccentricity}_{source_name}_{gas_name}"
                entry = {
                    "eccentricity_um": eccentricity,
                    "heat_source_W_m3": source,
                    "gas_conductivity_W_m_K": gas,
                }
                report["cases"][name] = entry
                try:
                    full, fs = run(name + "_full")
                    surrogate, ss = run(name + "_surrogate")
                    power = source * volume
                    q_full = np.array([float(fs[f"contact.sector_{i}.total_heat_rate"]) for i in range(8)])
                    q_surrogate = np.array(
                        [float(ss[f"contact.sector_{i}.total_heat_rate"]) for i in range(8)]
                    )
                    gaps = np.array([float(fs[f"contact.sector_{i}.minimum_gap"]) for i in range(8)])
                    if (
                        not np.isfinite(q_full).all()
                        or not np.isfinite(q_surrogate).all()
                        or not np.isfinite(gaps).all()
                    ):
                        raise AssertionError("Nonfinite contact diagnostic")
                    if int(fs["global_state_dofs"]) != 75 or int(ss["global_state_dofs"]) != 74:
                        raise AssertionError("Interior temperature DOF was not eliminated")
                    if not np.isnan(surrogate["temperature"][13]):
                        raise AssertionError("Surrogate unexpectedly reconstructed an interior temperature")
                    expected_gap = 100e-6 * np.cos(np.pi / 8) - eccentricity * 1e-6 * np.cos(
                        np.pi / 8 + np.arange(8) * np.pi / 4
                    )
                    if not np.allclose(gaps, expected_gap, atol=1e-15, rtol=0) or np.min(gaps) <= 1e-6:
                        raise AssertionError("Unexpected scan geometry or active minimum-gap floor")
                    for i in range(8):
                        if float(ss[f"contact.sector_{i}.minimum_gap"]) != gaps[i]:
                            raise AssertionError("Full and surrogate gap geometries differ")
                        for status in [fs, ss]:
                            if int(status[f"contact.sector_{i}.projected_contact_nodes"]) != 18:
                                raise AssertionError("Missing contact integration points")
                            if int(status[f"contact.sector_{i}.unprojected_contact_nodes"]) != 0:
                                raise AssertionError("Unprojected contact integration points")
                            if (
                                int(status["load_cutbacks"]) != 0
                                or status["used_backtracking_fallback"] != "false"
                            ):
                                raise AssertionError("Unexpected nonlinear retry or fallback")
                    full_energy = float(abs(cooling(full["temperature"]) - power))
                    if full_energy > 1e-9 or abs(q_full.sum() - power) > 1e-9:
                        raise AssertionError("Full FEM energy balance failed")
                    temperature_min = float(full["temperature"][surface].min())
                    temperature_max = float(full["temperature"][surface].max())
                    entry.update(
                        {
                            "full_surface_temperature_min_K": temperature_min,
                            "full_surface_temperature_max_K": temperature_max,
                            "within_training_temperature_envelope": bool(
                                temperature_min >= training["training_temperature_min_K"]
                                and temperature_max <= training["training_temperature_max_K"]
                            ),
                            "pellet_temperature_error_K": float(
                                abs(surrogate["temperature"][surface] - full["temperature"][surface]).max()
                            ),
                            "cladding_temperature_error_K": float(
                                abs(surrogate["temperature"][27:] - full["temperature"][27:]).max()
                            ),
                            "sector_heat_error_W": float(abs(q_surrogate - q_full).max()),
                            "sector_heat_relative_error": float(abs((q_surrogate - q_full) / q_full).max()),
                            "total_gap_heat_error_W": float(abs(q_surrogate.sum() - q_full.sum())),
                            "energy_error_W": float(abs(cooling(surrogate["temperature"]) - power)),
                            "full_energy_error_W": full_energy,
                            "full_sector_heat_W": q_full.tolist(),
                            "surrogate_sector_heat_W": q_surrogate.tolist(),
                            "minimum_gap_m": float(gaps.min()),
                            "full_nonlinear_iterations": int(fs["nonlinear_iterations_total"]),
                            "surrogate_nonlinear_iterations": int(ss["nonlinear_iterations_total"]),
                        }
                    )
                    entry["passed"] = all(
                        entry[key] <= 1e-3
                        for key in [
                            "pellet_temperature_error_K",
                            "cladding_temperature_error_K",
                            "sector_heat_error_W",
                            "total_gap_heat_error_W",
                            "energy_error_W",
                        ]
                    )
                    if not entry["passed"]:
                        failed.append(name)
                except (AssertionError, KeyError, ValueError) as error:
                    entry["passed"] = False
                    entry["error"] = str(error)
                    failed.append(name)
    report["passed"] = not failed
    report["failed_cases"] = failed
    (work / "scan_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    if failed:
        raise AssertionError("Parameter-scan failures: " + ", ".join(failed) + "; see scan_metrics.json")
    return report
