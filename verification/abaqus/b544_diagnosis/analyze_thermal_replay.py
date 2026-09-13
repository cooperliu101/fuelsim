#!/usr/bin/env python3
"""Independently reconstruct B544 prescribed-history heat reactions.

This reads native CSVs and uses the independent B532 NumPy geometry rules.
It does not call Fuelsim, solve an equation, fit a coefficient, or change inputs.
Total storage alone does not identify each nodal capacity weight. When present,
the separate 1e-20 W/(m K) conductivity job supplies that nodal isolation.
"""

import csv
from pathlib import Path
import sys

import numpy as np

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY.parent))
from analyze_b532_candidates import geometry, thermal_coefficients
from analyze_b544_mechanical_reactions import COORDINATES, ELEMENTS, center_geometry
from analyze_bulk_split import read_job


def read_heat_fields(prefix):
    data = read_job(prefix, False)
    data["rfl"] = np.zeros((10, 30))
    data["hfl"] = np.zeros((10, 8, 3))
    seen_nodes, seen_points = set(), set()
    with open(str(prefix) + "_nodal.csv", newline="") as stream:
        for row in csv.DictReader(stream):
            key = int(row["increment"]) - 1, int(row["node"]) - 1
            assert key not in seen_nodes
            seen_nodes.add(key)
            data["rfl"][key] = float(row["reaction_heat_flux_w"])
    with open(str(prefix) + "_integration.csv", newline="") as stream:
        for row in csv.DictReader(stream):
            key = int(row["increment"]) - 1, int(row["element"]) - 1
            assert key not in seen_points
            seen_points.add(key)
            data["hfl"][key] = [float(row[f"hfl{i}_w_m2"]) for i in (1, 2, 3)]
    assert seen_nodes == {(step, node) for step in range(10) for node in range(30)}
    assert seen_points == {(step, element) for step in range(10) for element in range(8)}
    assert all(np.all(np.isfinite(value)) for value in data.values())
    return data


def scalar_metrics(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    difference = actual - expected
    nonzero = expected != 0.0
    return {
        "L2_percent": 100 * np.linalg.norm(difference) / np.linalg.norm(expected),
        "peak_percent": 100 * abs(np.max(abs(actual)) - np.max(abs(expected))) / np.max(abs(expected)),
        "pointwise_percent": 100 * np.max(abs(difference[nonzero] / expected[nonzero])),
        "maximum_absolute": np.max(abs(difference)),
        "zero_count": int(np.sum(~nonzero)),
        "zero_absolute": np.max(abs(difference[~nonzero])) if np.any(~nonzero) else 0.0,
    }


def report(label, actual, expected):
    print(label, " ".join(f"{key}={value:.12g}" for key, value in scalar_metrics(actual, expected).items()))


def candidates(data):
    capacity_names = ("current", "reference", "old", "midpoint", "center", "current_equal", "linear_cp")
    capacity = {name: np.zeros((10, 30)) for name in capacity_names}
    conduction = {}
    flux = {}
    uniform = np.zeros((10, 30))
    hourglass = np.zeros((10, 30))
    weighted_temperature = np.zeros((10, 8))
    temperature_names = ("current", "reference", "midpoint", "arithmetic", "native")
    for configuration in ("current", "midpoint", "reference"):
        for temperature_name in temperature_names:
            name = f"{configuration}/{temperature_name}"
            conduction[name] = np.zeros((10, 30))
            flux[name] = np.zeros((10, 8, 3))
    for name in ("linear_k", "uniform_test_reference_push", "hourglass_center_metric", "hourglass_center_measure_metric"):
        conduction[name] = np.zeros((10, 30))

    for step in range(10):
        for element, nodes in enumerate(ELEMENTS):
            reference = COORDINATES[nodes]
            displacement = data["u"][step, nodes]
            old_displacement = data["u"][step - 1, nodes] if step else np.zeros((8, 3))
            temperature = data["nodal_temperature"][step, nodes]
            old_temperature = data["nodal_temperature"][step - 1, nodes] if step else np.full(8, 300.0)
            current = reference + displacement
            geometries = {
                name: geometry(coordinates)
                for name, coordinates in (
                    ("current", current),
                    ("reference", reference),
                    ("old", reference + old_displacement),
                    ("midpoint", reference + 0.5 * (displacement + old_displacement)),
                )
            }
            volume, gradient, weights, gamma = geometries["current"]
            assert abs(volume - data["volume"][step, element]) < 1e-12
            sampled_temperatures = {name: item[2] @ temperature / item[0] for name, item in geometries.items()}
            sampled_temperatures.update(arithmetic=np.mean(temperature), native=data["temperature"][step, element])
            weighted_temperature[step, element] = sampled_temperatures["current"]
            rate = (temperature - old_temperature) / 100000.0
            specific_heat = 500.0 + 1.25 * (np.clip(temperature, 300.0, 700.0) - 300.0)
            for name in capacity_names:
                if name == "center":
                    node_weights = np.full(8, center_geometry(current)[0] / 8)
                elif name == "current_equal":
                    node_weights = np.full(8, volume / 8)
                else:
                    node_weights = geometries["current" if name == "linear_cp" else name][2]
                cp = 500.0 + 1.25 * (temperature - 300.0) if name == "linear_cp" else specific_heat
                capacity[name][step, nodes] += node_weights * 2000.0 * cp * rate

            for configuration in ("current", "midpoint", "reference"):
                candidate_volume, candidate_gradient, _, candidate_gamma = geometries[configuration]
                coefficients = thermal_coefficients(candidate_volume, candidate_gradient)
                unit_uniform = candidate_volume * candidate_gradient @ (candidate_gradient.T @ temperature)
                unit_hourglass = candidate_gamma @ (coefficients * (candidate_gamma.T @ temperature))
                for temperature_name in temperature_names:
                    name = f"{configuration}/{temperature_name}"
                    sample = np.clip(sampled_temperatures[temperature_name], 300.0, 700.0)
                    conductivity = 10.0 + 0.02 * (sample - 300.0)
                    conduction[name][step, nodes] += conductivity * (unit_uniform + unit_hourglass)
                    flux[name][step, element] = -conductivity * candidate_gradient.T @ temperature
                    if name == "current/current":
                        uniform[step, nodes] += conductivity * unit_uniform
                        hourglass[step, nodes] += conductivity * unit_hourglass

            temperature_sample = sampled_temperatures["current"]
            conductivity = 10.0 + 0.02 * (np.clip(temperature_sample, 300.0, 700.0) - 300.0)
            unit_hourglass = gamma @ (thermal_coefficients(volume, gradient) * (gamma.T @ temperature))
            conduction["linear_k"][step, nodes] += (10.0 + 0.02 * (temperature_sample - 300.0)) * (
                volume * gradient @ (gradient.T @ temperature) + unit_hourglass
            )
            reference_gradient = geometries["reference"][1]
            average_deformation = np.eye(3) + displacement.T @ reference_gradient
            pushed_gradient = reference_gradient @ np.linalg.inv(average_deformation)
            conduction["uniform_test_reference_push"][step, nodes] += conductivity * (
                volume * pushed_gradient @ (gradient.T @ temperature) + unit_hourglass
            )
            center_measure, center_gradient = center_geometry(current)
            conduction["hourglass_center_metric"][step, nodes] += conductivity * (
                volume * gradient @ (gradient.T @ temperature)
                + gamma @ (thermal_coefficients(volume, center_gradient) * (gamma.T @ temperature))
            )
            conduction["hourglass_center_measure_metric"][step, nodes] += conductivity * (
                volume * gradient @ (gradient.T @ temperature)
                + gamma @ (thermal_coefficients(center_measure, center_gradient) * (gamma.T @ temperature))
            )
    return capacity, conduction, flux, uniform, hourglass, weighted_temperature


def original_probe_matrices():
    # These independent, constant-k steady probes predate B544. Recover every
    # thermal matrix entry directly from their prescribed +/- nodal temperatures.
    # Modal projection diagnoses the native matrix, without fitting its entries
    # into the geometry candidate used for the separate B544 comparison.
    cases = (
        ("b528_hex8_c3d8rt_operator", False),
        ("b528_hex8_c3d8rt_warped_operator", False),
        ("b528_hex8_c3d8rt_holdout_operator", False),
        ("b532_hex8_c3d8rt_finite_regular_operator", True),
        ("b532_hex8_c3d8rt_finite_warped_operator", True),
    )
    for prefix, finite in cases:
        records, seen = {}, set()
        with open(DIRECTORY.parent / f"{prefix}_nodal.csv", newline="") as stream:
            for row in csv.DictReader(stream):
                case = row.get("case", row.get("step"))
                node = int(row["node"]) - 1
                assert (case, node) not in seen
                seen.add((case, node))
                values = records.setdefault(case, (np.zeros(8), np.zeros((8, 3)), np.zeros(8)))
                values[0][node] = float(row["temperature_k"])
                values[1][node] = [float(row[name]) for name in ("ux_m", "uy_m", "uz_m")]
                values[2][node] = float(row["rfl_w"])
        assert len(records) == 65 and len(seen) == 520
        coordinates = []
        in_nodes = False
        for line in (DIRECTORY.parent / f"{prefix}_probe.inp").read_text().splitlines():
            if line.lower().startswith("*node"):
                in_nodes = True
                continue
            if in_nodes and line.startswith("*"):
                break
            if in_nodes:
                fields = line.split(",")
                if len(coordinates) < 8:
                    assert int(fields[0]) == len(coordinates) + 1
                    coordinates.append([float(value) for value in fields[1:4]])
        assert len(coordinates) == 8
        coordinates = np.array(coordinates)
        if finite:
            coordinates += records["BASE"][1]
        volume, gradient, _, gamma = geometry(coordinates)
        center_measure, center_gradient = center_geometry(coordinates)
        native_matrix = np.stack([
            (records[f"D{column:02d}_PLUS"][2] - records[f"D{column:02d}_MINUS"][2]) / 0.002 / 4.0
            for column in range(8)
        ], axis=1)
        uniform_matrix = volume * gradient @ gradient.T
        old_matrix = uniform_matrix + gamma @ np.diag(thermal_coefficients(volume, gradient)) @ gamma.T
        center_matrix = uniform_matrix + gamma @ np.diag(
            thermal_coefficients(center_measure, center_gradient)
        ) @ gamma.T
        projection = np.linalg.pinv(gamma)
        modal_matrix = projection @ (native_matrix - uniform_matrix) @ projection.T
        outside = native_matrix - uniform_matrix - gamma @ modal_matrix @ gamma.T
        print("probe=", prefix, "finite=", finite,
              "old_thermal_matrix_L2_percent=", 100 * np.linalg.norm(old_matrix - native_matrix) / np.linalg.norm(native_matrix),
              "center_thermal_matrix_L2_percent=", 100 * np.linalg.norm(center_matrix - native_matrix) / np.linalg.norm(native_matrix),
              "outside_modal_L2_percent=", 100 * np.linalg.norm(outside) / np.linalg.norm(native_matrix),
              "modal_offdiagonal_maximum=", np.max(abs(modal_matrix - np.diag(np.diag(modal_matrix)))))
        report(f"{prefix}_center_thermal_residual", 4 * center_matrix @ records["BASE"][0], records["BASE"][2])


def main():
    default = read_heat_fields(DIRECTORY / "b544_prescribed_default")
    weak = read_heat_fields(DIRECTORY / "b544_prescribed_weak")
    print("coverage_per_job=300 nodal / 80 integration / 10 energy records; all 10 accepted increments")
    for field in ("u", "nodal_temperature", "temperature", "volume", "rfl", "hfl"):
        print(f"default_weak_{field}_maximum_difference={np.max(abs(default[field] - weak[field])):.12g}")
    capacity, conduction, flux, uniform, hourglass, temperature = candidates(default)
    print("current_weighted_vs_native_temperature_maximum_K=", np.max(abs(temperature - default["temperature"])))
    print("native_total_RFL_W=", np.sum(default["rfl"], axis=1))
    for name, values in capacity.items():
        print(f"capacity_total_{name}_minus_native_RFL_W=", np.sum(values - default["rfl"], axis=1))
    print("Capacity totals alone do not establish individual nodal capacity weights.")

    left = np.where(COORDINATES[:, 0] == 0.0)[0]
    assert len(left) == 6 and np.all(default["nodal_temperature"][:, left] == 300.0)
    assert np.all(capacity["current"][:, left] == 0.0)
    print("left_60_samples_have_zero_nodal_temperature_rate_and_separate_conduction_only_evidence=True")
    for name, values in conduction.items():
        report(f"RFL_{name}", capacity["current"] + values, default["rfl"])
        report(f"left_conduction_{name}", values[:, left], default["rfl"][:, left])
        if name in flux:
            report(f"HFL_{name}", flux[name], default["hfl"])
    report("literal_production_linear_k_cp_RFL", capacity["linear_cp"] + conduction["linear_k"], default["rfl"])
    report("left_hourglass_after_subtracting_uniform_flux", hourglass[:, left], (default["rfl"] - uniform)[:, left])
    gap = capacity["current"] + conduction["current/current"] - default["rfl"]
    index = np.unravel_index(np.argmax(abs(gap)), gap.shape)
    print("baseline_maximum_error_increment_node=", index[0] + 1, index[1] + 1,
          "capacity_W=", capacity["current"][index], "uniform_W=", uniform[index],
          "hourglass_W=", hourglass[index], "native_RFL_W=", default["rfl"][index])

    isolated_prefix = DIRECTORY / "b544_prescribed_capacity_only"
    if Path(str(isolated_prefix) + "_nodal.csv").exists():
        isolated = read_heat_fields(isolated_prefix)
        for field in ("u", "nodal_temperature", "temperature", "volume"):
            print(f"capacity_job_{field}_maximum_difference={np.max(abs(default[field] - isolated[field])):.12g}")
        isolated_capacity, _, _, _, _, _ = candidates(isolated)
        for name, values in isolated_capacity.items():
            report(f"isolated_nodal_capacity_{name}", values, isolated["rfl"])
        isolated_conduction = default["rfl"] - isolated["rfl"]
        report("independently_isolated_conduction", conduction["current/current"], isolated_conduction)
        report("independently_isolated_hourglass", hourglass, isolated_conduction - uniform)
        report("center_measure_metric_isolated_conduction",
               conduction["hourglass_center_measure_metric"], isolated_conduction)
        print("capacity_job_HFL_maximum_W_per_m2=", np.max(abs(isolated["hfl"])))
    else:
        print("Separate native capacity-only job has not been extracted; nodal capacity identification remains pending.")
    original_probe_matrices()
    print("Scope: independent prescribed-history diagnostic, not a production full-field acceptance gate.")


if __name__ == "__main__":
    main()
