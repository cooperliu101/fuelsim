#!/usr/bin/env python3
"""Read native B544 prescribed histories and distinguish bulk/HG force operators.

This is a diagnostic reconstruction, not a production solve or an acceptance
gate. All stress components and reaction forces come from the native CSV files.
No material parameter or hourglass constant is fitted. The geometry helpers are
independent NumPy implementations used by the earlier B532/B533 probes.
"""

import csv
from pathlib import Path
import sys

import numpy as np

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY.parent))
from analyze_b532_candidates import geometry
from analyze_b533_candidates import finite_total_stiffness
from analyze_b544_mechanical_reactions import COORDINATES, ELEMENTS, center_geometry, grouped_metrics


def read_job(prefix, disconnected):
    node_count = 64 if disconnected else 30
    result = {
        "u": np.zeros((10, node_count, 3)),
        "rf": np.zeros((10, node_count, 3)),
        "nodal_temperature": np.zeros((10, node_count)),
        "stress": np.zeros((10, 8, 3, 3)),
        "temperature": np.zeros((10, 8)),
        "volume": np.zeros((10, 8)),
        "energy": np.zeros(10),
    }
    seen_nodes, seen_points, seen_energy = set(), set(), set()
    for row in csv.DictReader(open(str(prefix) + "_nodal.csv", newline="")):
        step, node = int(row["increment"]) - 1, int(row["node"]) - 1
        assert (step, node) not in seen_nodes
        seen_nodes.add((step, node))
        assert float(row["time_s"]) == (step + 1) * 100000
        if disconnected:
            element, local_node = int(row["element"]) - 1, int(row["local_node"]) - 1
            assert node == element * 8 + local_node
            assert int(row["source_node"]) == ELEMENTS[element][local_node] + 1
        result["u"][step, node] = [float(row["u%d_m" % component]) for component in range(1, 4)]
        result["rf"][step, node] = [float(row["rf%d_n" % component]) for component in range(1, 4)]
        result["nodal_temperature"][step, node] = float(row["temperature_k"])
    for row in csv.DictReader(open(str(prefix) + "_integration.csv", newline="")):
        step, element = int(row["increment"]) - 1, int(row["element"]) - 1
        assert (step, element) not in seen_points and int(row["integration_point"]) == 1
        seen_points.add((step, element))
        assert float(row["time_s"]) == (step + 1) * 100000
        if disconnected:
            assert int(row["local_csys_present"]) == 0
        for i in range(3):
            for j in range(3):
                a, b = sorted((i + 1, j + 1))
                result["stress"][step, element, i, j] = float(row["s%d%d_pa" % (a, b)])
        result["temperature"][step, element] = float(row["temperature_k"])
        result["volume"][step, element] = float(row["ivol_m3"])
    for row in csv.DictReader(open(str(prefix) + "_energy.csv", newline="")):
        step = int(row["increment"]) - 1
        assert step not in seen_energy
        seen_energy.add(step)
        result["energy"][step] = float(row["allae_j"])
    assert seen_nodes == {(step, node) for step in range(10) for node in range(node_count)}
    assert seen_points == {(step, element) for step in range(10) for element in range(8)}
    assert seen_energy == set(range(10))
    assert all(np.all(np.isfinite(value)) for value in result.values())
    return result


def operators(reference, displacement, previous, stress):
    current = reference + displacement
    volume, current_gradient, _, gamma = geometry(current)
    reference_volume, reference_gradient, _, reference_gamma = geometry(reference)
    deformation = np.eye(3) + displacement.T @ reference_gradient
    pushed_gradient = reference_gradient @ np.linalg.inv(deformation)
    _, center_gradient = center_geometry(current)
    midpoint = reference + 0.5 * (displacement + previous)
    _, midpoint_gradient, _, _ = geometry(midpoint)
    midpoint_deformation = np.eye(3) + (current - midpoint).T @ midpoint_gradient
    midpoint_pushed_gradient = midpoint_gradient @ np.linalg.inv(midpoint_deformation)
    pressure = np.trace(stress) / 3
    deviator = stress - pressure * np.eye(3)
    forces = {
        name: volume * (gradient @ deviator.T + current_gradient * pressure)
        for name, gradient in (
            ("current", current_gradient),
            ("center", center_gradient),
            ("reference_push", pushed_gradient),
            ("midpoint_push", midpoint_pushed_gradient),
        )
    }
    # The B533 helper includes .005 * (2e11 / 2.5) = 4e8 Pa.
    # The weak input specifies an absolute stiffness of 1 Pa, hence /4e8.
    weak_force = finite_total_stiffness(reference, displacement) / 4e8
    metric = np.sqrt(2) * reference_volume * np.sum(reference_gradient**2, axis=0) / 6
    transported = (reference_gamma.T @ displacement) @ deformation
    weak_energy = 0.5 * np.sum(metric * transported**2)
    return forces, weak_force, weak_energy, volume, current_gradient, gamma


def relative_l2(actual, expected):
    return 100 * np.linalg.norm(actual - expected) / np.linalg.norm(expected)


def analyze_disconnected(job):
    data = read_job(DIRECTORY / job, True)
    forces = {name: [] for name in ("current", "center", "reference_push", "midpoint_push")}
    bulk, gaps, affine_errors, modal_errors = [], [], [], []
    energies = np.zeros(10)
    for step in range(10):
        for element, connectivity in enumerate(ELEMENTS):
            indices = slice(8 * element, 8 * element + 8)
            displacement = data["u"][step, indices]
            previous = data["u"][step - 1, indices] if step else np.zeros((8, 3))
            reference = COORDINATES[connectivity]
            stress = data["stress"][step, element]
            candidates, weak, energy, volume, gradient, gamma = operators(reference, displacement, previous, stress)
            assert abs(volume - data["volume"][step, element]) < 1e-12
            native_bulk = data["rf"][step, indices] - weak
            bulk.append(native_bulk)
            gap = candidates["current"] - native_bulk
            gaps.append(gap)
            energies[step] += energy
            affine_stress = native_bulk.T @ (reference + displacement) / volume
            affine_errors.append(relative_l2(affine_stress, stress))
            coefficients = np.linalg.lstsq(gamma, gap, rcond=None)[0]
            modal_errors.append(relative_l2(gamma @ coefficients, gap))
            for name, candidate in candidates.items():
                forces[name].append(candidate)
    bulk = np.array(bulk)
    gap = np.array(gaps).reshape(10, 8, 8, 3)
    print("\njob=", job, "coverage=640 nodal / 80 integration / 10 energy samples; local stress frames absent")
    print("nodal_temperature_range_K=", data["nodal_temperature"].min(), data["nodal_temperature"].max())
    print("material_temperature_range_K=", data["temperature"].min(), data["temperature"].max())
    for name, values in forces.items():
        values = np.array(values)
        print("deviatoric_gradient=", name, "relative_L2_percent=", relative_l2(values, bulk),
              "maximum_absolute_force_error_N=", np.max(abs(values - bulk)))
    print("maximum_affine_stress_error_percent=", max(affine_errors))
    print("maximum_current_hourglass_subspace_residual_percent=", max(modal_errors))
    print("weak_energy_relative_error_percent_each_increment=", 100 * (energies / data["energy"] - 1))
    return data, gap


def analyze_connected():
    source = read_job(DIRECTORY.parent / "b544_hex8_c3d8rt_distorted_bending", False)
    weak = read_job(DIRECTORY / "b544_prescribed_weak", False)
    prediction = np.zeros_like(source["rf"])
    weak_prediction = np.zeros_like(source["rf"])
    for step in range(10):
        for element, connectivity in enumerate(ELEMENTS):
            displacement = source["u"][step, connectivity]
            previous = source["u"][step - 1, connectivity] if step else np.zeros((8, 3))
            candidates, weak_force, *_ = operators(
                COORDINATES[connectivity], displacement, previous, source["stress"][step, element]
            )
            # B544 default stiffness is .005 * G_initial = .005 * 8e7 = 4e5 Pa.
            prediction[step, connectivity] += candidates["reference_push"] + weak_force * 4e5
            weak_prediction[step, connectivity] += candidates["reference_push"] + weak_force
    supports = [0, 5, 10, 15, 20, 25]
    print("\noriginal_60_support_vectors_corrected_L2_peak_pointwise_percent=",
          [100 * value for value in grouped_metrics(prediction[:, supports], source["rf"][:, supports])])
    print("prescribed_300_weak_vectors_corrected_L2_peak_pointwise_percent=",
          [100 * value for value in grouped_metrics(weak_prediction, weak["rf"])])


def main():
    original, original_gap = analyze_disconnected("b544_disconnected_weak")
    isothermal, isothermal_gap = analyze_disconnected("b544_disconnected_weak_isothermal")
    print("\noriginal_vs_isothermal_maximum_U_difference_m=", np.max(abs(original["u"] - isothermal["u"])))
    print("original_vs_isothermal_maximum_gap_difference_N=", np.max(abs(original_gap - isothermal_gap)))
    # Actual native elastic table: E=2e8 at 300 K, E=1.2e8 at 700 K, nu=.25.
    # Endpoint holding is independently visible in the first native S/EE values.
    shear_ratio = np.clip(1 - 0.001 * (original["temperature"] - 300), 0.6, 1)
    scaled_gap = shear_ratio[:, :, None, None] * isothermal_gap
    print("current_G_ratio_scaled_gap_relative_L2_percent=", relative_l2(scaled_gap, original_gap))
    print("current_G_ratio_scaled_gap_maximum_difference_N=", np.max(abs(scaled_gap - original_gap)))
    print("Scope: these fixed U/T histories identify the force operator; this script does not establish")
    print("production solution accuracy, material-update equivalence, or arbitrary-loading acceptance.")
    analyze_connected()


if __name__ == "__main__":
    main()
