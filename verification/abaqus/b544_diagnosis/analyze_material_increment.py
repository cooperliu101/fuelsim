#!/usr/bin/env python3
"""Identify B544 material increments from all 80 native elastic histories.

This reads native U/T/S/EE and performs algebra only. It neither calls Fuelsim
nor fits material parameters. EE engineering shear components are divided by
two before tensor operations. Reported tensor metrics use the Frobenius norm;
they are diagnostics, not replacements for the production acceptance gates.

The native card has constant isotropic alpha=1e-5/K and nu=.25, with the stated
two-row E(T) table. This probe does not identify an arbitrary nonlinear-alpha
averaging convention or validate inelastic material integration by itself.
"""

import csv
from pathlib import Path

import numpy as np

from analyze_bulk_split import COORDINATES, ELEMENTS, geometry, read_job

DIRECTORY = Path(__file__).resolve().parent
ALPHA = 1e-5


def rotation(gradient):
    half_spin = 0.25 * (gradient - gradient.T)
    return (np.eye(3) + half_spin) @ np.linalg.inv(np.eye(3) - half_spin)


def deviator(tensor):
    return tensor - np.trace(tensor) * np.eye(3) / 3


def elasticity(temperature):
    # Native table holds the endpoint values outside [300,700] K.
    return np.clip(2e8 - 2e5 * (temperature - 300), 1.2e8, 2e8)


def hooke(elastic_strain, young_modulus):
    return 0.8 * young_modulus * elastic_strain + 0.4 * young_modulus * np.trace(elastic_strain) * np.eye(3)


def metrics(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    differences = np.linalg.norm(actual - expected, axis=(-2, -1))
    reference = np.linalg.norm(expected, axis=(-2, -1))
    nonzero = reference != 0
    return [
        float(100 * np.linalg.norm(differences) / np.linalg.norm(reference)),
        float(100 * differences.max() / reference.max()),
        float(100 * np.max(differences[nonzero] / reference[nonzero])),
    ]


def main():
    prefix = DIRECTORY.parent / "b544_hex8_c3d8rt_distorted_bending"
    data = read_job(prefix, False)
    displacement = np.concatenate((np.zeros((1, 30, 3)), data["u"]))
    temperature = np.concatenate((np.full((1, 8), 300.0), data["temperature"]))
    nodal_temperature = np.concatenate((np.full((1, 30), 300.0), data["nodal_temperature"]))
    stress = np.concatenate((np.zeros((1, 8, 3, 3)), data["stress"]))
    elastic = np.zeros_like(stress)
    seen = set()
    with open(str(prefix) + "_integration.csv", newline="") as stream:
        for row in csv.DictReader(stream):
            step, element = int(row["increment"]), int(row["element"]) - 1
            assert (step, element) not in seen
            seen.add((step, element))
            for i in range(3):
                for j in range(3):
                    a, b = sorted((i + 1, j + 1))
                    suffix = "" if a == b else "_engineering"
                    elastic[step, element, i, j] = float(row["ee%d%d%s" % (a, b, suffix)]) / (1 if a == b else 2)
    assert seen == {(step, element) for step in range(1, 11) for element in range(8)}
    assert np.all(np.isfinite(elastic))

    names = ("midpoint", "reference_push", "split", "split_current_thermal")
    predicted = {name: np.zeros_like(elastic) for name in names}
    needed = {name: [] for name in names}
    increments = {name: [] for name in names}
    recovered_modulus_errors, hooke_stress = [], []
    temperature_errors, theta_errors = [], []
    for step in range(1, 11):
        for element, connectivity in enumerate(ELEMENTS):
            reference = COORDINATES[connectivity]
            current_u = displacement[step, connectivity]
            old_u = displacement[step - 1, connectivity]
            increment = current_u - old_u
            current_volume, _, current_weights, _ = geometry(reference + current_u)
            _, midpoint_gradient, _, _ = geometry(reference + 0.5 * (current_u + old_u))
            _, reference_gradient, _, _ = geometry(reference)
            midpoint_deformation = np.eye(3) + 0.5 * (current_u + old_u).T @ reference_gradient
            pushed_gradient = reference_gradient @ np.linalg.inv(midpoint_deformation)
            midpoint_increment = increment.T @ midpoint_gradient
            pushed_increment = increment.T @ pushed_gradient
            split_increment = pushed_increment + np.eye(3) * (
                np.trace(midpoint_increment) - np.trace(pushed_increment)
            ) / 3

            material_temperature = temperature[step, element]
            temperature_errors.append(current_weights @ nodal_temperature[step, connectivity] / current_volume
                                      - material_temperature)
            young_modulus = elasticity(material_temperature)
            unit_stress = hooke(elastic[step, element], 1.0)
            recovered_modulus = np.sum(unit_stress * stress[step, element]) / np.sum(unit_stress**2)
            recovered_modulus_errors.append(recovered_modulus / young_modulus - 1)
            hooke_stress.append(hooke(elastic[step, element], young_modulus))

            original_temperature_increment = material_temperature - temperature[step - 1, element]
            nodal_temperature_increment = nodal_temperature[step, connectivity] - nodal_temperature[step - 1, connectivity]
            current_temperature_increment = current_weights @ nodal_temperature_increment / current_volume
            needed_temperature_increment = (
                np.trace(midpoint_increment) - np.trace(elastic[step, element] - elastic[step - 1, element])
            ) / (3 * ALPHA)
            theta_errors.append((original_temperature_increment - needed_temperature_increment,
                                 current_temperature_increment - needed_temperature_increment))
            for name, gradient in (("midpoint", midpoint_increment),
                                   ("reference_push", pushed_increment),
                                   ("split", split_increment),
                                   ("split_current_thermal", split_increment)):
                objective_rotation = rotation(gradient)
                spatial_increment = 0.5 * (gradient + gradient.T)
                thermal_temperature_increment = (current_temperature_increment if name == "split_current_thermal"
                                                 else original_temperature_increment)
                thermal = ALPHA * thermal_temperature_increment * np.eye(3)
                native_increment = (elastic[step, element]
                                    - objective_rotation @ elastic[step - 1, element] @ objective_rotation.T + thermal)
                needed[name].append(native_increment)
                increments[name].append(spatial_increment)
                predicted[name][step, element] = (
                    objective_rotation @ predicted[name][step - 1, element] @ objective_rotation.T
                    + spatial_increment - thermal
                )

    print("coverage=10 increments x 8 elements; each native S and EE has 6 independent tensor components")
    print("recovered_E_vs_table_maximum_relative_percent=", 100 * np.max(abs(np.asarray(recovered_modulus_errors))))
    print("native_S_from_EE_L2_peak_pointwise_percent=", metrics(hooke_stress, stress[1:].reshape(-1, 3, 3)))
    print("current_weighted_T_vs_native_TEMP_maximum_error_K=", np.max(abs(np.asarray(temperature_errors))))
    print("thermal_increment_old_vs_current_weights_maximum_error_K=", np.max(abs(np.asarray(theta_errors)), axis=0))
    for name in names:
        difference = np.asarray(increments[name]) - np.asarray(needed[name])
        accumulated_stress = [hooke(predicted[name][step, element], elasticity(temperature[step, element]))
                              for step in range(1, 11) for element in range(8)]
        print("\ncandidate=", name)
        print("increment_L2_peak_pointwise_percent=", metrics(increments[name], needed[name]))
        print("maximum_increment_tensor_difference=", np.max(np.linalg.norm(difference, axis=(-2, -1))))
        print("maximum_increment_deviator_difference=", max(np.linalg.norm(deviator(value)) for value in difference))
        print("maximum_increment_trace_difference=", max(abs(np.trace(value)) for value in difference))
        print("cumulative_EE_L2_peak_pointwise_percent=", metrics(predicted[name][1:], elastic[1:]))
        print("cumulative_S_L2_peak_pointwise_percent=", metrics(accumulated_stress, stress[1:].reshape(-1, 3, 3)))
    print("\nScope: constant isotropic alpha, native E(T) endpoint holding, and these prescribed U/T histories.")
    print("This is independent material-increment evidence, not a free-field production solve or a nonlinear-alpha probe.")


if __name__ == "__main__":
    main()
