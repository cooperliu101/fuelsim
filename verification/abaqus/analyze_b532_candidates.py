#!/usr/bin/env python3
"""Compare finite C3D8RT probe data with geometry-derived candidate rules."""

import csv
import math
import sys

import numpy as np


SIGNS = np.array(
    [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ]
)
RAW = np.array(
    [
        (1.0, -1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0, 1.0),
        (1.0, 1.0, -1.0, -1.0),
        (-1.0, 1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0, -1.0),
        (1.0, -1.0, 1.0, 1.0),
        (-1.0, -1.0, -1.0, -1.0),
    ]
)
GAUSS = 1.0 / math.sqrt(3.0)


def shape(natural):
    values = np.prod(1.0 + SIGNS * natural, axis=1) / 8.0
    derivatives = np.empty((8, 3))
    for node in range(8):
        for direction in range(3):
            other = [value for value in range(3) if value != direction]
            derivatives[node, direction] = (
                SIGNS[node, direction]
                * (1.0 + SIGNS[node, other[0]] * natural[other[0]])
                * (1.0 + SIGNS[node, other[1]] * natural[other[1]])
                / 8.0
            )
    return values, derivatives


def geometry(coordinates):
    volume = 0.0
    gradient_numerator = np.zeros((8, 3))
    weights = np.zeros(8)
    for xi in (-GAUSS, GAUSS):
        for eta in (-GAUSS, GAUSS):
            for zeta in (-GAUSS, GAUSS):
                values, derivatives = shape(np.array((xi, eta, zeta)))
                jacobian = coordinates.T @ derivatives
                determinant = np.linalg.det(jacobian)
                if determinant <= 0.0:
                    raise RuntimeError("nonpositive geometry")
                gradient = derivatives @ np.linalg.inv(jacobian)
                volume += determinant
                gradient_numerator += determinant * gradient
                weights += determinant * values
    average_gradient = gradient_numerator / volume
    gamma = RAW - average_gradient @ (coordinates.T @ RAW)
    return volume, average_gradient, weights, gamma


def thermal_coefficients(volume, gradient):
    inverse_effective = SIGNS.T @ gradient
    effective = np.linalg.inv(inverse_effective)
    metric = effective.T @ effective
    pivots = np.array(
        [
            metric[0, 0],
            metric[1, 1] - metric[0, 1] ** 2 / metric[0, 0],
            metric[2, 2]
            - (
                metric[1, 1] * metric[0, 2] ** 2
                - 2.0 * metric[0, 1] * metric[0, 2] * metric[1, 2]
                + metric[0, 0] * metric[1, 2] ** 2
            )
            / (metric[0, 0] * metric[1, 1] - metric[0, 1] ** 2),
        ]
    )
    inverse_length = 1.0 / pivots
    scale = volume / 192.0
    return scale * np.array(
        [
            inverse_length[0] + inverse_length[1],
            inverse_length[0] + inverse_length[2],
            inverse_length[0] + inverse_length[2],
            np.sum(inverse_length) / 3.0,
        ]
    )


def rotation_and_strain(displacement, midpoint_gradient):
    rate = displacement.T @ midpoint_gradient
    strain = 0.5 * (rate + rate.T)
    spin = 0.5 * (rate - rate.T)
    rotation = (np.eye(3) + 0.5 * spin) @ np.linalg.inv(np.eye(3) - 0.5 * spin)
    return rotation, rotation.T @ strain @ rotation


def candidate(reference, state, hourglass_rule, thermal_geometry, mechanics_geometry):
    temperature = state[:8]
    displacement = state[8:].reshape(3, 8).T
    current = reference + displacement
    midpoint = reference + 0.5 * displacement
    reference_volume, reference_gradient, reference_weights, reference_gamma = geometry(reference)
    current_volume, current_gradient, current_weights, current_gamma = geometry(current)
    _, midpoint_gradient, midpoint_weights, midpoint_gamma = geometry(midpoint)
    if mechanics_geometry == "uniform":
        mechanical_volume = current_volume
        mechanical_gradient = current_gradient
        strain_gradient = midpoint_gradient
    elif mechanics_geometry == "effective":
        deformation = np.eye(3) + displacement.T @ reference_gradient
        midpoint_deformation = np.eye(3) + 0.5 * displacement.T @ reference_gradient
        mechanical_volume = reference_volume * np.linalg.det(deformation)
        mechanical_gradient = reference_gradient @ np.linalg.inv(deformation)
        strain_gradient = reference_gradient @ np.linalg.inv(midpoint_deformation)
    elif mechanics_geometry == "mixed":
        mechanical_volume = current_volume
        mechanical_gradient = current_gradient
        midpoint_deformation = np.eye(3) + 0.5 * displacement.T @ reference_gradient
        strain_gradient = reference_gradient @ np.linalg.inv(midpoint_deformation)
    else:
        raise RuntimeError(mechanics_geometry)
    rotation, strain = rotation_and_strain(displacement, strain_gradient)
    material_temperature = np.dot(current_weights, temperature) / current_volume
    shear = 2.0e11 / 2.5
    lame = 2.0e11 * 0.25 / (1.25 * 0.5)
    elastic = strain - np.eye(3) * 1.2e-5 * (material_temperature - 300.0)
    stress_material = 2.0 * shear * elastic + lame * np.trace(elastic) * np.eye(3)
    stress = rotation @ stress_material @ rotation.T
    mechanical = mechanical_volume * mechanical_gradient @ stress.T

    if hourglass_rule == "none":
        gamma_force = reference_gamma
        gamma_increment = reference_gamma
        metric_gradient = reference_gradient
        metric_volume = 0.0
    elif hourglass_rule == "reference":
        gamma_force = reference_gamma
        gamma_increment = reference_gamma
        metric_gradient = reference_gradient
        metric_volume = reference_volume
    elif hourglass_rule == "midpoint":
        gamma_force = current_gamma
        gamma_increment = midpoint_gamma
        metric_gradient = reference_gradient
        metric_volume = reference_volume
    elif hourglass_rule == "current_metric":
        gamma_force = current_gamma
        gamma_increment = midpoint_gamma
        metric_gradient = current_gradient
        metric_volume = current_volume
    elif hourglass_rule == "reference_rotated":
        gamma_force = reference_gamma
        gamma_increment = reference_gamma
        metric_gradient = reference_gradient
        metric_volume = reference_volume
    elif hourglass_rule == "current_rotated":
        gamma_force = current_gamma
        gamma_increment = midpoint_gamma
        metric_gradient = reference_gradient
        metric_volume = reference_volume
    else:
        raise RuntimeError(hourglass_rule)
    modal = gamma_increment.T @ displacement
    if hourglass_rule in ("reference_rotated", "current_rotated"):
        modal = modal @ rotation.T
    metrics = math.sqrt(2.0) * metric_volume * np.sum(metric_gradient * metric_gradient, axis=0) / 6.0
    hourglass_force = np.zeros((8, 3))
    for component in range(3):
        hourglass_force[:, component] = (
            0.005 * shear * metrics[component] * gamma_force @ modal[:, component]
        )
    mechanical += hourglass_force

    if thermal_geometry == "reference":
        thermal_volume = reference_volume
        thermal_gradient = reference_gradient
        thermal_gamma = reference_gamma
    elif thermal_geometry == "current":
        thermal_volume = current_volume
        thermal_gradient = current_gradient
        thermal_gamma = current_gamma
    else:
        raise RuntimeError(thermal_geometry)
    coefficients = thermal_coefficients(thermal_volume, thermal_gradient)
    thermal = 4.0 * (
        thermal_volume * thermal_gradient @ (thermal_gradient.T @ temperature)
        + thermal_gamma @ (coefficients * (thermal_gamma.T @ temperature))
    )
    result = np.zeros(32)
    result[:8] = thermal
    result[8:] = mechanical.T.reshape(24)
    return result


def read(path):
    cases = {}
    with open(path, newline="") as source:
        for row in csv.DictReader(source):
            case = cases.setdefault(row["case"], (np.zeros(32), np.zeros(32)))
            node = int(row["node"]) - 1
            case[0][node] = float(row["temperature_k"])
            case[0][8 + node] = float(row["ux_m"])
            case[0][16 + node] = float(row["uy_m"])
            case[0][24 + node] = float(row["uz_m"])
            case[1][node] = float(row["rfl_w"])
            case[1][8 + node] = float(row["rf_x_n"])
            case[1][16 + node] = float(row["rf_y_n"])
            case[1][24 + node] = float(row["rf_z_n"])
    return cases


def relative(actual, expected, rows):
    return np.linalg.norm(actual[rows] - expected[rows]) / np.linalg.norm(expected[rows])


def main(path, geometry_name):
    coordinates = {
        "regular": np.array(
            [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
            dtype=float,
        ),
        "warped": np.array(
            [
                (0.00, 0.00, 0.00),
                (1.20, 0.10, -0.05),
                (1.10, 1.00, 0.10),
                (-0.10, 0.90, 0.00),
                (0.05, -0.05, 1.00),
                (1.15, 0.00, 1.20),
                (1.00, 1.10, 1.10),
                (-0.05, 1.00, 0.90),
            ]
        ),
    }[geometry_name]
    cases = read(path)
    expected = np.zeros((32, 32))
    for column in range(32):
        perturbation = 1.0e-3 if column < 8 else 1.0e-7
        expected[:, column] = (
            cases["D%02d_PLUS" % column][1] - cases["D%02d_MINUS" % column][1]
        ) / (2.0 * perturbation)
    for mechanics in ("uniform", "effective", "mixed"):
      for hourglass in ("none", "reference", "midpoint", "current_metric"):
        for thermal in ("current",):
            actual = np.zeros((32, 32))
            for column in range(32):
                perturbation = 1.0e-5 if column < 8 else 1.0e-7
                plus = cases["BASE"][0].copy()
                minus = cases["BASE"][0].copy()
                plus[column] += perturbation
                minus[column] -= perturbation
                actual[:, column] = (
                    candidate(coordinates, plus, hourglass, thermal, mechanics)
                    - candidate(coordinates, minus, hourglass, thermal, mechanics)
                ) / (2.0 * perturbation)
            base = candidate(coordinates, cases["BASE"][0], hourglass, thermal, mechanics)
            print(
                "%s/%s/%s base_t=%.6e base_u=%.6e Ktt=%.6e Ktu=%.6e Kut=%.6e Kuu=%.6e"
                % (
                    mechanics,
                    hourglass,
                    thermal,
                    relative(base, cases["BASE"][1], slice(0, 8)),
                    relative(base, cases["BASE"][1], slice(8, 32)),
                    relative(actual[:8, :8], expected[:8, :8], np.s_[:]),
                    relative(actual[:8, 8:], expected[:8, 8:], np.s_[:]),
                    relative(actual[8:, :8], expected[8:, :8], np.s_[:]),
                    relative(actual[8:, 8:], expected[8:, 8:], np.s_[:]),
                )
            )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
