#!/usr/bin/env python3
"""Separate B5.44 stress, assembly, and hourglass contributions to reactions."""

import csv
import sys

import numpy as np

from analyze_b532_candidates import geometry, shape
from analyze_b533_candidates import finite_total_stiffness


X_NODES, Y_NODES = 5, 2


def node_index(x, y, z):
    return z * X_NODES * Y_NODES + y * X_NODES + x


def reference_coordinate(x, y, z):
    x_coordinate = float(x)
    y_coordinate = float(y)
    z_coordinate = 0.5 * float(z)
    if x == 0 or x + 1 == X_NODES:
        return (x_coordinate, y_coordinate, z_coordinate)
    alternating = -1.0 if x % 2 == 0 else 1.0
    return (
        x_coordinate + 0.06 * (2.0 * y - 1.0) * (z - 1.0),
        y_coordinate + 0.04 * alternating * (z - 1.0),
        z_coordinate + 0.05 * alternating * (y - 0.5),
    )


COORDINATES = np.asarray(
    [reference_coordinate(x, y, z) for z in range(3) for y in range(Y_NODES) for x in range(X_NODES)]
)
ELEMENTS = []
for z in range(2):
    for x in range(4):
        ELEMENTS.append(
            [
                node_index(x, 0, z),
                node_index(x + 1, 0, z),
                node_index(x + 1, 1, z),
                node_index(x, 1, z),
                node_index(x, 0, z + 1),
                node_index(x + 1, 0, z + 1),
                node_index(x + 1, 1, z + 1),
                node_index(x, 1, z + 1),
            ]
        )


def relative(actual, expected):
    return np.linalg.norm(actual - expected) / np.linalg.norm(expected)


def grouped_metrics(actual, expected):
    differences = np.linalg.norm(actual - expected, axis=2)
    references = np.linalg.norm(expected, axis=2)
    return (
        np.linalg.norm(actual - expected) / np.linalg.norm(expected),
        np.max(differences) / np.max(references),
        np.max(differences[references != 0.0] / references[references != 0.0]),
    )


def center_geometry(coordinates):
    _, derivatives = shape(np.zeros(3))
    jacobian = coordinates.T @ derivatives
    determinant = np.linalg.det(jacobian)
    return 8.0 * determinant, derivatives @ np.linalg.inv(jacobian)


def main(nodal_path, integration_path):
    states = {}
    reactions = {}
    with open(nodal_path, newline="") as source:
        for row in csv.DictReader(source):
            increment = int(row["increment"])
            node = int(row["node"]) - 1
            states.setdefault(increment, np.zeros((len(COORDINATES), 3)))[node] = [
                float(row["u1_m"]),
                float(row["u2_m"]),
                float(row["u3_m"]),
            ]
            reactions.setdefault(increment, np.zeros((len(COORDINATES), 3)))[node] = [
                float(row["rf1_n"]),
                float(row["rf2_n"]),
                float(row["rf3_n"]),
            ]
    stresses = {}
    volumes = {}
    with open(integration_path, newline="") as source:
        for row in csv.DictReader(source):
            increment = int(row["increment"])
            element = int(row["element"]) - 1
            stresses.setdefault(increment, np.zeros((len(ELEMENTS), 3, 3)))[element] = [
                [float(row["s11_pa"]), float(row["s12_pa"]), float(row["s13_pa"])],
                [float(row["s12_pa"]), float(row["s22_pa"]), float(row["s23_pa"])],
                [float(row["s13_pa"]), float(row["s23_pa"]), float(row["s33_pa"])],
            ]
            volumes.setdefault(increment, np.zeros(len(ELEMENTS)))[element] = float(row["ivol_m3"])

    constrained = [node_index(0, y, z) for z in range(3) for y in range(Y_NODES)]
    all_predicted = []
    all_expected = []
    all_uniform = []
    all_center = []
    all_mixed = []
    all_effective = []
    for increment in sorted(states):
        assembled = np.zeros_like(states[increment])
        assembled_uniform = np.zeros_like(states[increment])
        assembled_center = np.zeros_like(states[increment])
        assembled_mixed = np.zeros_like(states[increment])
        assembled_effective = np.zeros_like(states[increment])
        for element, connectivity in enumerate(ELEMENTS):
            connectivity = np.asarray(connectivity)
            reference = COORDINATES[connectivity]
            displacement = states[increment][connectivity]
            current = reference + displacement
            current_volume, current_gradient, _, _ = geometry(current)
            if abs(current_volume - volumes[increment][element]) > 2.0e-6:
                raise RuntimeError("current volume does not match the Abaqus integration output")
            uniform = current_volume * current_gradient @ stresses[increment][element].T
            center_volume, center_gradient = center_geometry(current)
            center = center_volume * center_gradient @ stresses[increment][element].T
            reference_volume, reference_gradient, _, _ = geometry(reference)
            deformation = np.eye(3) + displacement.T @ reference_gradient
            pushed_gradient = reference_gradient @ np.linalg.inv(deformation)
            mixed = current_volume * pushed_gradient @ stresses[increment][element].T
            effective = reference_volume * np.linalg.det(deformation) * pushed_gradient @ stresses[increment][element].T
            hourglass = finite_total_stiffness(reference, displacement) * (8.0e7 / (2.0e11 / 2.5))
            assembled_uniform[connectivity] += uniform
            assembled_center[connectivity] += center + hourglass
            assembled_mixed[connectivity] += mixed + hourglass
            assembled_effective[connectivity] += effective + hourglass
            assembled[connectivity] += uniform + hourglass
        predicted = assembled[constrained]
        uniform = assembled_uniform[constrained]
        center = assembled_center[constrained]
        mixed = assembled_mixed[constrained]
        effective = assembled_effective[constrained]
        expected = reactions[increment][constrained]
        all_predicted.append(predicted)
        all_uniform.append(uniform)
        all_center.append(center)
        all_mixed.append(mixed)
        all_effective.append(effective)
        all_expected.append(expected)
        print(
            "increment=%d reconstructed=%.9e center=%.9e mixed=%.9e effective=%.9e without_hourglass=%.9e "
            "hourglass_relative=%.9e hourglass_norm=%.9e"
            % (increment, relative(predicted, expected), relative(center, expected), relative(mixed, expected),
               relative(effective, expected), relative(uniform, expected), relative(predicted - uniform, expected - uniform),
               np.linalg.norm(predicted - uniform))
        )
    reconstructed_metrics = grouped_metrics(np.asarray(all_predicted), np.asarray(all_expected))
    print(
        "all_reconstructed_l2=%.9e all_reconstructed_absolute_peak=%.9e "
        "all_reconstructed_pointwise=%.9e" % reconstructed_metrics
    )
    print("all_center=%.9e" % relative(np.asarray(all_center), np.asarray(all_expected)))
    print("all_mixed=%.9e" % relative(np.asarray(all_mixed), np.asarray(all_expected)))
    print("all_effective=%.9e" % relative(np.asarray(all_effective), np.asarray(all_expected)))
    print("all_without_hourglass=%.9e" % relative(np.asarray(all_uniform), np.asarray(all_expected)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
