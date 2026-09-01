#!/usr/bin/env python3
"""Compare final Fuelsim HEX20 fields with the tracked Abaqus B5.49 C3D20T result."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


class Metrics:
    def __init__(self):
        self.difference_squared = 0.0
        self.reference_squared = 0.0
        self.maximum_difference = 0.0
        self.maximum_reference = 0.0
        self.maximum_actual = 0.0
        self.maximum_pointwise_relative = 0.0
        self.maximum_pointwise_actual = 0.0
        self.maximum_pointwise_reference = 0.0
        self.maximum_absolute_difference = 0.0
        self.zero_reference_count = 0
        self.maximum_zero_reference_difference = 0.0
        self.count = 0

    def add(self, actual, reference):
        actual_values = np.atleast_1d(np.asarray(actual, dtype=float))
        reference_values = np.atleast_1d(np.asarray(reference, dtype=float))
        difference_values = actual_values - reference_values
        difference = float(np.linalg.norm(difference_values))
        reference_norm = float(np.linalg.norm(reference_values))
        self.difference_squared += difference * difference
        self.reference_squared += reference_norm * reference_norm
        self.maximum_difference = max(self.maximum_difference, difference)
        self.maximum_reference = max(self.maximum_reference, reference_norm)
        self.maximum_actual = max(self.maximum_actual, float(np.linalg.norm(actual_values)))
        self.maximum_absolute_difference = max(
            self.maximum_absolute_difference, float(np.max(np.abs(difference_values)))
        )
        if reference_norm == 0.0:
            self.zero_reference_count += 1
            self.maximum_zero_reference_difference = max(self.maximum_zero_reference_difference, difference)
        else:
            relative = difference / reference_norm
            if relative > self.maximum_pointwise_relative:
                self.maximum_pointwise_relative = relative
                self.maximum_pointwise_actual = float(np.linalg.norm(actual_values))
                self.maximum_pointwise_reference = reference_norm
        self.count += 1

    def row(self, name):
        return {
            "field": name,
            "count": self.count,
            "relative_l2_percent": 100.0 * math.sqrt(self.difference_squared / self.reference_squared),
            "relative_absolute_peak_percent": 100.0 * abs(self.maximum_actual - self.maximum_reference) / self.maximum_reference,
            "maximum_pointwise_relative_percent": 100.0 * self.maximum_pointwise_relative,
            "zero_reference_count": self.zero_reference_count,
            "maximum_zero_reference_absolute_difference": self.maximum_zero_reference_difference,
            "maximum_absolute_component_difference": self.maximum_absolute_difference,
            "maximum_pointwise_actual_norm": self.maximum_pointwise_actual,
            "maximum_pointwise_reference_norm": self.maximum_pointwise_reference,
        }


def decode_names(variable):
    return [b"".join(np.ma.getdata(row).tolist()).decode("ascii").strip("\x00 ") for row in variable[:]]


def read_csv(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def hex20_shapes(xi, eta, zeta):
    shape = np.zeros(20)
    signs = (
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    )
    for node, (sx, sy, sz) in enumerate(signs):
        ax, ay, az = 1.0 + sx * xi, 1.0 + sy * eta, 1.0 + sz * zeta
        shape[node] = 0.125 * ax * ay * az * (sx * xi + sy * eta + sz * zeta - 2.0)
    shape[8] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 - zeta)
    shape[9] = 0.25 * (1.0 - eta * eta) * (1.0 + xi) * (1.0 - zeta)
    shape[10] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 - zeta)
    shape[11] = 0.25 * (1.0 - eta * eta) * (1.0 - xi) * (1.0 - zeta)
    shape[12] = 0.25 * (1.0 - zeta * zeta) * (1.0 - xi) * (1.0 - eta)
    shape[13] = 0.25 * (1.0 - zeta * zeta) * (1.0 + xi) * (1.0 - eta)
    shape[14] = 0.25 * (1.0 - zeta * zeta) * (1.0 + xi) * (1.0 + eta)
    shape[15] = 0.25 * (1.0 - zeta * zeta) * (1.0 - xi) * (1.0 + eta)
    shape[16] = 0.25 * (1.0 - xi * xi) * (1.0 - eta) * (1.0 + zeta)
    shape[17] = 0.25 * (1.0 - eta * eta) * (1.0 + xi) * (1.0 + zeta)
    shape[18] = 0.25 * (1.0 - xi * xi) * (1.0 + eta) * (1.0 + zeta)
    shape[19] = 0.25 * (1.0 - eta * eta) * (1.0 - xi) * (1.0 + zeta)
    return shape


def equivalent_stress(components):
    mean = sum(components[:3]) / 3.0
    return math.sqrt(
        1.5
        * (
            sum((component - mean) ** 2 for component in components[:3])
            + 2.0 * sum(component * component for component in components[3:])
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fuelsim_exodus")
    parser.add_argument("abaqus_temperature")
    parser.add_argument("abaqus_displacement")
    parser.add_argument("abaqus_material")
    parser.add_argument("abaqus_mesh")
    parser.add_argument("output")
    arguments = parser.parse_args()

    database = Dataset(arguments.fuelsim_exodus)
    try:
        nodal_names = decode_names(database.variables["name_nod_var"])
        element_names = decode_names(database.variables["name_elem_var"])
        nodal = {
            name: np.asarray(database.variables["vals_nod_var%d" % (index + 1)][-1], dtype=float)
            for index, name in enumerate(nodal_names)
        }
        coordinates = np.column_stack(
            [np.asarray(database.variables[name][:], dtype=float) for name in ("coordx", "coordy", "coordz")]
        )
        connectivity = [
            np.asarray(database.variables["connect%d" % block][:], dtype=int) - 1 for block in (1, 2)
        ]
        element_variables = []
        for block in (1, 2):
            element_variables.append(
                {
                    name: np.asarray(
                        database.variables["vals_elem_var%deb%d" % (index + 1, block)][-1], dtype=float
                    )
                    for index, name in enumerate(element_names)
                }
            )
    finally:
        database.close()

    metrics = {
        name: Metrics()
        for name in (
            "temperature",
            "displacement_vector",
            "equivalent_stress",
            "equivalent_plastic_strain",
            "equivalent_creep_strain",
        )
    }
    maximum_coordinate_difference = 0.0
    for row in read_csv(arguments.abaqus_temperature):
        label = int(row["id"])
        maximum_coordinate_difference = max(
            maximum_coordinate_difference,
            float(
                np.max(
                    np.abs(
                        coordinates[label - 1]
                        - np.asarray([float(row[name]) for name in ("x", "y", "z")])
                    )
                )
            ),
        )
        metrics["temperature"].add(nodal["temperature"][label - 1], float(row["temperature"]))
    displacement = np.column_stack(
        [nodal["displacement_x"], nodal["displacement_y"], nodal["displacement_z"]]
    )
    for row in read_csv(arguments.abaqus_displacement):
        label = int(row["id"])
        reference = [float(row[name]) for name in ("displacement_x", "displacement_y", "displacement_z")]
        metrics["displacement_vector"].add(displacement[label - 1], reference)

    gauss = math.sqrt(3.0 / 5.0)
    shape_values = [
        hex20_shapes(xi, eta, zeta)
        for zeta in (-gauss, 0.0, gauss)
        for eta in (-gauss, 0.0, gauss)
        for xi in (-gauss, 0.0, gauss)
    ]
    current_coordinates = coordinates + displacement
    reference_by_element = {}
    for row in read_csv(arguments.abaqus_material):
        reference_by_element.setdefault(int(row["element"]), []).append(row)
    maximum_material_point_coordinate_difference = 0.0
    minimum_second_to_first_distance_ratio = math.inf
    for block in range(2):
        for local_element, element_nodes in enumerate(connectivity[block]):
            label = 1 + block * len(connectivity[0]) + local_element
            reference_rows = reference_by_element[label]
            reference_points = np.asarray(
                [[float(row[name]) for name in ("current_x", "current_y", "current_z")] for row in reference_rows]
            )
            actual_points = np.asarray([shape.dot(current_coordinates[element_nodes]) for shape in shape_values])
            distances = np.linalg.norm(actual_points[:, None, :] - reference_points[None, :, :], axis=2)
            pairs = []
            unused_actual, unused_reference = set(range(27)), set(range(27))
            for distance, actual_point, reference_point in sorted(
                (distances[actual, reference], actual, reference)
                for actual in range(27)
                for reference in range(27)
            ):
                if actual_point in unused_actual and reference_point in unused_reference:
                    pairs.append((actual_point, reference_point, distance))
                    unused_actual.remove(actual_point)
                    unused_reference.remove(reference_point)
            if unused_actual or unused_reference:
                raise RuntimeError("B5.49 integration-point coordinate matching failed")
            for actual_point, reference_point, distance in pairs:
                maximum_material_point_coordinate_difference = max(
                    maximum_material_point_coordinate_difference, distance
                )
                sorted_distances = np.sort(distances[actual_point])
                if sorted_distances[0] > 0.0:
                    minimum_second_to_first_distance_ratio = min(
                        minimum_second_to_first_distance_ratio, sorted_distances[1] / sorted_distances[0]
                    )
                variables = element_variables[block]
                stress = [
                    variables["stress_%s_q%d" % (component, actual_point)][local_element]
                    for component in ("xx", "yy", "zz", "xy", "yz", "xz")
                ]
                actual_fields = (
                    equivalent_stress(stress),
                    variables["equiv_plastic_q%d" % actual_point][local_element],
                    variables["equiv_creep_q%d" % actual_point][local_element],
                )
                reference_row = reference_rows[reference_point]
                reference_fields = tuple(
                    float(reference_row[name])
                    for name in (
                        "vonmises_stress",
                        "effective_plastic_strain",
                        "effective_creep_strain",
                    )
                )
                for name, actual, reference in zip(
                    ("equivalent_stress", "equivalent_plastic_strain", "equivalent_creep_strain"),
                    actual_fields,
                    reference_fields,
                ):
                    metrics[name].add(actual, reference)

    rows = [metrics[name].row(name) for name in metrics]
    with Path(arguments.output).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("b549_reference_coordinate_maximum_difference=%.12e" % maximum_coordinate_difference)
    print("b549_material_point_maximum_coordinate_difference=%.12e" % maximum_material_point_coordinate_difference)
    print("b549_material_point_minimum_second_to_first_distance_ratio=%.12e" % minimum_second_to_first_distance_ratio)
    passed = maximum_coordinate_difference < 1.0e-14 and minimum_second_to_first_distance_ratio > 10.0
    for row in rows:
        print(
            "b549_%s relative_l2=%.6g%% relative_peak=%.6g%% pointwise=%.6g%% max_abs=%.12e zeros=%d zero_max_abs=%.12e"
            % (
                row["field"],
                row["relative_l2_percent"],
                row["relative_absolute_peak_percent"],
                row["maximum_pointwise_relative_percent"],
                row["maximum_absolute_component_difference"],
                row["zero_reference_count"],
                row["maximum_zero_reference_absolute_difference"],
            )
        )
        print(
            "b549_%s pointwise_actual=%.12e pointwise_reference=%.12e"
            % (row["field"], row["maximum_pointwise_actual_norm"], row["maximum_pointwise_reference_norm"])
        )
        passed = passed and row["relative_l2_percent"] < 1.0 and row["relative_absolute_peak_percent"] < 1.0 and row["maximum_pointwise_relative_percent"] < 1.0
    if not passed:
        raise SystemExit("B5.49 comparison exceeds the one-percent acceptance boundary")


if __name__ == "__main__":
    main()
