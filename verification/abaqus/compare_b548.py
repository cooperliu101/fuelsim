#!/usr/bin/env python3
"""Compare final Fuelsim HEX20 results with a full-size Abaqus reference."""

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
        self.maximum_absolute_difference = 0.0
        self.maximum_absolute_reference = 0.0
        self.maximum_pointwise_relative = 0.0
        self.maximum_pointwise_actual = 0.0
        self.maximum_pointwise_reference = 0.0
        self.zero_reference_count = 0
        self.maximum_zero_reference_absolute_difference = 0.0
        self.count = 0

    def add(self, actual, reference):
        difference = abs(actual - reference)
        self.difference_squared += difference * difference
        self.reference_squared += reference * reference
        self.maximum_absolute_difference = max(self.maximum_absolute_difference, difference)
        self.maximum_absolute_reference = max(self.maximum_absolute_reference, abs(reference))
        if reference == 0.0:
            self.zero_reference_count += 1
            self.maximum_zero_reference_absolute_difference = max(
                self.maximum_zero_reference_absolute_difference, difference
            )
        else:
            relative = difference / abs(reference)
            if relative > self.maximum_pointwise_relative:
                self.maximum_pointwise_relative = relative
                self.maximum_pointwise_actual = actual
                self.maximum_pointwise_reference = reference
        self.count += 1

    def row(self, name, actual_maximum, reference_maximum):
        relative_l2 = (
            math.sqrt(self.difference_squared / self.reference_squared)
            if self.reference_squared > 0.0
            else (0.0 if self.difference_squared == 0.0 else math.nan)
        )
        relative_peak = (
            self.maximum_absolute_difference / self.maximum_absolute_reference
            if self.maximum_absolute_reference > 0.0
            else (0.0 if self.maximum_absolute_difference == 0.0 else math.nan)
        )
        return {
            "field": name,
            "count": self.count,
            "relative_l2_percent": 100.0 * relative_l2,
            "relative_absolute_peak_percent": 100.0 * relative_peak,
            "maximum_pointwise_relative_percent": 100.0 * self.maximum_pointwise_relative,
            "zero_reference_count": self.zero_reference_count,
            "maximum_zero_reference_absolute_difference": self.maximum_zero_reference_absolute_difference,
            "maximum_absolute_difference": self.maximum_absolute_difference,
            "maximum_pointwise_actual": self.maximum_pointwise_actual,
            "maximum_pointwise_reference": self.maximum_pointwise_reference,
            "fuelsim_maximum": actual_maximum,
            "abaqus_maximum": reference_maximum,
        }


def decode_names(variable):
    result = []
    for row in variable[:]:
        result.append(b"".join(np.ma.getdata(row).tolist()).decode("ascii").strip("\x00 "))
    return result


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


def quad8_shapes(xi, eta):
    shape, derivative_xi, derivative_eta = np.zeros(8), np.zeros(8), np.zeros(8)
    for node, (sx, sy) in enumerate(((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))):
        ax, ay, total = 1.0 + sx * xi, 1.0 + sy * eta, sx * xi + sy * eta - 1.0
        shape[node] = 0.25 * ax * ay * total
        derivative_xi[node] = 0.25 * sx * ay * (total + ax)
        derivative_eta[node] = 0.25 * sy * ax * (total + ay)
    shape[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta)
    derivative_xi[4], derivative_eta[4] = -xi * (1.0 - eta), -0.5 * (1.0 - xi * xi)
    shape[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta)
    derivative_xi[5], derivative_eta[5] = 0.5 * (1.0 - eta * eta), -(1.0 + xi) * eta
    shape[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta)
    derivative_xi[6], derivative_eta[6] = -xi * (1.0 + eta), 0.5 * (1.0 - xi * xi)
    shape[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta)
    derivative_xi[7], derivative_eta[7] = -0.5 * (1.0 - eta * eta), -(1.0 - xi) * eta
    return shape, derivative_xi, derivative_eta


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
    parser.add_argument("--case-name", default="B5.48")
    parser.add_argument("fuelsim_exodus")
    parser.add_argument("abaqus_nodal")
    parser.add_argument("abaqus_contact")
    parser.add_argument("abaqus_clad_points")
    parser.add_argument("abaqus_mesh")
    parser.add_argument("output")
    arguments = parser.parse_args()
    output_prefix = arguments.case_name.lower().replace(".", "")

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
        fuel_connectivity = np.asarray(database.variables["connect1"][:], dtype=int) - 1
        connectivity = np.asarray(database.variables["connect2"][:], dtype=int) - 1
        clad_variables = {
            name: np.asarray(database.variables["vals_elem_var%deb2" % (index + 1)][-1], dtype=float)
            for index, name in enumerate(element_names)
        }
        global_names = decode_names(database.variables["name_glo_var"])
        global_values = dict(zip(global_names, np.asarray(database.variables["vals_glo_var"][-1], dtype=float)))
    finally:
        database.close()

    displacement = np.column_stack(
        [nodal["displacement_x"], nodal["displacement_y"], nodal["displacement_z"]]
    )
    current_coordinates = coordinates + displacement
    abaqus_nodes = dict((int(row["id"]), row) for row in read_csv(arguments.abaqus_nodal))
    if len(abaqus_nodes) != len(coordinates):
        raise RuntimeError("%s nodal result count differs" % arguments.case_name)
    corner_nodes = set((np.vstack((fuel_connectivity[:, :8], connectivity[:, :8])) + 1).reshape(-1).tolist())

    metrics = {}
    maxima = {}
    for name in ("temperature_corner", "temperature_interpolated", "displacement_x", "displacement_y", "displacement_z", "radial_displacement", "tangential_displacement"):
        metrics[name] = Metrics()
        maxima[name] = [0.0, 0.0]
    axis_transverse_maximum = 0.0
    for label, reference_row in abaqus_nodes.items():
        index = label - 1
        actual_values = (
            nodal["temperature"][index],
            displacement[index, 0],
            displacement[index, 1],
            displacement[index, 2],
        )
        reference_values = tuple(float(reference_row[name]) for name in ("T", "disp_x", "disp_y", "disp_z"))
        temperature_name = "temperature_corner" if label in corner_nodes else "temperature_interpolated"
        metrics[temperature_name].add(actual_values[0], reference_values[0])
        maxima[temperature_name][0] = max(maxima[temperature_name][0], abs(actual_values[0]))
        maxima[temperature_name][1] = max(maxima[temperature_name][1], abs(reference_values[0]))
        for component, name in enumerate(("displacement_x", "displacement_y", "displacement_z"), 1):
            metrics[name].add(actual_values[component], reference_values[component])
            maxima[name][0] = max(maxima[name][0], abs(actual_values[component]))
            maxima[name][1] = max(maxima[name][1], abs(reference_values[component]))
        x, y = coordinates[index, :2]
        radius = math.hypot(x, y)
        if radius == 0.0:
            axis_transverse_maximum = max(axis_transverse_maximum, math.hypot(actual_values[1], actual_values[2]))
        else:
            actual_radial = (x * actual_values[1] + y * actual_values[2]) / radius
            reference_radial = (x * reference_values[1] + y * reference_values[2]) / radius
            actual_tangential = (-y * actual_values[1] + x * actual_values[2]) / radius
            reference_tangential = (-y * reference_values[1] + x * reference_values[2]) / radius
            for name, actual, reference in (
                ("radial_displacement", actual_radial, reference_radial),
                ("tangential_displacement", actual_tangential, reference_tangential),
            ):
                metrics[name].add(actual, reference)
                maxima[name][0] = max(maxima[name][0], abs(actual))
                maxima[name][1] = max(maxima[name][1], abs(reference))

    contact = read_csv(arguments.abaqus_contact)
    if len(contact) != 416:
        raise RuntimeError("%s contact result count differs" % arguments.case_name)
    contact_metrics = Metrics()
    contact_maxima = [0.0, 0.0]
    for row in contact:
        index = int(row["id"]) - 1
        actual = nodal["contact_pressure_fuel_cladding"][index]
        reference = float(row["contact_pressure"])
        if not math.isfinite(actual):
            raise RuntimeError("%s Fuelsim contact pressure is missing" % arguments.case_name)
        contact_metrics.add(actual, reference)
        contact_maxima[0] = max(contact_maxima[0], abs(actual))
        contact_maxima[1] = max(contact_maxima[1], abs(reference))
    metrics["contact_pressure"] = contact_metrics
    maxima["contact_pressure"] = contact_maxima

    with Path(arguments.abaqus_mesh).open() as source:
        mesh = json.load(source)
    element_records = dict((element["label"], element) for element in mesh["elements"])
    face_nodes = (
        (0, 1, 5, 4, 8, 13, 16, 12),
        (1, 2, 6, 5, 9, 14, 17, 13),
        (2, 3, 7, 6, 10, 15, 18, 14),
        (3, 0, 4, 7, 11, 12, 19, 15),
        (0, 3, 2, 1, 11, 10, 9, 8),
        (4, 5, 6, 7, 16, 17, 18, 19),
    )
    fuel_outer_faces = next(side_set["faces"] for side_set in mesh["side_sets"] if side_set["name"] == "FUEL_OUTER")
    contact_by_label = dict((int(row["id"]), row) for row in contact)
    gauss_surface = (-math.sqrt(3.0 / 5.0), 0.0, math.sqrt(3.0 / 5.0))
    weights_surface = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)
    abaqus_contact_force = 0.0
    abaqus_tangential_force = 0.0
    abaqus_heat_rate = 0.0
    abaqus_contact_area = 0.0
    for face in fuel_outer_faces:
        element_nodes = element_records[face["element"]]["nodes"]
        labels = [element_nodes[index] for index in face_nodes[face["exodus_side"] - 1]]
        rows = [contact_by_label[label] for label in labels]
        points = np.asarray([[float(row[name]) for name in ("current_x", "current_y", "current_z")] for row in rows])
        pressure = np.asarray([float(row["contact_pressure"]) for row in rows])
        shear = np.asarray([math.hypot(float(row["shear_1"]), float(row["shear_2"])) for row in rows])
        heat = np.asarray([float(row["heat_flow"]) for row in rows])
        for xi_index, xi in enumerate(gauss_surface):
            for eta_index, eta in enumerate(gauss_surface):
                shape, derivative_xi, derivative_eta = quad8_shapes(xi, eta)
                area = (
                    np.linalg.norm(np.cross(derivative_xi.dot(points), derivative_eta.dot(points)))
                    * weights_surface[xi_index]
                    * weights_surface[eta_index]
                )
                abaqus_contact_area += area
                abaqus_contact_force += shape.dot(pressure) * area
                abaqus_tangential_force += shape.dot(shear) * area
                abaqus_heat_rate += shape.dot(heat) * area
    for name, actual, reference in (
        ("recovered_contact_pressure_integral", global_values["contact_force_fuel_cladding"], abaqus_contact_force),
        ("recovered_contact_shear_integral", global_values["contact_tangential_force_fuel_cladding"], abaqus_tangential_force),
        ("recovered_contact_heat_integral", global_values["contact_heat_rate_fuel_cladding"], abs(abaqus_heat_rate)),
    ):
        scalar = Metrics()
        scalar.add(actual, reference)
        metrics[name] = scalar
        maxima[name] = [abs(actual), abs(reference)]

    point_rows = read_csv(arguments.abaqus_clad_points)
    if len(point_rows) != 512 * 27:
        raise RuntimeError("%s clad integration-point result count differs" % arguments.case_name)
    reference_by_element = {}
    for row in point_rows:
        reference_by_element.setdefault(int(row["element"]), []).append(row)
    for name in ("equivalent_stress", "equivalent_plastic_strain", "equivalent_creep_strain"):
        metrics[name] = Metrics()
        maxima[name] = [0.0, 0.0]
    gauss = math.sqrt(3.0 / 5.0)
    natural_points = [
        (xi, eta, zeta)
        for zeta in (-gauss, 0.0, gauss)
        for eta in (-gauss, 0.0, gauss)
        for xi in (-gauss, 0.0, gauss)
    ]
    shape_values = [hex20_shapes(*point) for point in natural_points]
    maximum_point_coordinate_difference = 0.0
    minimum_second_to_first_distance_ratio = math.inf
    for local_element, element_nodes in enumerate(connectivity):
        label = local_element + 641
        reference_rows = reference_by_element[label]
        reference_points = np.asarray(
            [[float(row[name]) for name in ("current_x", "current_y", "current_z")] for row in reference_rows]
        )
        actual_points = np.asarray([shape.dot(current_coordinates[element_nodes]) for shape in shape_values])
        distances = np.linalg.norm(actual_points[:, None, :] - reference_points[None, :, :], axis=2)
        pairs = []
        unused_actual, unused_reference = set(range(27)), set(range(27))
        for distance, actual_point, reference_point in sorted(
            (distances[actual, reference], actual, reference) for actual in range(27) for reference in range(27)
        ):
            if actual_point in unused_actual and reference_point in unused_reference:
                pairs.append((actual_point, reference_point, distance))
                unused_actual.remove(actual_point)
                unused_reference.remove(reference_point)
        if unused_actual or unused_reference:
            raise RuntimeError("%s integration-point coordinate matching failed" % arguments.case_name)
        for actual_point, reference_point, distance in pairs:
            maximum_point_coordinate_difference = max(maximum_point_coordinate_difference, distance)
            sorted_distances = np.sort(distances[actual_point])
            if sorted_distances[0] > 0.0:
                minimum_second_to_first_distance_ratio = min(
                    minimum_second_to_first_distance_ratio, sorted_distances[1] / sorted_distances[0]
                )
            stress = [clad_variables["stress_%s_q%d" % (component, actual_point)][local_element] for component in ("xx", "yy", "zz", "xy", "yz", "xz")]
            actual_fields = (
                equivalent_stress(stress),
                clad_variables["equiv_plastic_q%d" % actual_point][local_element],
                clad_variables["equiv_creep_q%d" % actual_point][local_element],
            )
            row = reference_rows[reference_point]
            reference_fields = tuple(
                float(row[name])
                for name in ("vonmises_stress", "effective_plastic_strain", "effective_creep_strain")
            )
            for name, actual, reference in zip(
                ("equivalent_stress", "equivalent_plastic_strain", "equivalent_creep_strain"),
                actual_fields,
                reference_fields,
            ):
                metrics[name].add(actual, reference)
                maxima[name][0] = max(maxima[name][0], abs(actual))
                maxima[name][1] = max(maxima[name][1], abs(reference))

    rows = [metrics[name].row(name, *maxima[name]) for name in metrics]
    fieldnames = list(rows[0].keys())
    with Path(arguments.output).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("%s_axis_transverse_analytic_zero_maximum_absolute=%.12e" % (output_prefix, axis_transverse_maximum))
    print("%s_abaqus_integrated_contact_area=%.12e" % (output_prefix, abaqus_contact_area))
    print("%s_integration_point_maximum_coordinate_difference=%.12e" %
          (output_prefix, maximum_point_coordinate_difference))
    print("%s_integration_point_minimum_second_to_first_distance_ratio=%.12e" %
          (output_prefix, minimum_second_to_first_distance_ratio))
    for row in rows:
        print(
            "%s_%s relative_l2=%.6g%% relative_peak=%.6g%% pointwise=%.6g%% max_abs=%.12e zeros=%d"
            % (
                output_prefix,
                row["field"],
                row["relative_l2_percent"],
                row["relative_absolute_peak_percent"],
                row["maximum_pointwise_relative_percent"],
                row["maximum_absolute_difference"],
                row["zero_reference_count"],
            )
        )


if __name__ == "__main__":
    main()
