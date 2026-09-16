#!/usr/bin/env python3
"""Compare complete production Exodus output with every native Abaqus frame."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, chartostring

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from compare_b548 import hex20_shapes


def columns(rows, names):
    return np.column_stack([rows[name] for name in names])


def read_reference(path):
    return np.genfromtxt(path, delimiter=",", names=True, encoding="ascii")


def metrics(name, actual, reference, labels, zero_tolerance, scope):
    actual, reference = np.asarray(actual), np.asarray(reference)
    if actual.shape != reference.shape or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(reference)):
        raise RuntimeError("Invalid field or coverage: " + name)
    vector = actual.ndim == 2
    norm = lambda value: np.linalg.norm(value, axis=1) if vector else np.abs(value)
    delta, ref, val = norm(actual - reference), norm(reference), norm(actual)
    zero = ref == 0.0
    nonzero = ~zero
    relative = np.zeros(len(ref))
    relative[nonzero] = delta[nonzero] / ref[nonzero]
    failed = np.flatnonzero(nonzero & (relative >= 0.001))
    failed_worst = failed[np.argmax(delta[failed])] if len(failed) else None
    worst = int(np.argmax(relative))
    absolute_worst = int(np.argmax(delta))
    maximum_ref = float(np.max(ref))
    l2 = float(np.linalg.norm(delta) / np.linalg.norm(ref)) if np.any(nonzero) else float("nan")
    peak = (float(np.max(delta)) if vector else abs(float(np.max(val)) - maximum_ref))
    peak = peak / maximum_ref if maximum_ref else float("nan")
    zero_difference = float(np.max(delta[zero], initial=0.0))
    passed = (not np.any(nonzero) or max(l2, peak, float(relative[worst])) < 0.001)
    passed = passed and zero_difference < zero_tolerance
    return {"scope": scope, "field": name, "samples": len(ref), "relative_l2_percent": 100 * l2,
            "relative_absolute_peak_percent": 100 * peak,
            "maximum_pointwise_relative_percent": 100 * float(relative[worst]) if np.any(nonzero) else float("nan"),
            "maximum_absolute_difference": float(delta[absolute_worst]),
            "worst_relative_sample": labels[worst], "worst_relative_reference_norm": float(ref[worst]),
            "worst_relative_absolute_difference": float(delta[worst]),
            "worst_absolute_sample": labels[absolute_worst],
            "zero_reference_samples": int(np.count_nonzero(zero)),
            "maximum_zero_reference_absolute_difference": zero_difference,
            "nonzero_pointwise_failures": len(failed),
            "largest_failed_reference_norm": float(np.max(ref[failed], initial=0.0)),
            "largest_failed_absolute_difference": float(np.max(delta[failed], initial=0.0)),
            "largest_difference_sample": labels[failed_worst] if failed_worst is not None else "",
            "zero_reference_absolute_failures": int(np.count_nonzero(zero & (delta >= zero_tolerance))),
            "zero_absolute_tolerance": zero_tolerance, "strict_0_1_percent_passed": int(passed)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--reference-directory", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-directory", type=Path)
    arguments = parser.parse_args()
    root = arguments.reference_directory
    output = arguments.output_directory or root
    output.mkdir(parents=True, exist_ok=True)
    nodes, contact, points = [read_reference(root / ("reference_" + kind + ".csv.gz"))
                              for kind in ("nodes", "contact", "points")]
    times = np.unique(nodes["time"])
    # ODB frame times are stored as float32 even with full-precision fields.
    if len(times) != 20 or not np.allclose(times, np.arange(1, 21) * 0.05, rtol=0, atol=5e-8):
        raise RuntimeError("Expected the complete twenty-increment Abaqus history")
    if (len(nodes), len(contact), len(points)) != (20 * 5969, 20 * 416, 20 * 1152 * 8):
        raise RuntimeError("Incomplete Abaqus reference history")
    nodes, contact, points = nodes.reshape(20, 5969), contact.reshape(20, 416), points.reshape(20, 1152, 8)
    gauss = 1 / np.sqrt(3)
    shapes = np.asarray([hex20_shapes(xi, eta, zeta) for zeta in (-gauss, gauss)
                         for eta in (-gauss, gauss) for xi in (-gauss, gauss)])
    collected = {}
    coverage = []
    previous_actual = np.zeros((512, 8, 2))
    previous_reference = np.zeros((512, 8, 2))

    def add(name, actual, reference, labels, zero_tolerance):
        collected.setdefault(name, [[], [], [], zero_tolerance])
        record = collected[name]
        record[0].append(actual)
        record[1].append(reference)
        record[2].append(labels)

    with Dataset(arguments.results) as database:
        database.set_auto_mask(False)
        actual_times = np.asarray(database["time_whole"][:])
        if len(actual_times) != 21 or not np.allclose(actual_times[1:], times, rtol=0, atol=5e-8):
            raise RuntimeError("Fuelsim and Abaqus must both complete the same twenty increments")
        nodal_names = {name: index + 1 for index, name in enumerate(chartostring(database["name_nod_var"][:]))}
        element_names = {name: index + 1 for index, name in enumerate(chartostring(database["name_elem_var"][:]))}
        global_names = {name: index for index, name in enumerate(chartostring(database["name_glo_var"][:]))}
        connectivity = np.concatenate([database["connect1"][:], database["connect2"][:]]) - 1
        coordinates = np.column_stack([database["coord" + axis][:] for axis in "xyz"])
        corners = np.unique(connectivity[:, :8])
        interpolated = np.setdiff1d(np.arange(5969), corners)
        if coordinates.shape != (5969, 3) or connectivity.shape != (1152, 20) or len(corners) != 1617:
            raise RuntimeError("Unexpected Fuelsim mesh size")
        components = ("xx", "yy", "zz", "xy", "yz", "xz")
        for frame, time in enumerate(times):
            def nodal(name):
                return np.asarray(database["vals_nod_var%d" % nodal_names[name]][frame + 1])

            def element(name):
                index = element_names[name]
                return np.concatenate([database["vals_elem_var%deb%d" % (index, block)][frame + 1]
                                       for block in (1, 2)])

            def material(prefix, suffixes):
                return np.stack([np.column_stack([element(prefix + suffix + "_q%d" % q) for suffix in suffixes])
                                 for q in range(8)], axis=1)

            reference_nodes, reference_contact, reference_points = nodes[frame], contact[frame], points[frame]
            if not np.array_equal(reference_nodes["node"], np.arange(1, 5970)):
                raise RuntimeError("Missing or duplicate reference node labels")
            if not np.all(reference_points["element"] == np.arange(1, 1153)[:, None]) or not np.all(
                    reference_points["point"] == np.arange(1, 9)[None, :]):
                raise RuntimeError("Missing or duplicate reference material-point labels")
            for record in (reference_nodes, reference_contact, reference_points):
                if not np.all(np.abs(record["time"] - time) < 1e-12):
                    raise RuntimeError("Reference frame ordering differs")
            if not np.all(element("material_point_count") == 8):
                raise RuntimeError("Expected exactly eight active material points")
            for q in range(8, 27):
                if not np.all(np.isnan(element("stress_xx_q%d" % q))):
                    raise RuntimeError("Inactive integration points must be NaN")
            node_labels = ["t=%.12g node=%d" % (time, node) for node in range(1, 5970)]
            for name, indices in (("temperature_corner", corners), ("temperature_interpolated", interpolated)):
                add(name, nodal("temperature")[indices], reference_nodes["temperature"][indices],
                    [node_labels[index] for index in indices], 1e-10)
            for name, fields, refs, tolerance in (
                ("displacement", ["displacement_" + axis for axis in "xyz"], ["u" + axis for axis in "xyz"], 1e-12),
                ("reaction", ["reaction_force_" + axis for axis in "xyz"], ["r" + axis for axis in "xyz"], 1e-5)):
                add(name, np.column_stack([nodal(field) for field in fields]), columns(reference_nodes, refs),
                    node_labels, tolerance)
            add("heat_reaction", nodal("reaction_heat_flux")[corners], reference_nodes["heat_reaction"][corners],
                [node_labels[index] for index in corners], 1e-7)
            secondary = reference_contact["node"].astype(int) - 1
            actual_secondary = np.flatnonzero(np.isfinite(nodal("contact_projected_fuel_cladding")))
            if len(np.unique(secondary)) != 416 or not np.array_equal(secondary, actual_secondary):
                raise RuntimeError("Missing or duplicate contact nodes")
            contact_labels = [node_labels[index] for index in secondary]
            for name, reference_name, tolerance in (("gap", "gap", 1e-10), ("pressure", "pressure", 1e-4)):
                add("contact_" + name, nodal("contact_" + name + "_fuel_cladding")[secondary],
                    reference_contact[reference_name], contact_labels, tolerance)
            for name, field, prefix, sign, tolerance in (
                ("normal_force", "normal_force_", "n", -1, 1e-5),
                ("tangential_force", "tangential_force_", "t", -1, 1e-5),
                ("slip", "total_slip_", "slip", 1, 1e-12)):
                actual = sign * np.column_stack([nodal("contact_" + field + axis + "_fuel_cladding")[secondary]
                                                  for axis in "xyz"])
                add("contact_" + name, actual, columns(reference_contact, [prefix + axis for axis in "xyz"]),
                    contact_labels, tolerance)
            actual_heat = database["vals_glo_var"][frame + 1, global_names["contact_heat_rate_fuel_cladding"]]
            add("contact_heat_rate", np.asarray([actual_heat]), np.asarray([abs(np.sum(reference_contact["heat_rate"]))]),
                ["t=%.12g" % time], 1e-7)
            # Determine numbering solely from Abaqus geometry. Fuelsim results
            # never choose or alter the matching material-point permutation.
            reference_current = coordinates + columns(reference_nodes, ["u" + axis for axis in "xyz"])
            expected = np.einsum("qn,enc->eqc", shapes, reference_current[connectivity])
            native = np.stack([reference_points[axis] for axis in "xyz"], axis=2)
            distances = np.linalg.norm(expected[:, :, None, :] - native[:, None, :, :], axis=3)
            permutation = np.argmin(distances, axis=2)
            if not np.all(np.sort(permutation, axis=1) == np.arange(8)):
                raise RuntimeError("Abaqus point matching is not bijective")
            match_error = float(np.max(np.min(distances, axis=2)))
            if match_error > 1e-8:
                raise RuntimeError("Abaqus material-point coordinate reconstruction differs: %.17g" % match_error)
            reference_points = np.take_along_axis(reference_points, permutation, axis=1)
            point_labels = ["t=%.12g element=%d point=%d" % (time, element_index + 1, reference_points["point"][element_index, q])
                            for element_index in range(1152) for q in range(8)]
            for name, prefix, ref_prefix, tolerance in (("stress", "stress_", "s_", 1e-4),
                                                       ("elastic_strain", "elastic_", "e_", 1e-12),
                                                       ("plastic_strain", "plastic_", "p_", 1e-12),
                                                       ("creep_strain", "creep_", "c_", 1e-12)):
                actual = material(prefix, components)
                reference = np.stack([reference_points[ref_prefix + component] for component in components], axis=2)
                for region, subset in (("fuel", slice(0, 640)), ("cladding", slice(640, 1152))):
                    labels = np.asarray(point_labels).reshape(1152, 8)[subset].ravel().tolist()
                    add(region + "_" + name, actual[subset].reshape(-1, 6), reference[subset].reshape(-1, 6), labels, tolerance)
                    if name == "stress":
                        def mises(stress):
                            deviator = stress.copy()
                            deviator[:, :, :3] -= np.mean(stress[:, :, :3], axis=2)[:, :, None]
                            return np.sqrt(1.5 * (np.sum(deviator[:, :, :3] ** 2, axis=2)
                                                  + 2 * np.sum(deviator[:, :, 3:] ** 2, axis=2)))
                        add(region + "_equivalent_stress", mises(actual[subset]).ravel(),
                            mises(reference[subset]).ravel(), labels, tolerance)
            current_actual, current_reference = [], []
            for name in ("equiv_plastic", "equiv_creep"):
                actual = np.column_stack([element(name + "_q%d" % q) for q in range(8)])[640:]
                add("cladding_" + name, actual.ravel(), reference_points[name][640:].ravel(), point_labels[640 * 8:], 1e-12)
                current_actual.append(actual)
                current_reference.append(reference_points[name][640:])
            current_actual, current_reference = np.stack(current_actual, axis=2), np.stack(current_reference, axis=2)
            actual_increments = current_actual - previous_actual
            reference_increments = current_reference - previous_reference
            previous_actual, previous_reference = current_actual, current_reference
            coverage.append({"time": time, "temperature_dofs": len(corners), "displacement_nodes": 5969,
                             "contact_nodes": 416, "material_points": 9216,
                             "abaqus_coordinate_matching_maximum_m": match_error,
                             "fuelsim_simultaneous_plastic_creep_points": int(np.count_nonzero(
                                 (actual_increments[:, :, 0] > 0) & (actual_increments[:, :, 1] > 0))),
                             "abaqus_simultaneous_plastic_creep_points": int(np.count_nonzero(
                                 (reference_increments[:, :, 0] > 0) & (reference_increments[:, :, 1] > 0))),
                             "projected_contact_nodes": int(np.sum(nodal("contact_projected_fuel_cladding")[secondary])),
                             "sliding_contact_nodes": int(np.sum(nodal("contact_sliding_fuel_cladding")[secondary]))})
    rows, by_frame = [], []
    for name, (actual, reference, labels, tolerance) in collected.items():
        rows.append(metrics(name, np.concatenate(actual), np.concatenate(reference), sum(labels, []), tolerance, "all_frames"))
        rows.append(metrics(name, actual[-1], reference[-1], labels[-1], tolerance, "final"))
        for frame, time in enumerate(times):
            record = metrics(name, actual[frame], reference[frame], labels[frame], tolerance, "frame")
            by_frame.append({"time": time, **record})
    for filename, records in (("comparison.csv", rows), ("comparison_by_frame.csv", by_frame), ("coverage.csv", coverage)):
        with (output / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    passed = all(row["strict_0_1_percent_passed"] for row in rows)
    print(json.dumps({"completed_increments": 20, "complete_coverage": True, "strict_0_1_percent_passed": passed}))
    for row in rows:
        if row["scope"] == "all_frames":
            print("%s L2=%.8g%% peak=%.8g%% pointwise=%.8g%% zero_abs=%.8g passed=%d" %
                  (row["field"], row["relative_l2_percent"], row["relative_absolute_peak_percent"],
                   row["maximum_pointwise_relative_percent"], row["maximum_zero_reference_absolute_difference"],
                   row["strict_0_1_percent_passed"]))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
