#!/usr/bin/env python3
"""Read recorded fields only; compare explicit HEX8 thermal weight candidates.

This is diagnostic algebra, not a production residual or a replacement solver.
Every native candidate uses the same Abaqus nodal temperature/displacement.
"""

import csv
import hashlib
import json
from pathlib import Path

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIGINAL = HERE.parent / "b523_hex8_c3d8t_integrated_path"
NOHEAT = HERE / "b523_integrated_no_friction_heat"
ACTUAL = ROOT / "build/blackbox/b523/verification/fuelsim/transient_b523_integrated_contact_results.e"
SIGNS = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                  [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)
NATIVE_POINTS = [1, 2, 4, 3, 5, 6, 8, 7]
BOUNDARY = [1, 4, 7, 10, 15, 18, 21, 24]


def read_csv(path):
    with path.open(newline="") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def write_tsv(name, rows):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def geometry(path):
    coordinates, connectivity = {}, {}
    section = ""
    for line in path.read_text().splitlines():
        if line.startswith("**"):
            continue
        if line.startswith("*"):
            section = line.split(",")[0].lower()
        elif section in ("*node", "*element"):
            values = [float(value) for value in line.split(",")]
            label = int(values[0])
            if section == "*node":
                coordinates[label] = values[1:]
            else:
                connectivity[label] = np.array(values[1:], dtype=int) - 1
    assert sorted(coordinates) == list(range(1, 25))
    assert sorted(connectivity) == list(range(1, 5))
    return np.array([coordinates[node] for node in sorted(coordinates)]), connectivity


def shape(natural):
    factors = 1.0 + SIGNS * natural
    values = np.prod(factors, axis=1) / 8.0
    derivatives = np.empty((8, 3))
    for direction in range(3):
        other = [axis for axis in range(3) if axis != direction]
        derivatives[:, direction] = SIGNS[:, direction] * np.prod(factors[:, other], axis=1) / 8.0
    return values, derivatives


def metrics(coordinates, natural):
    values, derivatives = shape(natural)
    jacobian = coordinates.T @ derivatives
    determinant = np.linalg.det(jacobian)
    assert determinant > 0.0
    return values, derivatives @ np.linalg.inv(jacobian), determinant


def load_native(prefix):
    nodes = read_csv(Path(str(prefix) + "_nodal.csv"))
    points = read_csv(Path(str(prefix) + "_integration.csv"))
    energy = read_csv(Path(str(prefix) + "_energy.csv"))
    assert len(nodes) == 480 and len(points) == 640 and len(energy) == 20
    return ({(int(row["increment"]), int(row["node"])): row for row in nodes},
            {(int(row["increment"]), int(row["element"]), int(row["integration_point"])): row for row in points},
            energy)


def load_actual():
    with netCDF4.Dataset(ACTUAL) as dataset:
        names = netCDF4.chartostring(dataset["name_nod_var"][:]).tolist()
        fields = {name: np.array(dataset["vals_nod_var%d" % (names.index(name) + 1)][:])
                  for name in ["temperature", "displacement_x", "displacement_y", "displacement_z", "reaction_heat_flux"]}
        times = np.array(dataset["time_whole"][:])
        assert np.max(np.abs(times - np.arange(21) * 0.02)) < 1.0e-12
        return fields


def main():
    reference, elements = geometry(Path(str(NOHEAT) + ".inp"))
    actual = load_actual()
    all_candidates, energy_rows, point_rows, decomposition, volume_rows = [], [], [], [], []
    capacity_rules = ["current_gauss", "current_corner", "reference_gauss_volume_ratio",
                      "reference_corner_volume_ratio", "current_row_sum", "native_ivol"]
    conduction_rules = ["current_physical", "current_volume_ratio", "midpoint_physical",
                        "midpoint_volume_ratio", "native_hfl_current", "native_hfl_midpoint"]
    for case, prefix in [("original_friction_heat", ORIGINAL), ("no_friction_heat", NOHEAT)]:
        nodes, points, energies = load_native(prefix)
        old_temperatures = np.full(24, 300.0)
        old_positions = reference.copy()
        old_friction = 0.0
        for increment in range(1, 21):
            temperatures = np.array([nodes[increment, node]["temperature_k"] for node in range(1, 25)])
            positions = reference + np.array([[nodes[increment, node]["u%d_m" % direction]
                                               for direction in range(1, 4)] for node in range(1, 25)])
            native_reactions = np.array([nodes[increment, node]["reaction_heat_flux_w"] for node in range(1, 25)])
            capacities = {rule: np.zeros(24) for rule in capacity_rules}
            conductions = {rule: np.zeros(24) for rule in conduction_rules}
            for element, indices in elements.items():
                ref = reference[indices]
                current = positions[indices]
                midpoint = 0.5 * (current + old_positions[indices])
                local_temperature = temperatures[indices]
                rate = (local_temperature - old_temperatures[indices]) / 0.02
                specific_heat = 1.0 + 0.001 * (local_temperature - 300.0)
                ref_metrics = [metrics(ref, sign / np.sqrt(3.0)) for sign in SIGNS]
                cur_metrics = [metrics(current, sign / np.sqrt(3.0)) for sign in SIGNS]
                mid_metrics = [metrics(midpoint, sign / np.sqrt(3.0)) for sign in SIGNS]
                ratio = sum(metric[2] for metric in cur_metrics) / sum(metric[2] for metric in ref_metrics)
                row_sum = sum(metric[0] * metric[2] for metric in cur_metrics)
                if case == "no_friction_heat":
                    reference_volume = sum(metric[2] for metric in ref_metrics)
                    average_gradient = sum(current.T @ metric[1] * metric[2] for metric in ref_metrics) / reference_volume
                    native_weights = [points[increment, element, point]["ivol_m3"] for point in NATIVE_POINTS]
                    volume_rows.append(dict(increment=increment, time_s=increment * 0.02, element=element,
                                            reference_volume_m3=reference_volume,
                                            current_volume_m3=sum(metric[2] for metric in cur_metrics),
                                            native_ivol_sum_m3=sum(native_weights), current_reference_ratio=ratio,
                                            determinant_average_f=np.linalg.det(average_gradient),
                                            maximum_ivol_weight_difference_m3=max(abs(native_weights[index] - metric[2] * ratio)
                                                                                 for index, metric in enumerate(ref_metrics))))
                for local_node, global_index in enumerate(indices):
                    native = points[increment, element, NATIVE_POINTS[local_node]]
                    assert abs(native["temperature_k"] - local_temperature[local_node]) < 1.0e-6
                    weights = {"current_gauss": cur_metrics[local_node][2],
                               "current_corner": metrics(current, SIGNS[local_node])[2],
                               "reference_gauss_volume_ratio": ref_metrics[local_node][2] * ratio,
                               "reference_corner_volume_ratio": metrics(ref, SIGNS[local_node])[2] * ratio,
                               "current_row_sum": row_sum[local_node], "native_ivol": native["ivol_m3"]}
                    for rule, weight in weights.items():
                        capacities[rule][global_index] += weight * 100.0 * specific_heat[local_node] * rate[local_node]
                    if case == "no_friction_heat":
                        point_rows.append(dict(increment=increment, time_s=increment * 0.02,
                                               element=element, node=global_index + 1, **weights))
                    conductivity_base = 15.0 if element <= 2 else 10.0
                    conductivity = conductivity_base * (1.0 + 0.001 * (local_temperature[local_node] - 300.0))
                    native_flux = np.array([native["hfl%d_w_m2" % direction] for direction in range(1, 4)])
                    for name, metric in [("current", cur_metrics[local_node]), ("midpoint", mid_metrics[local_node])]:
                        gradient = metric[1]
                        temperature_gradient = local_temperature @ gradient
                        physical = conductivity * (gradient @ temperature_gradient)
                        conductions[name + "_physical"][indices] += metric[2] * physical
                        conductions[name + "_volume_ratio"][indices] += ref_metrics[local_node][2] * ratio * physical
                        conductions["native_hfl_" + name][indices] -= native["ivol_m3"] * (gradient @ native_flux)
            energy = energies[increment - 1]
            friction_power = (energy["allfd_j"] - old_friction) / 0.02
            energy_rows.append(dict(case=case, increment=increment, time_s=increment * 0.02,
                                   native_ivol_storage_w=float(sum(capacities["native_ivol"])),
                                   native_boundary_w=energy["boundary_heat_rate_w"],
                                   friction_power_w=friction_power,
                                   thermal_balance_w=float(sum(capacities["native_ivol"])) - energy["boundary_heat_rate_w"]
                                                     - (friction_power if case == "original_friction_heat" else 0.0)))
            for node in BOUNDARY:
                index = node - 1
                for capacity_rule in capacity_rules:
                    for conduction_rule in conduction_rules:
                        prediction = capacities[capacity_rule][index] + conductions[conduction_rule][index]
                        all_candidates.append(dict(case=case, increment=increment, time_s=increment * 0.02, node=node,
                                                   capacity_rule=capacity_rule, conduction_rule=conduction_rule,
                                                   capacity_w=capacities[capacity_rule][index],
                                                   conduction_w=conductions[conduction_rule][index],
                                                   predicted_rfl_w=prediction, native_rfl_w=native_reactions[index],
                                                   prediction_minus_native_w=prediction - native_reactions[index]))
                if case == "no_friction_heat":
                    full_difference = actual["reaction_heat_flux"][increment, index] - native_reactions[index]
                    measure_difference = capacities["current_gauss"][index] - capacities["native_ivol"][index]
                    decomposition.append(dict(increment=increment, time_s=increment * 0.02, node=node,
                                              actual_rfl_w=actual["reaction_heat_flux"][increment, index],
                                              native_rfl_w=native_reactions[index], full_difference_w=full_difference,
                                              native_ivol_capacity_w=capacities["native_ivol"][index],
                                              current_gauss_capacity_w=capacities["current_gauss"][index],
                                              current_physical_conduction_w=conductions["current_physical"][index],
                                              current_volume_ratio_conduction_w=conductions["current_volume_ratio"][index],
                                              candidate_minus_native_w=capacities["reference_gauss_volume_ratio"][index]
                                                                       + conductions["current_volume_ratio"][index]
                                                                       - native_reactions[index],
                                              capacity_measure_difference_w=measure_difference,
                                              remaining_difference_w=full_difference - measure_difference))
            old_temperatures, old_positions, old_friction = temperatures, positions, energy["allfd_j"]
    write_tsv("thermal_energy_audit.tsv", energy_rows)
    write_tsv("thermal_weight_candidates.tsv", point_rows)
    summary = []
    for case in ["original_friction_heat", "no_friction_heat"]:
        for capacity_rule in capacity_rules:
            for conduction_rule in conduction_rules:
                selected = [row for row in all_candidates if row["case"] == case
                            and row["capacity_rule"] == capacity_rule and row["conduction_rule"] == conduction_rule]
                peak = max(selected, key=lambda row: abs(row["prediction_minus_native_w"]))
                summary.append(dict(case=case, capacity_rule=capacity_rule, conduction_rule=conduction_rule,
                                    samples=len(selected), maximum_absolute_error_w=abs(peak["prediction_minus_native_w"]),
                                    maximum_absolute_increment=peak["increment"], maximum_absolute_node=peak["node"],
                                    relative_l2=np.linalg.norm([row["prediction_minus_native_w"] for row in selected])
                                                / np.linalg.norm([row["native_rfl_w"] for row in selected])))
    write_tsv("thermal_reaction_candidates.tsv", summary)
    write_tsv("thermal_reaction_decomposition.tsv", decomposition)
    write_tsv("thermal_volume_audit.tsv", volume_rows)
    paths = [ACTUAL, Path(__file__)]
    for prefix in (ORIGINAL, NOHEAT):
        paths.extend([Path(str(prefix) + suffix) for suffix in (".inp", "_nodal.csv", "_integration.csv", "_energy.csv")])
    provenance = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    (HERE / "thermal_audit_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    for capacity_rule in capacity_rules:
        for conduction_rule in conduction_rules:
            values = [abs(row["prediction_minus_native_w"]) for row in all_candidates
                      if row["case"] == "no_friction_heat" and row["capacity_rule"] == capacity_rule
                      and row["conduction_rule"] == conduction_rule]
            print(capacity_rule, conduction_rule, "maximum_absolute_w", max(values))
    print("final hot-node decomposition", [row for row in decomposition if row["increment"] == 20 and row["node"] >= 15])
    print("energy closure", {case: max(abs(row["thermal_balance_w"]) for row in energy_rows if row["case"] == case)
                             for case in ["original_friction_heat", "no_friction_heat"]})


if __name__ == "__main__":
    main()
