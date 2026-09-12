"""Check native references against independent uniform/operator equations.

This reads exported Abaqus results only. It does not solve Fuelsim or infer a
full-field equivalence between the CAX4T and radial GPS discretizations.
"""

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MAXIMUM_RELATIVE = {}
MAXIMUM_ZERO_ABSOLUTE = {}
RELATIVE_TOLERANCE = 1e-6


def rows(case, kind):
    with (ROOT / (case + "_" + kind + ".csv")).open() as stream:
        result = list(csv.DictReader(stream))
    for row in result:
        for key, value in row.items():
            if key != "field" and not math.isfinite(float(value)):
                raise AssertionError("Nonfinite reference value: " + case + ":" + key)
    return result


def check(name, actual, expected, zero_tolerance):
    difference = abs(actual - expected)
    if expected == 0:
        MAXIMUM_ZERO_ABSOLUTE[name] = max(MAXIMUM_ZERO_ABSOLUTE.get(name, 0.0), difference)
        if difference > zero_tolerance:
            raise AssertionError("%s: absolute difference %.17g exceeds %.17g" % (name, difference, zero_tolerance))
    else:
        error = difference / abs(expected)
        MAXIMUM_RELATIVE[name] = max(MAXIMUM_RELATIVE.get(name, 0.0), error)
        if error >= RELATIVE_TOLERANCE:
            raise AssertionError("%s: relative difference %.17g exceeds %.17g" % (name, error, RELATIVE_TOLERANCE))


def frames(records, key):
    result = {}
    for row in records:
        time = float(row["time"])
        index = tuple(int(row[name]) for name in key)
        if index in result.setdefault(time, {}):
            raise AssertionError("Duplicate reference record")
        result[time][index] = row
    if len(result) != 10:
        raise AssertionError("Missing accepted reference increments")
    return result


def uniform(case, finite):
    nodes = frames(rows(case, "nodes"), ("node",))
    points = frames(rows(case, "points"), ("element", "point"))
    if set(nodes) != set(points):
        raise AssertionError("Node and material histories differ")
    for increment, time in enumerate(sorted(nodes), 1):
        check("uniform_time", time, increment / 10.0, 1e-16)
        ns, ps = nodes[time], points[time]
        if set(ns) != {(1,), (2,), (3,), (4,)} or set(ps) != {(1, q) for q in range(1, 5)}:
            raise AssertionError("Incomplete single-annulus node or material-point coverage")
        temperature = 600.0 + 5.0 * increment
        load = 2.8274333882308138 * increment
        axial_length = 0.01 + (float(ns[3,]["uz"]) if finite else 0.0)
        stress_force = sum(float(row["stress_zz"]) * float(row["volume"]) for row in ps.values()) / axial_length
        prefix = "finite_" if finite else "small_"
        check(prefix + "axial_stress_resultant", stress_force, load, 1e-8)
        check(prefix + "bottom_force", sum(float(ns[n,]["rf_z"]) for n in (1, 2)), -load, 1e-8)
        check(prefix + "uniform_heat_capacity",
              sum(float(row["reaction_heat"]) for row in ns.values()),
              1e6 * 50.0 * sum(float(row["volume"]) for row in ps.values()), 1e-10)
        for row in ns.values():
            check(prefix + "temperature", float(row["temperature"]), temperature, 1e-11)
            check(prefix + "radial_external_reaction", float(row["rf_r"]), 0.0, 1e-8)
            if not finite:
                scale = increment / 10.0
                check("small_radial_displacement", float(row["ur"]), float(row["r"]) * 2e-4 * scale, 1e-13)
                check("small_axial_displacement", float(row["uz"]), float(row["z"]) * 1.5e-3 * scale, 1e-13)
        for row in ps.values():
            check(prefix + "heat_flux_r", float(row["heat_flux_r"]), 0.0, 1e-7)
            check(prefix + "heat_flux_z", float(row["heat_flux_z"]), 0.0, 1e-7)
            if not finite:
                for component, factor in (("rr", 0.0), ("zz", 1.0), ("hoop", 0.0), ("rz", 0.0)):
                    check("small_stress_" + component, float(row["stress_" + component]), factor * increment * 1e5, 1e-4)
                for component, factor in (("rr", -0.3), ("zz", 1.0), ("hoop", -0.3), ("rz", 0.0)):
                    check("small_elastic_" + component, float(row["elastic_" + component]), factor * increment * 1e-4, 1e-13)
    print(case + "=passed")


def contact_operator():
    case = "gps_two_slice_contact"
    nodes = frames(rows(case, "nodes"), ("node",))
    points = frames(rows(case, "points"), ("element", "point"))
    native = {}
    for row in rows(case, "contact"):
        short, pair = row["field"].split()
        key = (float(row["time"]), pair, short, int(row["node"]), int(row["component"]))
        if key in native:
            raise AssertionError("Duplicate native contact value")
        native[key] = float(row["value"])
    displacements = (0.0, 0.0, 1e-8, 2e-6, 1.999e-6, -2e-6, -2e-6, 2e-5, 2e-5, 2.001e-5, 2.4e-5)
    summaries = []
    for layer, height, secondary, primary, hot, cold in (
        ("ONE", 0.01, (2, 3), (5, 8), (1, 4), (6, 7)),
        ("TWO", 0.02, (10, 11), (13, 16), (9, 12), (14, 15))):
        pair = "SECONDARY_" + layer + "/PRIMARY_" + layer
        area = 2 * math.pi * 0.005 * height
        allowable_elastic_slip = 1e-5 * height
        elastic, total = 0.0, 0.0
        for increment, time in enumerate(sorted(nodes), 1):
            check("contact_time", time, float(increment), 1e-16)
            if set(nodes[time]) != {(n,) for n in range(1, 17)}:
                raise AssertionError("Incomplete two-slice node coverage")
            if set(points[time]) != {(e, q) for e in range(1, 5) for q in range(1, 5)}:
                raise AssertionError("Incomplete two-slice material-point coverage")
            pressure = 0.0 if increment in (6, 7) else 1e7
            gap = 1e-5 if pressure == 0 else -1e-5
            if pressure == 0:
                elastic = 0.0
            else:
                delta = displacements[increment] - displacements[increment - 1]
                total += delta
                elastic = max(-allowable_elastic_slip, min(allowable_elastic_slip, elastic + delta))
            traction = 0.3 * pressure * elastic / allowable_elastic_slip
            get = lambda short, node, component=0: native[time, pair, short, node, component]
            for node in secondary:
                check("contact_pressure", get("CPRESS", node), pressure, 1e-4)
                check("contact_gap", get("COPEN", node), gap, 1e-13)
                check("contact_signed_traction", get("CSHEAR1", node), traction, 1e-4)
                check("contact_accumulated_slip", get("CSLIP1", node), total, 1e-13)
                flux = 1e4 * (float(nodes[time][node,]["temperature"])
                              - float(nodes[time][primary[secondary.index(node)],]["temperature"]))
                check("contact_heat_flux_law", get("HFL", node), flux, 1e-7)
            for component in (0, 1):
                check("contact_normal_force_conservation",
                      sum(get("CNORMF", node, component) for node in secondary + primary), 0.0, 1e-8)
                check("contact_tangential_force_conservation",
                      sum(get("CSHEARF", node, component) for node in secondary + primary), 0.0, 1e-8)
            normal_force = sum(get("CNORMF", node, 0) for node in secondary)
            tangential_force = sum(get("CSHEARF", node, 1) for node in secondary)
            heat = 0.5 * area * sum(get("HFL", node) for node in secondary)
            check("contact_normal_reference_area", normal_force, -pressure * area, 1e-8)
            check("contact_friction_reference_area", tangential_force, -traction * area, 1e-8)
            hot_heat = sum(float(nodes[time][node,]["reaction_heat"]) for node in hot)
            cold_heat = sum(float(nodes[time][node,]["reaction_heat"]) for node in cold)
            check("contact_heat_reaction_balance", hot_heat + cold_heat, 0.0, 1e-8)
            check("contact_heat_reference_area", heat, hot_heat, 1e-8)
            # Exact linear-radial FE conduction stiffness, not the logarithmic continuum solution.
            inner_stiffness = 2 * math.pi * 10 * height * 0.0045 / 0.001
            outer_stiffness = 2 * math.pi * 10 * height * 0.005505 / 0.00099
            exact_discrete_heat = 100 / (1 / inner_stiffness + 1 / (1e4 * area) + 1 / outer_stiffness)
            check("contact_discrete_heat_network", heat, exact_discrete_heat, 1e-8)
            summaries.append([time, layer, height, area, allowable_elastic_slip,
                              get("CPRESS", secondary[0]), get("COPEN", secondary[0]),
                              get("CSHEAR1", secondary[0]), get("CSLIP1", secondary[0]),
                              normal_force, tangential_force, heat])
    with (ROOT / "gps_two_slice_contact_summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time", "layer", "reference_height", "secondary_reference_area", "maximum_elastic_slip",
                         "pressure", "gap", "signed_internal_traction", "accumulated_slip",
                         "secondary_normal_force_r", "secondary_tangential_force_z", "heat_rate_secondary_to_primary"])
        writer.writerows(sorted(summaries))
    print(case + "_operator=passed")


uniform("gps_uniform_axial", False)
uniform("gps_uniform_axial_finite", True)
contact_operator()
for name, value in sorted(MAXIMUM_RELATIVE.items()):
    print(name + "_maximum_relative_error=%.17g" % value)
for name, value in sorted(MAXIMUM_ZERO_ABSOLUTE.items()):
    print(name + "_maximum_zero_reference_absolute_difference=%.17g" % value)
print("gps_native_reference_checks=passed")
