"""Read native Abaqus CSVs; identify material temperature without a Fuelsim solve."""
import argparse
import csv
import hashlib
import math
from pathlib import Path

from analyze_plastic import analyze as analyze_plastic
from analyze_creep import analyze as analyze_creep


CASES = ("elastic_moduli", "elastic_moduli_finite", "plastic_yield", "plastic_yield_finite",
         "plastic_hardening", "plastic_hardening_finite", "creep_coefficient",
         "creep_coefficient_finite", "creep_exponent", "creep_exponent_finite")
RULES = ("corner", "gauss_interpolated", "element_average")
POINT_NODE = {1: 1, 2: 2, 3: 4, 4: 3}
NODE_SIGNS = ((-1, -1), (1, -1), (1, 1), (-1, 1))
GATE = 0.001


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_provenance(directory):
    extractor_hash = hash_file(directory / "extract.py")
    for case in CASES:
        text = (directory / (case + "_provenance.txt")).read_text(encoding="utf-8-sig")
        values = dict(line.split("=", 1) for line in text.splitlines() if "=" in line)
        if values.get("input_sha256", "").lower() != hash_file(directory / (case + ".inp")):
            raise AssertionError(case + ": input hash differs from native execution")
        if values.get("extractor_sha256", "").lower() != extractor_hash:
            raise AssertionError(case + ": extractor hash differs from native execution")
        if values.get("analysis_exit_code") != "0" or values.get("cpus") != "1":
            raise AssertionError(case + ": native execution contract failed")
        if b"THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in (directory / (case + ".sta")).read_bytes():
            raise AssertionError(case + ": no native completion marker")
        if b"Abaqus 3DEXPERIENCE R2018x" not in (directory / (case + ".dat")).read_bytes():
            raise AssertionError(case + ": native version differs from this qualification")
        if case.startswith("creep_"):
            messages = (directory / (case + ".msg")).read_bytes()
            expected = b"IMPLICIT" if case.endswith("_finite") else b"EXPLICIT"
            other = b"EXPLICIT" if case.endswith("_finite") else b"IMPLICIT"
            if messages.count(expected + b" TIME INTEGRATION WILL BE USED") != 2 or \
                    other + b" TIME INTEGRATION WILL BE USED" in messages:
                raise AssertionError(case + ": native creep time integration differs from the checked relation")
        complete = (directory / (case + "_export_complete.txt")).read_text()
        if "position=INTEGRATION_POINT" not in complete or "component_order=11,22,33,12" not in complete:
            raise AssertionError(case + ": material-point extraction not verified")
        counts = ({"FIRST": 4, "SWAP": 4} if case.startswith("elastic_") else
                  {"LOAD_FIXED_T": 5, "LOAD_SWAPPED_T": 5} if case.startswith("plastic_") else
                  {"PRELOAD": 1, "RELAX_FAST": 4, "RELAX_SLOW": 4})
        expected_frames = {(step, frame) for step, count in counts.items() for frame in range(count + 1)}
        if "accepted_frames=%d\n" % sum(counts.values()) not in complete:
            raise AssertionError(case + ": native increment count changed")
        for kind, label in (("nodes", "node"), ("points", "point")):
            frames = {}
            for row in read_csv(directory / (case + "_" + kind + ".csv")):
                key = (row["step"], int(row["frame"]))
                identity = int(row[label])
                if identity in frames.setdefault(key, set()):
                    raise AssertionError(case + ": duplicate exported row")
                frames[key].add(identity)
            if set(frames) != expected_frames or any(ids != {1, 2, 3, 4} for ids in frames.values()):
                raise AssertionError(case + ": incomplete native frame/node/point coverage")


def tensor_inner(a, b):
    return sum(a[i] * b[i] for i in range(3)) + 2 * a[3] * b[3]


def temperature_candidates(point, temperatures):
    sx, sy = NODE_SIGNS[POINT_NODE[point] - 1]
    xi, eta = sx / math.sqrt(3), sy / math.sqrt(3)
    interpolated = sum(temperatures[n + 1] * (1 + nx * xi) * (1 + ny * eta) / 4
                       for n, (nx, ny) in enumerate(NODE_SIGNS))
    return {"corner": temperatures[POINT_NODE[point]], "gauss_interpolated": interpolated,
            "element_average": sum(temperatures.values()) / 4}


def analyze_elastic(directory):
    metrics = []
    details = []
    for case in CASES[:2]:
        nodes = {}
        for row in read_csv(directory / (case + "_nodes.csv")):
            key = (row["step"], int(row["frame"]))
            node = int(row["node"])
            if node in nodes.setdefault(key, {}):
                raise AssertionError("Duplicate elastic nodal sample")
            temperature = float(row["temperature"])
            if not math.isfinite(temperature):
                raise AssertionError(case + ": nonfinite native temperature")
            initial = {1: 600.0, 2: 650.0, 3: 700.0, 4: 750.0}
            final = {1: 750.0, 2: 700.0, 3: 650.0, 4: 600.0}
            # The native default temperature boundary changes immediately at
            # the start of SWAP; frame zero still holds the preceding state.
            expected = final[node] if row["step"] == "SWAP" and int(row["frame"]) > 0 else initial[node]
            if abs(temperature - expected) > 1e-6:
                raise AssertionError(case + ": prescribed elastic temperature differs from native output")
            nodes[key][node] = temperature
        samples = {}
        errors = {(parameter, rule): [] for parameter in ("young_modulus", "poisson_ratio") for rule in RULES}
        for row in read_csv(directory / (case + "_points.csv")):
            frame = int(row["frame"])
            if not frame:
                continue
            key = (row["step"], frame)
            point = int(row["point"])
            if int(row["element"]) != 1 or point in samples.setdefault(key, set()):
                raise AssertionError("Duplicate or unexpected elastic point")
            samples[key].add(point)
            temperatures = nodes[key]
            if set(temperatures) != {1, 2, 3, 4}:
                raise AssertionError("Incomplete elastic temperature coverage")
            stress = [float(row["s_" + c]) for c in ("rr", "zz", "hoop", "rz")]
            strain = [float(row["ee_" + c]) for c in ("rr", "zz", "hoop", "rz")]
            if not all(math.isfinite(value) for value in stress + strain):
                raise AssertionError(case + ": nonfinite native stress or elastic strain")
            stress_trace, strain_trace = sum(stress[:3]), sum(strain[:3])
            ds = [stress[i] - (stress_trace / 3 if i < 3 else 0) for i in range(4)]
            de = [strain[i] - (strain_trace / 3 if i < 3 else 0) for i in range(4)]
            if abs(strain_trace) < 1e-6 or tensor_inner(de, de) < 1e-10:
                raise AssertionError("Elastic sample cannot independently identify bulk and shear moduli")
            bulk = stress_trace / (3 * strain_trace)
            shear = tensor_inner(ds, de) / (2 * tensor_inner(de, de))
            young = 9 * bulk * shear / (3 * bulk + shear)
            poisson = (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear))
            if not math.isfinite(young) or not 0 < poisson < 0.5:
                raise AssertionError("Invalid inferred elastic parameters")
            detail = {"case": case, "step": row["step"], "frame": frame, "time": row["time"],
                      "point": point, "paired_node": POINT_NODE[point],
                      "inferred_young_modulus": young, "inferred_poisson_ratio": poisson}
            for rule, temperature in temperature_candidates(point, temperatures).items():
                expected_young = 1e9 - 2e6 * (temperature - 600)
                expected_poisson = 0.2 + 0.001 * (temperature - 600)
                detail[rule + "_temperature"] = temperature
                for name, measured, expected in (("young_modulus", young, expected_young),
                                                 ("poisson_ratio", poisson, expected_poisson)):
                    errors[name, rule].append((abs(measured - expected), abs(measured - expected) / abs(expected)))
                if rule == "corner":
                    k = expected_young / (3 * (1 - 2 * expected_poisson))
                    g = expected_young / (2 * (1 + expected_poisson))
                    predicted = [2 * g * de[i] + (k * strain_trace if i < 3 else 0) for i in range(4)]
                    norm = math.sqrt(tensor_inner(stress, stress))
                    if max(abs(a - b) for a, b in zip(stress, predicted)) / norm >= GATE:
                        raise AssertionError("Corner rule fails full elastic stress reconstruction")
            details.append(detail)
        if len(samples) != 8 or any(points != {1, 2, 3, 4} for points in samples.values()):
            raise AssertionError(case + ": expected all eight increments and four material points")
        for (parameter, rule), values in errors.items():
            maximum = max(relative for _, relative in values)
            accepted = maximum < GATE
            if accepted != (rule == "corner"):
                raise AssertionError(case + ": ambiguous or rejected elastic temperature rule " + rule)
            metrics.append({"case": case, "parameter": parameter, "rule": rule, "samples": len(values),
                            "max_absolute_error": max(absolute for absolute, _ in values),
                            "max_relative_error": maximum, "pass": accepted})
    write_csv(directory / "elastic_identification.csv", details)
    return metrics


def write_csv(path, rows):
    columns = list(dict.fromkeys(column for row in rows for column in row))
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    directory = parser.parse_args().directory
    verify_provenance(directory)
    rows = analyze_elastic(directory) + analyze_plastic(directory) + analyze_creep(directory)
    write_csv(directory / "metrics.csv", rows)
    for row in rows:
        parameter = row.get("parameter", "constitutive_relation")
        print("%s %s %s samples=%s max_relative_error=%.9g%% pass=%s" %
              (row["case"], parameter, row["rule"], row["samples"],
               100 * row["max_relative_error"], row["pass"]))
    print("PASS: all ten native cases identify paired-corner temperatures; alternatives rejected")


if __name__ == "__main__":
    main()
