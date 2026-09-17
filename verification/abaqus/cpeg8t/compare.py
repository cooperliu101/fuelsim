"""Run unchanged production cards and compare exported fields with native Abaqus.

Zero masks below follow the prescribed constraints and rectangular symmetry;
they do not clip reference denominators or omit any physical field nodes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np


CASES = {
    "prescribed": "small_prescribed",
    "free_controls": "free_controls",
    "finite_prescribed": "finite_reference_node_probe",
    "contact": "contact",
    "finite_contact": "finite_contact",
}


def verify_references(source, manifest="reference.sha256"):
    for line in (source / manifest).read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        if hashlib.sha256((source / filename).read_bytes()).hexdigest() != expected:
            raise AssertionError("Native reference checksum mismatch: " + filename)


def compare(actual, reference, zeros, absolute, label, report):
    actual, reference = np.asarray(actual), np.asarray(reference)
    zeros = np.broadcast_to(np.asarray(zeros, dtype=bool), actual.shape)
    if actual.shape != reference.shape or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(reference)):
        raise AssertionError(label + ": missing or nonfinite field values")
    errors = actual - reference
    result = {"count": int(actual.size), "maximum_absolute": float(np.max(np.abs(errors)))}
    report[label] = result
    if np.any(zeros):
        zero_error = float(max(np.max(np.abs(actual[zeros])), np.max(np.abs(reference[zeros]))))
        result["zero_absolute"] = zero_error
        if zero_error > absolute:
            raise AssertionError(f"{label}: analytical zero {zero_error} exceeds {absolute}")
    if np.any(~zeros):
        a, r = actual[~zeros], reference[~zeros]
        if np.any(r == 0):
            raise AssertionError(label + ": unexpected zero reference needs an explicit analytical mask")
        result.update(relative_l2=float(np.linalg.norm(a-r)/np.linalg.norm(r)),
                      relative_peak=float(abs(np.max(np.abs(a))-np.max(np.abs(r)))/np.max(np.abs(r))),
                      maximum_pointwise_relative=float(np.max(np.abs((a-r)/r))))
        if max(result[k] for k in ("relative_l2", "relative_peak", "maximum_pointwise_relative")) >= .005:
            raise AssertionError(f"{label}: exceeds 0.5%: {result}")
    report[label] = result


def check(case, directory, source):
    native = json.loads((source / (CASES[case] + "_fields.json")).read_text())
    step = list(native["steps"])[-1]
    fields = native["steps"][step][-1]["fields"]
    contact = case in ("contact", "finite_contact")
    count = 16 if contact else 8
    corners = [1, 2, 3, 4, 9, 10, 11, 12] if contact else [1, 2, 3, 4]
    report = {}
    with netCDF4.Dataset(directory / (case + "_results.e")) as output:
        nodal_names = list(netCDF4.chartostring(output["name_nod_var"][:]))
        def nodal(name):
            return np.asarray(output["vals_nod_var%d" % (nodal_names.index(name)+1)][-1])
        for name, key, component, nodes, zero, absolute in [
            ("temperature", "NT11", None, corners, lambda n: False, 1e-9),
            ("displacement_x", "U", 0, range(1, count+1), lambda n: True, 1e-12),
            ("displacement_y", "U", 1, range(1, count+1), lambda n: not contact or n in (1,2,5), 1e-12),
            ("reaction_x", "RF", 0, range(1, count+1), lambda n: contact and n in (5,7,13,15), 1e-8),
            ("reaction_y", "RF", 1, range(1, count+1), lambda n: contact and n not in (1,2,5,11,12,15), 1e-8),
            ("heat_reaction", "RFL11", None, corners,
             lambda n: case == "prescribed" or (contact and n not in (1,2,11,12)), 1e-8),
        ]:
            values = {v["nodeLabel"]: v["data"] for v in fields[key]["values"]}
            ref = [values[n] if component is None else values[n][component] for n in nodes]
            compare(nodal(name)[np.array(list(nodes))-1], ref, [zero(n) for n in nodes], absolute, name, report)
        element_names = list(netCDF4.chartostring(output["name_elem_var"][:]))
        for component, name in enumerate(("xx", "yy", "zz", "xy")):
            actual, reference = [], []
            for value in fields["S"]["values"]:
                variable = element_names.index("stress_%s_q%d" % (name, value["integrationPoint"]-1))+1
                block = value["elementLabel"]
                actual.append(float(output["vals_elem_var%deb%d" % (variable, block)][-1, 0]))
                reference.append(value["data"][component])
            compare(actual, reference, name == "xy" or (case == "free_controls" and name == "zz"),
                    1e-7, "stress_"+name, report)
        controls = [17, 18] if contact else [100]
        sections = ["lower", "upper"] if contact else ["solid"]
        globals_ = dict(zip(netCDF4.chartostring(output["name_glo_var"][:]), output["vals_glo_var"][-1]))
        histories = []
        for node in controls:
            regions = [r for r in native["history"][step] if r.startswith("Node ") and r.endswith("."+str(node))]
            if len(regions) != 1:
                raise AssertionError("Missing native reference node history")
            histories.append({k: v[-1][1] for k, v in native["history"][step][regions[0]].items()})
        for name, key in [("u3", "U3"), ("rotation_x", "UR1"), ("rotation_y", "UR2"),
                          ("axial_force", "RF3"), ("moment_x", "RM1"), ("moment_y", "RM2")]:
            reference, actual = [], []
            for section, history in zip(sections, histories):
                prefix = "section_" + section + "_"
                x, y = globals_[prefix+"origin_x"], globals_[prefix+"origin_y"]
                value = history[key]
                if key == "U3":
                    value += y*history["UR1"] - x*history["UR2"]
                elif key == "RM1":
                    value -= y*history["RF3"]
                elif key == "RM2":
                    value += x*history["RF3"]
                reference.append(value)
                actual.append(globals_[prefix+name])
            zero = (contact and key not in ("RF3", "RM1")) or (case == "free_controls" and key in ("RF3", "RM1", "RM2"))
            if case == "contact" and key == "RM1":
                zero = True
            compare(actual, reference, zero, 1e-9, name, report)
        if contact:
            globals_ = dict(zip(netCDF4.chartostring(output["name_glo_var"][:]), output["vals_glo_var"][-1]))
            for name, prefix in [("pressure", "CPRESS"), ("gap", "COPEN")]:
                native_key = next(k for k in fields if k.startswith(prefix))
                values = np.array([v["data"] for v in fields[native_key]["values"]])
                compare(values, np.full_like(values, values.mean()), False, 1e-12, "native_uniform_"+name, report)
                compare([globals_["contact_0_q%d_%s" % (q, name)] for q in range(3)],
                        [values.mean()]*3, False, 1e-12, "contact_"+name, report)
            native_key = next(k for k in fields if k.startswith("CNORMF"))
            force = sum(v["data"][1] for v in fields[native_key]["values"] if v["nodeLabel"] in (9,10,13))
            compare([sum(globals_["contact_0_q%d_force" % q] for q in range(3))], [force], False,
                    1e-8, "contact_force", report)
            heat = sum(v["data"] for v in fields["RFL11"]["values"] if v["nodeLabel"] in (11,12))
            compare([sum(value for key,value in globals_.items()
                         if key.startswith('contact_0_q') and key.endswith('_heat_rate'))], [heat], False,
                    1e-8, "contact_heat_rate", report)
    (directory / "comparison.json").write_text(json.dumps(report, indent=2)+"\n")
    maximum = max((v.get("maximum_pointwise_relative", 0) for v in report.values()), default=0)
    print(f"{case}: all nodes, 9 material points per element and section controls passed; maximum relative error={maximum:.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=CASES)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--serial-work", type=Path)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    verify_references(source)
    args.work.mkdir(parents=True, exist_ok=True)
    mesh = "contact.e" if "contact" in args.case else "rectangle.e"
    for filename in (args.case+".fsi", mesh):
        shutil.copyfile(source / filename, args.work / filename)
        if (source / filename).read_bytes() != (args.work / filename).read_bytes():
            raise AssertionError("Input must remain byte-for-byte unchanged")
    command = [str(args.executable.resolve()), "-i", args.case+".fsi"]
    if args.mpiexec:
        command = [str(args.mpiexec), "-n", "2"] + command
    process = subprocess.run(command,
                             cwd=args.work, capture_output=True, text=True)
    (args.work / "production.log").write_text(process.stdout+process.stderr)
    if process.returncode or "completed=true" not in process.stdout:
        raise RuntimeError(process.stdout+process.stderr)
    check(args.case, args.work, source)
    if args.serial_work:
        with netCDF4.Dataset(args.work / (args.case+"_results.e")) as parallel, \
                netCDF4.Dataset(args.serial_work / (args.case+"_results.e")) as serial:
            element_names = list(netCDF4.chartostring(serial["name_elem_var"][:]))
            zero_shear_names = {"stress_xy_q%d" % q for q in range(9)}
            for key in serial.variables:
                if not key.startswith(("vals_nod_var", "vals_elem_var", "vals_glo_var", "time_whole")):
                    continue
                actual, reference = np.asarray(parallel[key][:]), np.asarray(serial[key][:])
                if not np.array_equal(np.isnan(actual), np.isnan(reference)):
                    raise AssertionError("MPI changes field coverage: "+key)
                valid = np.isfinite(reference)
                if args.case == "finite_contact" and key.startswith("vals_elem_var"):
                    variable = int(key.removeprefix("vals_elem_var").split("eb")[0])-1
                    if element_names[variable] in zero_shear_names:
                        # User-authorized exception: this rectangular case has
                        # analytically zero shear stress. All other MPI gates
                        # retain their original relative and absolute limits.
                        maximum = max(np.max(np.abs(actual[valid])), np.max(np.abs(reference[valid])),
                                      np.max(np.abs(actual[valid]-reference[valid])))
                        if not np.isfinite(maximum) or maximum > 1e-7:
                            raise AssertionError("MPI analytical zero shear stress exceeds 1e-7 Pa: "+key)
                        continue
                if not np.allclose(actual[valid], reference[valid], rtol=1e-12, atol=1e-12):
                    raise AssertionError("MPI differs from serial exported fields: "+key)
        print("Two-rank and one-rank fields, all material histories and contact outputs agree")
