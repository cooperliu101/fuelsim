#!/usr/bin/env python3
# Unified workflow for PCMI parallel scaling and MOOSE/Fuelsim consistency checks.

import argparse
import json
import os
import re
import shlex
import statistics
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import netCDF4 as nc
import numpy as np


DEFAULT_RANKS = [1, 2, 4, 8]


def decode_name(raw_name: np.ndarray) -> str:
    return bytes(raw_name).decode("utf-8", errors="ignore").replace("\x00", "").strip()


def nodal_name_map(dataset: nc.Dataset) -> Dict[str, int]:
    names = dataset.variables["name_nod_var"][:]
    mapping = {}
    for index, raw_name in enumerate(names, start=1):
        name = decode_name(raw_name)
        if name:
            mapping[name] = index
    return mapping


def nodal_field(dataset: nc.Dataset, mapping: Dict[str, int], field: str,
                step: int) -> np.ndarray:
    if field not in mapping:
        raise KeyError(f"Missing nodal field '{field}' in Exodus file")
    var_name = f"vals_nod_var{mapping[field]}"
    return np.array(dataset.variables[var_name][step], dtype=float)


def select_contact_field(
    mapping: Dict[str, int], requested: Optional[str], solver: str
) -> Optional[str]:
    if requested is not None:
        if requested not in mapping:
            raise KeyError(
                f"Requested {solver} contact pressure field '{requested}' is missing"
            )
        return requested
    candidates = sorted(
        name for name in mapping if name.startswith("contact_pressure")
    )
    if len(candidates) > 1:
        raise RuntimeError(
            f"Multiple {solver} contact pressure fields found; select one with "
            f"--{solver}-contact: {', '.join(candidates)}"
        )
    return candidates[0] if candidates else None


def exodus_contact_field(
    path: Path, requested: Optional[str], solver: str
) -> Optional[str]:
    with nc.Dataset(path) as dataset:
        return select_contact_field(nodal_name_map(dataset), requested, solver)


def compute_metrics(actual: np.ndarray, reference: np.ndarray, mask=None) -> Dict[str, object]:
    if mask is not None:
        actual = np.asarray(actual, dtype=float)[mask]
        reference = np.asarray(reference, dtype=float)[mask]
    else:
        actual = np.asarray(actual, dtype=float)
        reference = np.asarray(reference, dtype=float)

    if actual.shape != reference.shape:
        raise ValueError(
            f"Metric arrays have different shapes: {actual.shape} vs {reference.shape}"
        )

    total_count = int(actual.size)
    actual_nonfinite_count = int(np.count_nonzero(~np.isfinite(actual)))
    reference_nonfinite_count = int(np.count_nonzero(~np.isfinite(reference)))
    finite_pairs = np.isfinite(actual) & np.isfinite(reference)
    excluded_nonfinite_count = int(np.count_nonzero(~finite_pairs))
    actual = actual[finite_pairs]
    reference = reference[finite_pairs]

    if actual.size == 0:
        return {
            "count": total_count,
            "compared_count": 0,
            "actual_nonfinite_count": actual_nonfinite_count,
            "reference_nonfinite_count": reference_nonfinite_count,
            "excluded_nonfinite_count": excluded_nonfinite_count,
            "relative_l2": float("nan"),
            "relative_abs_peak": float("nan"),
            "max_pointwise_relative": float("nan"),
            "max_abs_diff": float("nan"),
            "zero_reference_count": 0,
            "max_zero_reference_abs_diff": float("nan"),
        }

    diff = actual - reference
    abs_diff = np.abs(diff)
    abs_ref = np.abs(reference)

    ref_norm2 = float(np.dot(reference, reference))
    rel_l2 = float(np.sqrt(float(np.dot(diff, diff) / ref_norm2))) if ref_norm2 > 0.0 else float("nan")

    max_ref = float(np.max(abs_ref))
    rel_abs_peak = float(np.max(abs_diff) / max_ref) if max_ref > 0.0 else float("nan")

    finite = abs_ref > 0.0
    if np.any(finite):
        rel = abs_diff[finite] / abs_ref[finite]
        max_pointwise = float(np.max(rel))
    else:
        max_pointwise = float("nan")

    zero_mask = (abs_ref == 0.0)
    zero_count = int(np.count_nonzero(zero_mask))
    max_zero_abs = float(np.max(abs_diff[zero_mask])) if zero_count > 0 else 0.0

    return {
        "count": total_count,
        "compared_count": int(actual.size),
        "actual_nonfinite_count": actual_nonfinite_count,
        "reference_nonfinite_count": reference_nonfinite_count,
        "excluded_nonfinite_count": excluded_nonfinite_count,
        "relative_l2": rel_l2,
        "relative_abs_peak": rel_abs_peak,
        "max_pointwise_relative": max_pointwise,
        "max_abs_diff": float(np.max(abs_diff)),
        "zero_reference_count": zero_count,
        "max_zero_reference_abs_diff": max_zero_abs,
    }


def replace_fuelsim_mesh_with_absolute_path(original: str, input_path: Path) -> str:
    lines = original.splitlines(keepends=True)
    mesh_start = None
    for index, line in enumerate(lines):
        if re.match(r"^\s*\[Mesh\]\s*$", line):
            mesh_start = index
            break
    if mesh_start is None:
        raise ValueError(f"Missing [Mesh] section in fuelsim input: {input_path}")

    mesh_end = mesh_start + 1
    while mesh_end < len(lines) and not re.match(r"^\s*\[\]\s*$", lines[mesh_end]):
        mesh_end += 1
    for index in range(mesh_start + 1, mesh_end):
        match = re.match(r"^(\s*)file\s*=\s*(.*?)\s*$", lines[index].rstrip("\n"))
        if match is None:
            continue
        configured = Path(match.group(2))
        resolved = configured if configured.is_absolute() else input_path.resolve().parent / configured
        lines[index] = f"{match.group(1)}file = {resolved.resolve()}\n"
        return "".join(lines)
    raise ValueError(f"Missing Mesh.file in fuelsim input: {input_path}")


def configure_fuelsim_outputs(
    original: str, output_exodus: Optional[Path]
) -> str:
    lines = original.splitlines(keepends=True)
    outputs_start = None
    for index, line in enumerate(lines):
        if re.match(r"^\s*\[Outputs\]\s*$", line):
            outputs_start = index
            break

    if outputs_start is None:
        if original and not original.endswith("\n"):
            original += "\n"
        block = "\n[Outputs]\n  console = false\n"
        if output_exodus is not None:
            block += f"  exodus = {output_exodus}\n  exodus_interval = 1\n"
        return original + block + "[]\n"

    outputs_end = outputs_start + 1
    while outputs_end < len(lines) and not re.match(r"^\s*\[\]\s*$", lines[outputs_end]):
        outputs_end += 1

    rewritten = lines[: outputs_start + 1]
    console_written = False
    file_keys = {
        "csv",
        "exodus",
        "exodus_interval",
        "history",
        "history_interval",
        "checkpoint",
        "checkpoint_interval",
    }
    for line in lines[outputs_start + 1 : outputs_end]:
        match = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        key = match.group(1) if match is not None else None
        if key in file_keys:
            continue
        if key == "console":
            if not console_written:
                rewritten.append("  console = false\n")
                console_written = True
            continue
        rewritten.append(line)
    if not console_written:
        rewritten.append("  console = false\n")
    if output_exodus is not None:
        rewritten.append(f"  exodus = {output_exodus}\n")
        rewritten.append("  exodus_interval = 1\n")
    rewritten.extend(lines[outputs_end:])
    return "".join(rewritten)


def compare_exodus_fields(
    fuelsim_path: Path,
    moose_path: Path,
    fuelsim_contact: Optional[str],
    moose_contact: Optional[str],
) -> Dict[str, object]:
    fuelsim = nc.Dataset(fuelsim_path)
    moose = nc.Dataset(moose_path)
    try:
        f_map = nodal_name_map(fuelsim)
        m_map = nodal_name_map(moose)

        f_time = float(fuelsim.variables["time_whole"][:][-1])
        m_time = float(moose.variables["time_whole"][:][-1])
        f_step = int(fuelsim.variables["vals_nod_var1"].shape[0] - 1)
        m_step = int(moose.variables["vals_nod_var1"].shape[0] - 1)

        f_x = np.array(fuelsim.variables["coordx"][:], dtype=float)
        f_y = np.array(fuelsim.variables["coordy"][:], dtype=float)
        m_x = np.array(moose.variables["coordx"][:], dtype=float)
        m_y = np.array(moose.variables["coordy"][:], dtype=float)

        if f_x.shape != m_x.shape:
            raise RuntimeError(
                f"Fuelsim and MOOSE node count differ: {f_x.shape[0]} vs {m_x.shape[0]}"
            )

        coord_max_err = {
            "max_abs_x_diff": float(np.max(np.abs(f_x - m_x))),
            "max_abs_y_diff": float(np.max(np.abs(f_y - m_y))),
        }

        result: Dict[str, object] = {
            "time": {"fuelsim": f_time, "moose": m_time, "time_delta": f_time - m_time},
            "steps": {"fuelsim": f_step + 1, "moose": m_step + 1},
            "coord_diff": coord_max_err,
            "fields": {},
        }

        def compare_pair(name_f: str, name_m: str, out_label: str) -> None:
            data_f = nodal_field(fuelsim, f_map, name_f, f_step)
            data_m = nodal_field(moose, m_map, name_m, m_step)
            result["fields"][out_label] = compute_metrics(data_f, data_m)

        compare_pair("temperature", "T", "temperature")
        compare_pair("displacement_r", "disp_x", "displacement_r")
        compare_pair("displacement_z", "disp_y", "displacement_z")

        # Contact pressure fields are optional but useful for friction validation.
        f_contact_name = select_contact_field(f_map, fuelsim_contact, "fuelsim")
        m_contact_name = select_contact_field(m_map, moose_contact, "moose")

        if f_contact_name is not None and m_contact_name is not None:
            compare_pair(f_contact_name, m_contact_name, "contact_pressure")
        else:
            result["fields"]["contact_pressure"] = {
                "skipped": True,
                "fuelsim_field": f_contact_name,
                "moose_field": m_contact_name,
                "reason": "missing_contact_pressure_field",
            }

        return result
    finally:
        fuelsim.close()
        moose.close()


def compare_equivalence(
    reference_path: Path,
    candidate_path: Path,
    fields: List[Tuple[str, str]],
) -> Dict[str, object]:
    ref = nc.Dataset(reference_path)
    cur = nc.Dataset(candidate_path)
    try:
        ref_map = nodal_name_map(ref)
        cur_map = nodal_name_map(cur)
        if len(ref.variables["coordx"][:]) != len(cur.variables["coordx"][:]):
            raise RuntimeError("Nodal count mismatch for equivalence comparison")

        ref_step = int(ref.variables["vals_nod_var1"].shape[0] - 1)
        cur_step = int(cur.variables["vals_nod_var1"].shape[0] - 1)

        out = {}
        for f_ref, f_cur in fields:
            a = nodal_field(cur, cur_map, f_cur, cur_step)
            b = nodal_field(ref, ref_map, f_ref, ref_step)
            out[f"{f_cur}_vs_{f_ref}"] = compute_metrics(a, b)
        return out
    finally:
        ref.close()
        cur.close()


def prepare_fuelsim_case(
    original: Path,
    output_exodus: Optional[Path],
    workspace: Path,
    filename: str,
) -> Path:
    workspace.mkdir(parents=True, exist_ok=True)
    prepared = workspace / filename
    text = original.read_text()
    text = replace_fuelsim_mesh_with_absolute_path(text, original)
    replaced = configure_fuelsim_outputs(text, output_exodus)
    prepared.write_text(replaced)
    return prepared


def run_command(cmd: List[str], env: Dict[str, str], log_path: Path) -> Dict[str, object]:
    start = time.perf_counter()
    with log_path.open("w") as log_file:
        proc = subprocess.run(
            cmd,
            env=env,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
    wall = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed, exit={proc.returncode}, log={log_path}")
    return {
        "command": shlex.join(cmd),
        "return_code": proc.returncode,
        "wall_seconds": wall,
        "log": str(log_path),
    }


def make_base_env(thread_count: int, extra_env: Dict[str, str]) -> Dict[str, str]:
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": str(thread_count),
        "OPENBLAS_NUM_THREADS": str(thread_count),
        "MKL_NUM_THREADS": str(thread_count),
        "NUMEXPR_NUM_THREADS": str(thread_count),
        "MPIR_CVAR_CH4_NETMOD": os.environ.get("MPIR_CVAR_CH4_NETMOD", "ofi"),
        "FI_PROVIDER": os.environ.get("FI_PROVIDER", "tcp"),
    })
    env.update(extra_env)
    return env


def run_rank_suite(
    rank: int,
    fuelsim_bin: Path,
    moose_bin: Path,
    mpiexec: str,
    fuelsim_case: Path,
    moose_input: Path,
    workspace: Path,
    thread_count: int,
    extra_env: Dict[str, str],
    extra_opts: List[str],
    moose_args: List[str],
    timing_repeats: int,
    run_fuelsim: bool,
    run_moose: bool,
    run_timing: bool,
) -> Dict[str, object]:
    out_dir = workspace / f"r{rank}"
    out_dir.mkdir(parents=True, exist_ok=True)
    entry: Dict[str, object] = {"rank": rank, "workspace": str(out_dir)}

    fuelsim_exodus = out_dir / f"fuelsim_r{rank}.e"
    moose_exodus = out_dir / f"mooseF_{rank}.e"
    env = make_base_env(thread_count, extra_env)
    mpiexec_cmd = [mpiexec, *extra_opts, "-n", str(rank)]

    if run_fuelsim:
        prepared_case = prepare_fuelsim_case(
            fuelsim_case, fuelsim_exodus, out_dir, "case_validation.fsi"
        )
        fuelsim_cmd = mpiexec_cmd + [str(fuelsim_bin), "-i", str(prepared_case)]
        entry["fuelsim"] = {
            "command": shlex.join(fuelsim_cmd),
            "run": run_command(
                fuelsim_cmd, env=env, log_path=out_dir / "fuelsim.log"
            ),
            "exodus": str(fuelsim_exodus),
            "timing_runs": [],
        }
    else:
        entry["fuelsim"] = {
            "run": {
                "wall_seconds": None,
                "command": "skip",
                "return_code": 0,
                "log": str(out_dir / "fuelsim.log"),
            },
            "exodus": str(fuelsim_exodus),
            "timing_runs": [],
        }

    if run_moose:
        moose_cmd = mpiexec_cmd + [
            str(moose_bin),
            "-i",
            str(moose_input),
            *moose_args,
            f"Outputs/file_base={out_dir / f'mooseF_{rank}'}",
            "Outputs/csv=true",
            "Outputs/exodus=true",
            "Outputs/console=false",
        ]
        entry["moose"] = {
            "command": shlex.join(moose_cmd),
            "run": run_command(
                moose_cmd, env=env, log_path=out_dir / "moose.log"
            ),
            "exodus": str(moose_exodus),
            "timing_runs": [],
        }
    else:
        entry["moose"] = {
            "run": {
                "wall_seconds": None,
                "command": "skip",
                "return_code": 0,
                "log": str(out_dir / "moose.log"),
            },
            "exodus": str(moose_exodus),
            "timing_runs": [],
        }

    if timing_repeats > 0 and run_timing:
        timing_case = None
        timing_case = prepare_fuelsim_case(
            fuelsim_case, None, out_dir, "case_timing.fsi"
        )
        for trial in range(1, timing_repeats + 1):
            fuelsim_timing_cmd = mpiexec_cmd + [
                str(fuelsim_bin), "-i", str(timing_case)
            ]
            entry["fuelsim"]["timing_runs"].append(
                run_command(
                    fuelsim_timing_cmd,
                    env=env,
                    log_path=out_dir / f"fuelsim_timing_{trial}.log",
                )
            )
            moose_timing_cmd = mpiexec_cmd + [
                str(moose_bin),
                "-i",
                str(moose_input),
                *moose_args,
                "Outputs/csv=false",
                "Outputs/exodus=false",
                "Outputs/console=false",
            ]
            entry["moose"]["timing_runs"].append(
                run_command(
                    moose_timing_cmd,
                    env=env,
                    log_path=out_dir / f"moose_timing_{trial}.log",
                )
            )

    for solver in ("fuelsim", "moose"):
        seconds = [
            run["wall_seconds"] for run in entry[solver]["timing_runs"]
        ]
        subsequent = seconds[1:]
        representative = statistics.median(subsequent or seconds) if seconds else None
        entry[solver]["timing_summary"] = {
            "first_seconds": seconds[0] if seconds else None,
            "subsequent_seconds": subsequent,
            "representative_seconds": representative,
            "representative_basis": (
                "median_of_subsequent_trials" if subsequent else
                "first_trial_only" if seconds else "not_measured"
            ),
            "file_output_disabled": True if seconds else None,
        }

    return entry


def fmt_float(value: float) -> str:
    if value is None or isinstance(value, str):
        return str(value)
    if isinstance(value, (float, int)) and np.isnan(value):
        return "nan"
    if isinstance(value, (float, int)) and np.isinf(value):
        return "inf"
    return f"{value:.6g}"


def print_metrics_table(results: Dict[str, object]) -> None:
    print("\n# Timing summary (file output disabled)")
    for rank in sorted(int(k) for k in results["runs"].keys()):
        entry = results["runs"][str(rank)]
        fs_summary = entry["fuelsim"]["timing_summary"]
        ms_summary = entry["moose"]["timing_summary"]
        print(
            f"r={rank:>2d} "
            f"fuelsim_first={fmt_float(fs_summary['first_seconds'])} s, "
            f"fuelsim_representative={fmt_float(fs_summary['representative_seconds'])} s, "
            f"moose_first={fmt_float(ms_summary['first_seconds'])} s, "
            f"moose_representative={fmt_float(ms_summary['representative_seconds'])} s"
        )

    baseline_rank = sorted(results["runs"].keys(), key=lambda item: int(item))[0]
    baseline = results["runs"][baseline_rank]
    b_fs = baseline["fuelsim"]["timing_summary"]["representative_seconds"]
    b_ms = baseline["moose"]["timing_summary"]["representative_seconds"]
    print(f"\n# Speedup summary (wall clock, compared to r={baseline_rank})")
    print("r\tFuelsim_speedup\tFuelsim_eff(%)\tMoose_speedup\tMoose_eff(%)")
    for rank in sorted(int(k) for k in results["runs"].keys()):
        entry = results["runs"][str(rank)]
        fs = entry["fuelsim"]["timing_summary"]["representative_seconds"]
        ms = entry["moose"]["timing_summary"]["representative_seconds"]
        fs_speed = b_fs / fs if b_fs is not None and fs not in (None, 0.0) else float("nan")
        ms_speed = b_ms / ms if b_ms is not None and ms not in (None, 0.0) else float("nan")
        relative_rank = rank / int(baseline_rank)
        print(
            f"{rank}\t{fmt_float(fs_speed)}\t\t{fmt_float(fs_speed / relative_rank * 100)}\t\t"
            f"{fmt_float(ms_speed)}\t\t{fmt_float(ms_speed / relative_rank * 100)}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run fuelsim + MOOSE for multiple MPI ranks, compare final Exodus"
            " fields, and print speedup/consistency metrics."
        )
    )
    parser.add_argument(
        "--fuelsim-case",
        required=True,
        type=Path,
        help="fuelsim input card path",
    )
    parser.add_argument(
        "--moose-input",
        required=True,
        type=Path,
        help="MOOSE input file path",
    )
    parser.add_argument(
        "--fuelsim-bin",
        type=Path,
        default=Path("build/fuelsim"),
        help="fuelsim executable",
    )
    parser.add_argument(
        "--moose-bin",
        type=Path,
        default=Path("/home/cooper/projects/july/july-opt"),
        help="MOOSE executable",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/pcmi_parallel_suite"),
        help="Workspace for logs and Exodus outputs",
    )
    parser.add_argument(
        "--ranks",
        nargs="+",
        type=int,
        default=DEFAULT_RANKS,
        help="MPI rank list",
    )
    parser.add_argument(
        "--thread-count",
        type=int,
        default=1,
        help="OpenMP/OpenBLAS/MKL/NumExpr thread count",
    )
    parser.add_argument(
        "--timing-repeats",
        type=int,
        default=2,
        help=(
            "Number of separate file-output-disabled timing trials per rank; "
            "the first and subsequent trials are reported separately"
        ),
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Only compare existing outputs in workspace/rN",
    )
    parser.add_argument(
        "--reuse-validation-results",
        action="store_true",
        help=(
            "Reuse existing fuelsim_rN.e and mooseF_N.e files for field "
            "comparison while still running the requested timing trials"
        ),
    )
    parser.add_argument(
        "--mpiexec",
        type=str,
        default="mpiexec",
        help="MPI launcher executable",
    )
    parser.add_argument(
        "--fuelsim-contact",
        default=None,
        help="Fuelsim contact pressure field name in nodal Exodus",
    )
    parser.add_argument(
        "--moose-contact",
        default=None,
        help="MOOSE contact pressure field name in nodal Exodus",
    )
    parser.add_argument(
        "--extra-env",
        action="append",
        default=[],
        help="Extra environment variable e.g. KEY=VALUE",
    )
    parser.add_argument(
        "--mpiexec-options",
        nargs="*",
        default=[],
        help="Extra options passed to mpiexec before -n (quoted as needed)",
    )
    parser.add_argument(
        "--moose-arg",
        action="append",
        default=[],
        help="Additional MOOSE command-line override; may be repeated",
    )
    parser.add_argument(
        "--no-cross-rank",
        action="store_true",
        help="Skip fuelsim/moose same-case cross-rank comparisons",
    )
    args = parser.parse_args()

    extra_env: Dict[str, str] = {}
    for item in args.extra_env:
        if "=" not in item:
            raise ValueError(f"Invalid --extra-env entry: {item}")
        key, value = item.split("=", 1)
        extra_env[key] = value

    if not args.ranks or any(rank <= 0 for rank in args.ranks):
        raise ValueError("--ranks must contain positive MPI rank counts")
    if args.thread_count <= 0:
        raise ValueError("--thread-count must be positive")
    if args.timing_repeats < 0:
        raise ValueError("--timing-repeats cannot be negative")
    if args.skip_run and args.reuse_validation_results:
        raise ValueError(
            "--skip-run and --reuse-validation-results are mutually exclusive"
        )
    if not args.fuelsim_case.exists():
        raise FileNotFoundError(f"fuelsim case not found: {args.fuelsim_case}")
    if not args.moose_input.exists():
        raise FileNotFoundError(f"moose input not found: {args.moose_input}")
    if not args.skip_run:
        if not args.fuelsim_bin.exists():
            raise FileNotFoundError(f"fuelsim executable not found: {args.fuelsim_bin}")
        if not args.moose_bin.exists():
            raise FileNotFoundError(f"moose executable not found: {args.moose_bin}")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Keep runs in a stable order.
    ranks = sorted(set(args.ranks))
    results = {
        "settings": {
            "fuelsim_case": str(args.fuelsim_case),
            "moose_input": str(args.moose_input),
            "fuelsim_bin": str(args.fuelsim_bin),
            "moose_bin": str(args.moose_bin),
            "ranks": ranks,
            "thread_count": args.thread_count,
            "timing_repeats": args.timing_repeats,
            "skip_run": args.skip_run,
            "reuse_validation_results": args.reuse_validation_results,
            "no_cross_rank": args.no_cross_rank,
        },
        "runs": {},
    }

    for rank in ranks:
        print(f"\n==== Running / comparing rank {rank} ====")
        entry = run_rank_suite(
            rank=rank,
            fuelsim_bin=args.fuelsim_bin,
            moose_bin=args.moose_bin,
            fuelsim_case=args.fuelsim_case,
            moose_input=args.moose_input,
            workspace=output_dir,
            thread_count=args.thread_count,
            extra_env=extra_env,
            extra_opts=args.mpiexec_options,
            moose_args=args.moose_arg,
            timing_repeats=0 if args.skip_run else args.timing_repeats,
            mpiexec=args.mpiexec,
            run_fuelsim=not args.skip_run and not args.reuse_validation_results,
            run_moose=not args.skip_run and not args.reuse_validation_results,
            run_timing=not args.skip_run,
        )
        results["runs"][str(rank)] = entry

        if (not args.skip_run) and rank in [1]:
            print(f"saved to {entry['fuelsim']['run']['log']}, {entry['moose']['run']['log']}")

        if not Path(entry["fuelsim"]["exodus"]).exists():
            raise FileNotFoundError(
                f"fuelsim output missing: {entry['fuelsim']['exodus']}"
            )
        if not Path(entry["moose"]["exodus"]).exists():
            raise FileNotFoundError(
                f"MOOSE output missing: {entry['moose']['exodus']}"
            )

        entry["comparison"] = compare_exodus_fields(
            fuelsim_path=Path(entry["fuelsim"]["exodus"]),
            moose_path=Path(entry["moose"]["exodus"]),
            fuelsim_contact=args.fuelsim_contact,
            moose_contact=args.moose_contact,
        )

    if len(ranks) >= 2 and not args.no_cross_rank:
        base_rank = str(ranks[0])
        base_fu = Path(results["runs"][base_rank]["fuelsim"]["exodus"])
        base_mo = Path(results["runs"][base_rank]["moose"]["exodus"])

        for rank in ranks[1:]:
            entry = results["runs"][str(rank)]
            candidate_fu = Path(entry["fuelsim"]["exodus"])
            candidate_mo = Path(entry["moose"]["exodus"])
            fuelsim_fields = [
                ("temperature", "temperature"),
                ("displacement_r", "displacement_r"),
                ("displacement_z", "displacement_z"),
            ]
            moose_fields = [
                ("T", "T"),
                ("disp_x", "disp_x"),
                ("disp_y", "disp_y"),
            ]
            base_fu_contact = exodus_contact_field(
                base_fu, args.fuelsim_contact, "fuelsim"
            )
            candidate_fu_contact = exodus_contact_field(
                candidate_fu, args.fuelsim_contact, "fuelsim"
            )
            if base_fu_contact is not None and candidate_fu_contact is not None:
                fuelsim_fields.append((base_fu_contact, candidate_fu_contact))
            base_mo_contact = exodus_contact_field(
                base_mo, args.moose_contact, "moose"
            )
            candidate_mo_contact = exodus_contact_field(
                candidate_mo, args.moose_contact, "moose"
            )
            if base_mo_contact is not None and candidate_mo_contact is not None:
                moose_fields.append((base_mo_contact, candidate_mo_contact))
            fuelsim_equiv = compare_equivalence(
                reference_path=base_fu,
                candidate_path=candidate_fu,
                fields=fuelsim_fields,
            )
            moose_equiv = compare_equivalence(
                reference_path=base_mo,
                candidate_path=candidate_mo,
                fields=moose_fields,
            )
            entry["equivalence"] = {
                "fuelsim_vs_rank1": fuelsim_equiv,
                "moose_vs_rank1": moose_equiv,
            }

    print_metrics_table(results)

    report = output_dir / "pcmi_parallel_report.json"
    report.write_text(json.dumps(results, indent=2))
    print(f"\nreport written: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
