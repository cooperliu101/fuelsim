"""Compare production Exodus arrays, excluding creation/provenance text records."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    reader = Path(__file__).resolve().parents[1] / "compare.py"
    spec = importlib.util.spec_from_file_location("rz_comparison", reader)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.result_database(args.baseline) as old, module.result_database(args.candidate) as new:
        if old.variables.keys() != new.variables.keys():
            raise ValueError("Result variable sets differ")
        numeric = [name for name, value in old.variables.items() if value.dtype.kind != "S"]
        text = [name for name, value in old.variables.items()
                if value.dtype.kind == "S" and name not in ("qa_records", "info_records")]
        differences = [name for name in numeric
                       if not np.array_equal(old.variables[name], new.variables[name], equal_nan=True)]
        differences += [name for name in text
                        if not np.array_equal(old.variables[name], new.variables[name])]
        report = {"frames": len(new.variables["time_whole"]), "numeric_arrays": len(numeric),
                  "text_arrays": len(text), "different_arrays": differences,
                  "all_identical": not differences}
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))
    return 0 if report["all_identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
