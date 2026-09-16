#!/usr/bin/env python3
"""Check this case's absolute-tolerance boundaries without running either solver."""

import numpy as np

from compare import metrics


def check(label, name, actual, reference, tolerance, strict, qualified, absolute_accepted):
    row = metrics(name, actual, reference, [str(i) for i in range(len(reference))], tolerance, "check")
    expected = {"strict_0_1_percent_passed": strict, "qualified_passed": qualified,
                "nonzero_pointwise_absolute_accepted": absolute_accepted}
    for key, value in expected.items():
        if row[key] != value:
            raise RuntimeError("%s: %s was %s, expected %s" % (label, key, row[key], value))
    print("PASS " + label)


def main():
    for name, tolerance in (("displacement", 1e-12), ("reaction", 1e-5), ("contact_slip", 1e-12)):
        for factor, accepted in ((0.5, 1), (1.0, 0), (2.0, 0)):
            # A unit reference keeps the aggregate checks below their thresholds.
            reference = [[factor * tolerance, 0, 0], [1, 0, 0]]
            check("%s nonzero absolute boundary %g" % (name, factor), name,
                  [[0, 0, 0], [1, 0, 0]], reference, tolerance, 0, accepted, accepted)
            check("%s zero reference boundary %g" % (name, factor), name,
                  reference, [[0, 0, 0], [1, 0, 0]], tolerance, accepted, accepted, 0)
        check(name + " relative pass with larger absolute difference", name,
              [[1.0001, 0, 0]], [[1, 0, 0]], tolerance, 1, 1, 0)

    check("unapproved field retains its relative check", "fuel_stress",
          [[0, 0, 0], [1, 0, 0]], [[1e-20, 0, 0], [1, 0, 0]], 1e-4, 0, 0, 0)
    check("absolute check uses the vector difference norm", "displacement",
          [[0.8e-12, 0, 0], [1, 0, 0]], [[0, 0.8e-12, 0], [1, 0, 0]], 1e-12, 0, 0, 0)
    check("relative threshold equality still fails", "displacement",
          [[1001, 0, 0], [1e6, 0, 0]], [[1000, 0, 0], [1e6, 0, 0]], 1e-12, 0, 0, 0)
    check("entire near-zero field uses its absolute check", "contact_slip",
          [[0, 0, 0]], [[1e-19, 0, 0]], 1e-12, 0, 1, 1)
    check("entire near-zero field rejects excessive absolute difference", "contact_slip",
          [[2e-12, 0, 0]], [[1e-19, 0, 0]], 1e-12, 0, 0, 0)
    check("field scale at the absolute threshold retains aggregate relative checks", "contact_slip",
          [[0.5e-12, 0, 0]], [[1e-12, 0, 0]], 1e-12, 0, 0, 1)
    check("unapproved entire near-zero field retains its relative check", "fuel_stress",
          [[0, 0, 0]], [[1e-20, 0, 0]], 1e-4, 0, 0, 0)

    reference = np.tile([1e-10, 0, 0], (100, 1))
    actual = reference.copy()
    actual[0, 1] = 5e-13
    check("aggregate peak failure cannot pass by absolute tolerance", "displacement",
          actual, reference, 1e-12, 0, 0, 1)
    reference = np.tile([1e-30, 0, 0], (100, 1))
    reference[0, 0] = 1e-10
    actual = reference.copy()
    actual[1:, 1] = 5e-14
    check("aggregate L2 failure cannot pass by absolute tolerance", "displacement",
          actual, reference, 1e-12, 0, 0, 99)

    print("All 30 tolerance checks passed.")


if __name__ == "__main__":
    main()
