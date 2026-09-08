"""Separate contact-transfer differences from finite-overclosure effects."""
import numpy as np
from analyze_h20_28_finite_transfer import extract

base = extract("h20_28_finite_transfer_probe")
shallow = extract("h20_28_finite_transfer_shallow_probe")
small = extract("h20_28_hex20_refined_primary_sts_default")
np.testing.assert_array_equal(base[3], shallow[3])
np.testing.assert_array_equal(base[4], shallow[4])
print("base_overclosure_m=1e-4")
print("shallow_overclosure_m=1e-6")
for index, name in enumerate(["secondary_average", "area", "primary_transfer"]):
    difference = shallow[index] - base[index]
    print(name + "_depth_relative_difference=" +
          str(np.linalg.norm(difference) / np.linalg.norm(base[index])))
    print(name + "_depth_maximum_absolute_difference=" + str(abs(difference).max()))
print("shallow_vs_small_primary_transfer_relative=" +
      str(np.linalg.norm(shallow[2] - small[2]) / np.linalg.norm(small[2])))
midpoint = (np.rint(shallow[3] * 16).astype(int) % 2).sum(axis=1) == 1
print("shallow_corner_minimum_midpoint_weight=" + str(shallow[2][0, midpoint].min()))
print("base_secondary_average_condition_number=" + str(np.linalg.cond(base[0])))
print("shallow_secondary_average_condition_number=" + str(np.linalg.cond(shallow[0])))
