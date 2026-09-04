import csv
import math
import sys


if len(sys.argv) != 7:
    raise RuntimeError(
        "usage: compare_b61_c3d20t_replay.py <replay.csv> <step1_integration.csv> "
        "<step2_integration.csv> <step3_integration.csv> <step4_integration.csv> <step5_integration.csv>"
    )


def read_rows(path, key_names):
    with open(path, newline="") as source:
        return {
            tuple(int(row[name]) for name in key_names): row
            for row in csv.DictReader(source)
        }


def metrics(calculated, reference):
    differences = [left - right for left, right in zip(calculated, reference)]
    relative_l2 = 100.0 * math.sqrt(
        sum(value * value for value in differences) / sum(value * value for value in reference)
    )
    calculated_peak = max(abs(value) for value in calculated)
    reference_peak = max(abs(value) for value in reference)
    relative_peak = 100.0 * abs(calculated_peak - reference_peak) / reference_peak
    maximum_pointwise = 100.0 * max(
        abs(difference) / abs(value)
        for difference, value in zip(differences, reference)
        if value != 0.0
    )
    return relative_l2, relative_peak, maximum_pointwise, max(abs(value) for value in differences)


replay = read_rows(sys.argv[1], ("step", "element", "integration_point"))
fields = ("vonmises_stress", "peeq", "ceeq")
worst_l2 = {field: 0.0 for field in fields}
for step in range(1, 6):
    reference = read_rows(sys.argv[step + 1], ("element", "integration_point"))
    keys = sorted(reference)
    for field in fields:
        calculated_values = [float(replay[(step,) + key][field]) for key in keys]
        reference_values = [float(reference[key][field]) for key in keys]
        values = metrics(calculated_values, reference_values)
        worst_l2[field] = max(worst_l2[field], values[0])
        prefix = "step_%d_%s" % (step, field)
        print(
            "%s_relative_l2_percent=%.12g\n"
            "%s_relative_peak_percent=%.12g\n"
            "%s_maximum_pointwise_percent=%.12g\n"
            "%s_maximum_absolute_difference=%.12g"
            % ((prefix, values[0], prefix, values[1], prefix, values[2], prefix, values[3]))
        )
for field in fields:
    print("five_step_%s_maximum_relative_l2_percent=%.12g" % (field, worst_l2[field]))
