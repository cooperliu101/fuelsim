from __future__ import print_function

import sys

from odbAccess import openOdb


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b44_swapped_operator.py <job.odb> <operator.csv>")

contact_labels = [7, 11, 15, 23, 25, 27, 33, 35, 37]
steps = [("BASE", "base", 0, 0.0)]
for index in range(1, 7):
    steps.append(("T%d_PLUS" % index, "tile", index, 1.0e-6))
    steps.append(("T%d_MINUS" % index, "tile", index, -1.0e-6))
for index in range(1, 10):
    steps.append(("S%d_PLUS" % index, "secondary", index, 1.0e-6))
    steps.append(("S%d_MINUS" % index, "secondary", index, -1.0e-6))


def field_with_prefix(frame, prefix):
    matches = [frame.fieldOutputs[name] for name in frame.fieldOutputs.keys() if name.strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError("field %s has %d matches, expected one" % (prefix, len(matches)))
    return matches[0]


def nodal_values(frame, prefix):
    result = {}
    for value in field_with_prefix(frame, prefix).values:
        if value.nodeLabel in contact_labels:
            result[value.nodeLabel] = value.dataDouble
    if len(result) != len(contact_labels):
        raise RuntimeError("%s has %d values, expected %d" % (prefix, len(result), len(contact_labels)))
    return result


def component(value, index):
    return value[index] if hasattr(value, "__len__") and index < len(value) else 0.0


odb = openOdb(path=sys.argv[1], readOnly=True)
output = open(sys.argv[2], "wb")
output.write(
    "step,input_kind,input_index,delta_m,output_index,output_label,copen_m,cpress_pa,"
    "cnormf_x_n,cnormf_y_n,cnormf_z_n\n"
)
for step_name, kind, input_index, delta in steps:
    frames = odb.steps[step_name].frames
    if len(frames) != 2:
        raise RuntimeError("step %s has %d frames, expected two" % (step_name, len(frames)))
    frame = frames[-1]
    opening = nodal_values(frame, "COPEN")
    pressure = nodal_values(frame, "CPRESS")
    force = nodal_values(frame, "CNORMF")
    for output_index, label in enumerate(contact_labels, 1):
        value = force[label]
        output.write(
            "%s,%s,%d,%.16g,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g\n"
            % (
                step_name,
                kind,
                input_index,
                delta,
                output_index,
                label,
                opening[label],
                pressure[label],
                component(value, 0),
                component(value, 1),
                component(value, 2),
            )
        )
output.close()
odb.close()
