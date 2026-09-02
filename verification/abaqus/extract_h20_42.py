from __future__ import print_function

import csv
import sys

from odbAccess import openOdb


if len(sys.argv) != 4:
    raise RuntimeError("usage: extract_h20_42.py <job.odb> <manifest.csv> <output.csv>")


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


states = {}
manifest = open(sys.argv[2], "rb")
reader = csv.DictReader(manifest)
for row in reader:
    states[int(row["copy"])] = row["state"]
manifest.close()

contact_nodes = set([
    71, 74, 75, 78, 82, 83, 86, 90, 92, 94, 97, 99, 102, 103, 106,
    107, 110, 114, 116, 118, 121, 122, 125, 126, 129, 133, 135, 137, 140,
])
odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    frame = odb.steps["LOAD"].frames[-1]
    by_prefix = {}
    for prefix in ("COPEN", "CPRESS"):
        values = {}
        matches = [frame.fieldOutputs[name] for name in frame.fieldOutputs.keys()
                   if name.strip().startswith(prefix)]
        if not matches:
            raise RuntimeError("field %s is missing" % prefix)
        for field in matches:
            for value in field.values:
                original_label = value.nodeLabel % 1000
                if original_label not in contact_nodes:
                    continue
                copy_index = value.nodeLabel // 1000
                key = (copy_index, original_label)
                if key in values:
                    raise RuntimeError("field %s contains duplicate node %s" % (prefix, key))
                values[key] = data(value)
        by_prefix[prefix] = values

    output = open(sys.argv[3], "wb")
    writer = csv.writer(output)
    writer.writerow(["state", "label", "copen", "cpress"])
    for copy_index in sorted(states):
        for label in sorted(contact_nodes):
            key = (copy_index, label)
            writer.writerow([states[copy_index], label, by_prefix["COPEN"][key], by_prefix["CPRESS"][key]])
    output.close()
finally:
    odb.close()
