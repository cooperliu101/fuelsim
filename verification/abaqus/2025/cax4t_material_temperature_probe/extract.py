"""Run with Abaqus Python. Export native material points without extrapolation."""
from __future__ import print_function

import csv
import sys
from odbAccess import openOdb
from abaqusConstants import INTEGRATION_POINT, NODAL


def payload(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field(frame, name, nodal=False, required=True):
    # Test the Python key list: older Abaqus repository membership checks can
    # leave a pending KeyError for legitimately absent optional output fields.
    if name not in frame.fieldOutputs.keys():
        if required:
            raise RuntimeError("Missing required field " + name)
        return {}
    output = frame.fieldOutputs[name]
    result = {}
    for value in output.values:
        expected = NODAL if nodal else INTEGRATION_POINT
        # COORD can contain both native nodal and material-point entries.
        if name == "COORD" and value.position == NODAL and not nodal:
            continue
        if value.position != expected:
            raise RuntimeError("Unexpected output position for " + name)
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate field entry for " + name)
        result[key] = payload(value)
    if name in ("S", "EE", "PE", "CE"):
        expected = tuple(name + suffix for suffix in ("11", "22", "33", "12"))
        if tuple(output.componentLabels) != expected:
            raise RuntimeError("Unexpected component ordering for " + name)
    return result


def tensor(values, key, strain=False):
    if key not in values:
        return [float("nan")] * 4
    result = list(values[key])
    if strain:
        result[3] *= 0.5
    return result


def output_file(path):
    return open(path, "w", newline="") if sys.version_info[0] < 3 else open(path, "w", newline="")


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract.py job.odb output_prefix")
odb = openOdb(sys.argv[1], readOnly=True)
prefix = sys.argv[2]
streams = [output_file(prefix + suffix) for suffix in ("_nodes.csv", "_points.csv")]
nodes, points = [csv.writer(stream, lineterminator="\n") for stream in streams]
nodes.writerow(["step", "frame", "time", "node", "r", "z", "temperature", "ur", "uz"])
points.writerow(["step", "frame", "time", "element", "point", "r", "z"]
                + [p + c for p in ("s_", "ee_", "pe_") for c in ("rr", "zz", "hoop", "rz")]
                + ["peeq"] + ["ce_" + c for c in ("rr", "zz", "hoop", "rz")] + ["ceeq"])
accepted_frames = 0
try:
    coordinates = {}
    for instance in odb.rootAssembly.instances.values():
        for node in instance.nodes:
            if node.label in coordinates:
                raise RuntimeError("Duplicate source node")
            coordinates[node.label] = node.coordinates
    if set(coordinates) != set((1, 2, 3, 4)):
        raise RuntimeError("Expected exactly four nodes")
    elapsed = 0.0
    for step_name, step in odb.steps.items():
        for frame_index, frame in enumerate(step.frames):
            temperatures = field(frame, "NT11", True)
            displacements = field(frame, "U", True)
            if set(temperatures) != set(coordinates) or set(displacements) != set(coordinates):
                raise RuntimeError("Incomplete nodal coverage")
            time = elapsed + frame.frameValue
            for node in sorted(coordinates):
                nodes.writerow([step_name, frame_index, time, node] + list(coordinates[node][:2])
                               + [temperatures[node]] + list(displacements[node][:2]))
            required = ["S", "EE"]
            if "plastic_" in prefix:
                required += ["PE", "PEEQ"]
            if "creep_" in prefix:
                required += ["CE", "CEEQ"]
            fields = dict((name, field(frame, name, required=name in required))
                          for name in ("S", "EE", "PE", "PEEQ", "CE", "CEEQ", "COORD"))
            expected = set((1, q) for q in (1, 2, 3, 4))
            for name in required:
                if set(fields[name]) != expected:
                    raise RuntimeError("Incomplete integration-point coverage in " + name)
            for key in sorted(expected):
                coord = list(fields["COORD"][key][:2]) if key in fields["COORD"] else [float("nan")] * 2
                points.writerow([step_name, frame_index, time, key[0], key[1]] + coord
                                + tensor(fields["S"], key) + tensor(fields["EE"], key, True)
                                + tensor(fields["PE"], key, True) + [fields["PEEQ"].get(key, float("nan"))]
                                + tensor(fields["CE"], key, True) + [fields["CEEQ"].get(key, float("nan"))])
            if frame_index:
                accepted_frames += 1
        elapsed += step.timePeriod
finally:
    for stream in streams:
        stream.close()
    odb.close()
with open(prefix + "_export_complete.txt", "w") as stream:
    stream.write("accepted_frames=%d\n" % accepted_frames)
    stream.write("position=INTEGRATION_POINT\ncomponent_order=11,22,33,12\nstrain_shear=tensor\n")
print("exported_accepted_frames=%d" % accepted_frames)
