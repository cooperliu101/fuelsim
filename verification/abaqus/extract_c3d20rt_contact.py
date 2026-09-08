from __future__ import print_function
import csv
import sys
from odbAccess import openOdb

def data(v):
    try:
        return v.dataDouble
    except Exception:
        return v.data

def nodal(frame, prefix):
    keys = [k for k in frame.fieldOutputs.keys() if k.strip() == prefix or k.strip().startswith(prefix + " ")]
    if len(keys) != 1:
        raise RuntimeError("Expected one field for %s, got %s" % (prefix, keys))
    return dict((v.nodeLabel, data(v)) for v in frame.fieldOutputs[keys[0]].values)

secondary = [int(value) for value in sys.argv[3].split(",")]
odb = openOdb(sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2] + "_nodes.csv", "wb") as nf, open(sys.argv[2] + "_contact.csv", "wb") as cf, \
         open(sys.argv[2] + "_points.csv", "wb") as pf:
        nodes, contact, points = csv.writer(nf), csv.writer(cf), csv.writer(pf)
        nodes.writerow(["time", "node", "temperature", "ux", "uy", "uz", "rx", "ry", "rz", "heat_reaction"])
        contact.writerow(["time", "node", "gap", "pressure", "nx", "ny", "nz", "tx", "ty", "tz",
                          "slipx", "slipy", "slipz", "heat_rate_w"])
        points.writerow(["time", "element", "point", "sxx", "syy", "szz", "sxy", "syz", "sxz"])
        elapsed = 0.0
        for step in odb.steps.values():
            for index in range(1, len(step.frames)):
                frame = step.frames[index]
                time = elapsed + frame.frameValue
                fields = dict((k, nodal(frame, k)) for k in ["NT11", "U", "RF", "RFL11", "COPEN", "CPRESS",
                                                              "CNORMF", "CSHEARF", "CSLIP1", "CSLIP2", "HFLA", "CTANDIR1", "CTANDIR2"])
                for node in sorted(fields["U"]):
                    nodes.writerow([time, node, fields["NT11"][node]] + list(fields["U"][node]) +
                                   list(fields["RF"][node]) + [fields["RFL11"].get(node, 0.0)])
                for node in secondary:
                    slip = [fields["CSLIP1"][node] * fields["CTANDIR1"][node][i] +
                            fields["CSLIP2"][node] * fields["CTANDIR2"][node][i] for i in range(3)]
                    contact.writerow([time, node, fields["COPEN"][node], fields["CPRESS"][node]] +
                                     list(fields["CNORMF"][node]) + list(fields["CSHEARF"][node]) +
                                     slip + [fields["HFLA"].get(node, 0.0)])
                for value in frame.fieldOutputs["S"].values:
                    tensor = data(value)
                    points.writerow([time, value.elementLabel, value.integrationPoint] +
                                    [tensor[i] for i in [0, 1, 2, 3, 5, 4]])
            elapsed += step.timePeriod
finally:
    odb.close()
