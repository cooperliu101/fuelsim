from __future__ import print_function

import csv
import gzip
import sys
from odbAccess import openOdb


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def values(frame, name, nodal=False):
    if name not in list(frame.fieldOutputs.keys()):
        return {}
    result = {}
    for value in frame.fieldOutputs[name].values:
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate output entry in " + name)
        result[key] = data(value)
    return result


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_results.py job.odb prefix")
odb = openOdb(sys.argv[1], readOnly=True)
try:
    with gzip.GzipFile(sys.argv[2]+"_nodes.csv.gz", "wb", mtime=0) as nf, gzip.GzipFile(sys.argv[2]+"_points.csv.gz", "wb", mtime=0) as pf, gzip.GzipFile(sys.argv[2]+"_contact.csv.gz", "wb", mtime=0) as cf:
        nodes, points, contact = [csv.writer(f) for f in (nf, pf, cf)]
        nodes.writerow(["time", "node", "temperature", "ur", "uz", "rf_r", "rf_z", "reaction_heat"])
        components = ["rr", "zz", "hoop", "rz"]
        points.writerow(["time", "element", "point"] + [p+c for p in ["stress_", "elastic_", "plastic_", "creep_"] for c in components] + ["equiv_plastic", "equiv_creep", "r", "z", "volume"])
        contact.writerow(["time", "node", "pressure", "gap", "shear", "slip", "normal_r", "normal_z", "tangential_r", "tangential_z", "heat"])
        elapsed = 0.0
        for step in odb.steps.values():
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = elapsed + frame.frameValue
                fields = dict((k, values(frame, k, True)) for k in ["NT11", "U", "RF", "RFL11"])
                for n in sorted(fields["U"]):
                    nodes.writerow([time, n, fields["NT11"].get(n, float("nan"))]+list(fields["U"][n][:2])+list(fields["RF"][n][:2])+[fields["RFL11"].get(n, float("nan"))])
                fields = dict((k, values(frame, k)) for k in ["S", "EE", "PE", "CE", "PEEQ", "CEEQ", "COORD", "IVOL"])
                for key in sorted(fields["S"]):
                    row = [time, key[0], key[1]]
                    for k in ["S", "EE", "PE", "CE"]:
                        v = list(fields[k].get(key, (0.0,)*4))
                        if k != "S": v[3] *= .5
                        row += v
                    row += [fields["PEEQ"].get(key, 0.0), fields["CEEQ"].get(key, 0.0)]
                    row += list(fields["COORD"][key][:2])+[fields["IVOL"][key]]
                    points.writerow(row)
                fields = {}
                for k in ["CPRESS", "COPEN", "CSHEAR1", "CSLIP1", "CNORMF", "CSHEARF", "HFL"]:
                    keys = [n for n in frame.fieldOutputs.keys() if n.split()[0] == k and "S_FUEL_RIGHT/S_CLAD_LEFT" in n]
                    if len(keys) > 1: raise RuntimeError("Ambiguous interface field "+k)
                    fields[k] = values(frame, keys[0], True) if keys else {}
                if not fields["COPEN"]: raise RuntimeError("Missing interface reference")
                for n in sorted(fields["COPEN"]):
                    contact.writerow([time, n]+[fields[k].get(n, 0.0) for k in ["CPRESS", "COPEN", "CSHEAR1", "CSLIP1"]]+list(fields["CNORMF"].get(n, (0.0, 0.0))[:2])+list(fields["CSHEARF"].get(n, (0.0, 0.0))[:2])+[fields["HFL"].get(n, 0.0)])
            elapsed += step.timePeriod
finally:
    odb.close()
