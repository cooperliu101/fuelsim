from __future__ import print_function
import csv
import sys
from odbAccess import openOdb


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field(frame, name, nodal=True):
    result = {}
    for value in frame.fieldOutputs[name].values:
        key = value.nodeLabel if nodal else (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate output: " + name)
        result[key] = data(value)
    return result


if len(sys.argv) != 3:
    raise RuntimeError("usage: extract_b8_rz.py job.odb prefix")
odb = openOdb(sys.argv[1], readOnly=True)
try:
    with open(sys.argv[2]+"_nodes.csv", "wb") as nf, open(sys.argv[2]+"_points.csv", "wb") as pf,\
            open(sys.argv[2]+"_contact.csv", "wb") as cf:
        nodes, points, contact = [csv.writer(f, lineterminator="\n") for f in (nf, pf, cf)]
        nodes.writerow(["time", "node", "temperature", "ur", "uz", "rf_r", "rf_z", "reaction_heat"])
        components = ["rr", "zz", "hoop", "rz"]
        points.writerow(["time", "element", "point"]+
                        [p+c for p in ["stress_", "elastic_", "plastic_", "creep_"] for c in components]+
                        ["equiv_plastic", "equiv_creep"])
        contact.writerow(["time", "node", "pressure", "gap", "shear", "slip", "normal_r", "normal_z",
                          "tangential_r", "tangential_z"])
        elapsed = 0.0
        for step in odb.steps.values():
            for i in range(1, len(step.frames)):
                frame = step.frames[i]
                time = elapsed+frame.frameValue
                fields = dict((k, field(frame, k)) for k in ["NT11", "U", "RF", "RFL11"])
                for n in sorted(fields["U"]):
                    nodes.writerow([time, n, fields["NT11"][n]]+list(fields["U"][n][:2])+
                                   list(fields["RF"][n][:2])+[fields["RFL11"][n]])
                stress, elastic = field(frame, "S", False), field(frame, "EE", False)
                for key in sorted(stress):
                    strain = list(elastic[key])
                    strain[3] *= 0.5
                    points.writerow([time, key[0], key[1]]+list(stress[key])+strain+[0.0]*10)
                fields = {}
                for short in ["CPRESS", "COPEN", "CSHEAR1", "CSLIP1", "CNORMF", "CSHEARF"]:
                    keys = [k for k in frame.fieldOutputs.keys()
                            if k.split()[0] == short and "S_INNER/S_OUTER" in k]
                    if len(keys) != 1:
                        raise RuntimeError("Missing or ambiguous contact output " + short)
                    fields[short] = field(frame, keys[0])
                for n in sorted(fields["COPEN"]):
                    contact.writerow([time, n]+[fields[k][n] for k in ["CPRESS", "COPEN", "CSHEAR1", "CSLIP1"]]+
                                     list(fields["CNORMF"][n][:2])+list(fields["CSHEARF"][n][:2]))
            elapsed += step.timePeriod
finally:
    odb.close()
