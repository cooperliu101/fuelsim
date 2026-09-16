"""Read every accepted frame from the independent Abaqus 2025 calculation."""

import csv
import gzip
import json
import sys

from odbAccess import openOdb


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def nodal(frame, prefix, labels, contact=False):
    names = [name for name in frame.fieldOutputs.keys()
             if name.strip() == prefix or (contact and name.strip().startswith(prefix + " "))]
    if len(names) != 1:
        raise RuntimeError("Expected one %s field, got %s" % (prefix, names))
    result = {}
    for value in frame.fieldOutputs[names[0]].values:
        if value.nodeLabel not in labels:
            continue
        if value.nodeLabel in result:
            raise RuntimeError("Duplicate %s node %d" % (prefix, value.nodeLabel))
        result[value.nodeLabel] = data(value)
    if set(result) != labels:
        raise RuntimeError("Incomplete %s nodal coverage: %d of %d" % (prefix, len(result), len(labels)))
    return result


def integration(frame, name, expected):
    result = {}
    for value in frame.fieldOutputs[name].values:
        if str(value.position) != "INTEGRATION_POINT":
            if name == "COORD" and str(value.position) == "NODAL":
                continue
            raise RuntimeError("Expected material-point field %s" % name)
        key = (value.elementLabel, value.integrationPoint)
        if key in result:
            raise RuntimeError("Duplicate %s point %s" % (name, key))
        result[key] = data(value)
    if set(result) != expected:
        raise RuntimeError("Incomplete %s material-point coverage: %d of %d" % (name, len(result), len(expected)))
    return result


def main():
    if len(sys.argv) != 4:
        raise RuntimeError("usage: abaqus python extract.py job.odb mesh.json output_prefix")
    with open(sys.argv[2]) as stream:
        mesh = json.load(stream)
    nodes = {node["label"] for node in mesh["nodes"]}
    elements = {element["label"]: element for element in mesh["elements"]}
    clad = {element["label"] for block in mesh["blocks"] if block["name"] == "CLAD"
            for element in block["elements"]}
    corners = {node for element in elements.values() for node in element["nodes"][:8]}
    face_nodes = ((0, 1, 5, 4, 8, 13, 16, 12), (1, 2, 6, 5, 9, 14, 17, 13),
                  (2, 3, 7, 6, 10, 15, 18, 14), (3, 0, 4, 7, 11, 12, 19, 15),
                  (0, 3, 2, 1, 11, 10, 9, 8), (4, 5, 6, 7, 16, 17, 18, 19))
    faces = next(side["faces"] for side in mesh["side_sets"] if side["name"] == "FUEL_OUTER")
    secondary = {elements[face["element"]]["nodes"][index] for face in faces
                 for index in face_nodes[face["exodus_side"] - 1]}
    if (len(nodes), len(corners), len(elements), len(clad), len(secondary)) != (5969, 1617, 1152, 512, 416):
        raise RuntimeError("Unexpected medium mesh counts")
    all_points = {(element, point) for element in elements for point in range(1, 9)}
    clad_points = {(element, point) for element in clad for point in range(1, 9)}
    prefix = sys.argv[3]
    odb = openOdb(sys.argv[1], readOnly=True)
    try:
        step = odb.steps["PATH"]
        if abs(step.frames[-1].frameValue - 1.0) > 1e-12:
            raise RuntimeError("Abaqus did not reach one second")
        with gzip.open(prefix + "_nodes.csv.gz", "wt", newline="") as nf, \
             gzip.open(prefix + "_contact.csv.gz", "wt", newline="") as cf, \
             gzip.open(prefix + "_points.csv.gz", "wt", newline="") as pf, \
             open(prefix + "_frames.csv", "w", newline="") as ff:
            writers = [csv.writer(stream, lineterminator="\n") for stream in (nf, cf, pf, ff)]
            nw, cw, pw, fw = writers
            nw.writerow(["time", "node", "temperature", "ux", "uy", "uz", "rx", "ry", "rz", "heat_reaction"])
            cw.writerow(["time", "node", "gap", "pressure", "nx", "ny", "nz", "tx", "ty", "tz",
                         "slipx", "slipy", "slipz", "heat_rate"])
            components = ("xx", "yy", "zz", "xy", "yz", "xz")
            pw.writerow(["time", "element", "point", "x", "y", "z", "ivol"] +
                        [prefix + component for prefix in ("s_", "e_", "p_", "c_") for component in components] +
                        ["equiv_plastic", "equiv_creep"])
            fw.writerow(["time", "increment", "nodes", "contact_nodes", "material_points"])
            for frame_index in range(1, len(step.frames)):
                frame = step.frames[frame_index]
                time = frame.frameValue
                fields = {name: nodal(frame, name, nodes) for name in ("NT11", "U", "RF")}
                thermal_reaction = nodal(frame, "RFL11", corners)
                for node in sorted(nodes):
                    nw.writerow([time, node, fields["NT11"][node]] + list(fields["U"][node]) +
                                list(fields["RF"][node]) + [thermal_reaction.get(node, 0.0)])
                fields = {name: nodal(frame, name, secondary, contact=True)
                          for name in ("COPEN", "CPRESS", "CNORMF", "CSHEARF", "CSLIP1", "CSLIP2",
                                       "CTANDIR1", "CTANDIR2")}
                # HFLA is the native nodal heat transfer, while HFL is a flux density.
                heat = nodal(frame, "HFLA", secondary.intersection(corners), contact=True)
                for node in sorted(secondary):
                    slip = [fields["CSLIP1"][node] * fields["CTANDIR1"][node][i] +
                            fields["CSLIP2"][node] * fields["CTANDIR2"][node][i] for i in range(3)]
                    cw.writerow([time, node, fields["COPEN"][node], fields["CPRESS"][node]] +
                                list(fields["CNORMF"][node]) + list(fields["CSHEARF"][node]) +
                                slip + [heat.get(node, 0.0)])
                fields = {name: integration(frame, name, all_points) for name in ("S", "EE", "COORD", "IVOL")}
                for name in ("PE", "CE", "PEEQ", "CEEQ"):
                    fields[name] = integration(frame, name, clad_points)
                for key in sorted(all_points):
                    row = [time, key[0], key[1]] + list(fields["COORD"][key]) + [fields["IVOL"][key]]
                    for name in ("S", "EE", "PE", "CE"):
                        vector = fields[name].get(key, (0.0,) * 6)
                        row += [vector[index] * (0.5 if name != "S" and index >= 3 else 1.0)
                                for index in (0, 1, 2, 3, 5, 4)]
                    row += [fields[name].get(key, 0.0) for name in ("PEEQ", "CEEQ")]
                    pw.writerow(row)
                fw.writerow([time, frame.incrementNumber, len(nodes), len(secondary), len(all_points)])
                print("time=%.17g nodes=%d contact_nodes=%d material_points=%d" %
                      (time, len(nodes), len(secondary), len(all_points)))
    finally:
        odb.close()


if __name__ == "__main__":
    main()
