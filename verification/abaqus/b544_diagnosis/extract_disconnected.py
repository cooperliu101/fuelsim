from __future__ import print_function

import csv
import sys
from odbAccess import openOdb

SOURCE_NODES = (
    (1, 2, 7, 6, 11, 12, 17, 16),
    (2, 3, 8, 7, 12, 13, 18, 17),
    (3, 4, 9, 8, 13, 14, 19, 18),
    (4, 5, 10, 9, 14, 15, 20, 19),
    (11, 12, 17, 16, 21, 22, 27, 26),
    (12, 13, 18, 17, 22, 23, 28, 27),
    (13, 14, 19, 18, 23, 24, 29, 28),
    (14, 15, 20, 19, 24, 25, 30, 29),
)


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


def field(frame, name):
    if name in frame.fieldOutputs:
        return frame.fieldOutputs[name]
    matches = [key for key in frame.fieldOutputs.keys() if key.strip().startswith(name)]
    if len(matches) != 1:
        raise RuntimeError("Expected one native field %s" % name)
    return frame.fieldOutputs[matches[0]]


def nodal(frame, name):
    result = dict((value.nodeLabel, data(value)) for value in field(frame, name).values if value.nodeLabel)
    if set(result.keys()) != set(range(1, 65)):
        raise RuntimeError("Incomplete native node coverage for %s" % name)
    return result


def points(frame, name):
    result = dict(((value.elementLabel, value.integrationPoint), value)
                  for value in field(frame, name).values if value.elementLabel and value.integrationPoint)
    if set(result.keys()) != set((element, 1) for element in range(1, 9)):
        raise RuntimeError("Incomplete native material-point coverage for %s" % name)
    return result


def orientation(value):
    try:
        matrix = value.localCoordSystemDouble
    except Exception:
        matrix = value.localCoordSystem
    if matrix is None or len(matrix) == 0:
        return [0] + [""] * 9
    return [1] + [component for row in matrix for component in row]


def write(writer, values):
    writer.writerow(["%.17g" % value if isinstance(value, float) else value for value in values])


if len(sys.argv) != 5:
    raise RuntimeError("usage: extract_disconnected.py job.odb nodal.csv integration.csv energy.csv")

odb = openOdb(path=sys.argv[1], readOnly=True)
try:
    step = odb.steps["PATH"]
    if len(step.frames) != 11:
        raise RuntimeError("Expected initial frame and ten accepted increments")
    histories = [region.historyOutputs["ALLAE"].data for region in step.historyRegions.values()
                 if "ALLAE" in region.historyOutputs]
    if len(histories) != 1 or len(histories[0]) != 11:
        raise RuntimeError("Expected one complete native ALLAE history")
    with open(sys.argv[2], "wb") as node_file, open(sys.argv[3], "wb") as point_file, open(sys.argv[4], "wb") as energy_file:
        nodes_writer, points_writer, energy_writer = csv.writer(node_file), csv.writer(point_file), csv.writer(energy_file)
        nodes_writer.writerow(["increment", "time_s", "element", "local_node", "node", "source_node",
                               "temperature_k", "u1_m", "u2_m", "u3_m", "rf1_n", "rf2_n", "rf3_n"])
        points_writer.writerow(["increment", "time_s", "element", "integration_point", "temperature_k",
                                "x_m", "y_m", "z_m", "s11_pa", "s22_pa", "s33_pa", "s12_pa", "s13_pa", "s23_pa",
                                "ivol_m3", "local_csys_present", "local_csys_11", "local_csys_12", "local_csys_13",
                                "local_csys_21", "local_csys_22", "local_csys_23", "local_csys_31", "local_csys_32", "local_csys_33"])
        energy_writer.writerow(["increment", "time_s", "allae_j"])
        for increment in range(1, 11):
            frame = step.frames[increment]
            if abs(frame.frameValue - increment * 100000.0) > 1e-7:
                raise RuntimeError("Unexpected native increment time")
            temperatures, displacements, reactions = nodal(frame, "NT11"), nodal(frame, "U"), nodal(frame, "RF")
            if tuple(field(frame, "S").componentLabels) != ("S11", "S22", "S33", "S12", "S13", "S23"):
                raise RuntimeError("Unexpected native stress component order")
            stresses, volumes, coordinates, material_temperatures = (
                points(frame, name) for name in ("S", "IVOL", "COORD", "TEMP"))
            for element in range(1, 9):
                for local in range(1, 9):
                    node = 8 * (element - 1) + local
                    write(nodes_writer, [increment, frame.frameValue, element, local, node, SOURCE_NODES[element - 1][local - 1],
                                         temperatures[node]] + list(displacements[node]) + list(reactions[node]))
                key = (element, 1)
                write(points_writer, [increment, frame.frameValue, element, 1, data(material_temperatures[key])]
                      + list(data(coordinates[key])) + list(data(stresses[key])) + [data(volumes[key])]
                      + orientation(stresses[key]))
            time, value = histories[0][increment]
            if abs(time - frame.frameValue) > 1e-7:
                raise RuntimeError("Native ALLAE history time differs from field output")
            write(energy_writer, [increment, frame.frameValue, value])
    open(sys.argv[2] + ".ok", "wb").close()
finally:
    odb.close()
