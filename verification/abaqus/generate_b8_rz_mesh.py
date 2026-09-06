"""Geometry only: deformable coaxial rings with nonmatching sliding surfaces."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

root = Path(__file__).resolve().parent
coordinates = []
blocks = []
boundaries = []
for name, radii, heights in [
        ("inner", (0.001, 0.002), (0.00125, 0.00175, 0.00225)),
        ("outer", (0.002001, 0.003001), (0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006))]:
    node_start = len(coordinates) + 1
    element_start = sum(len(b[1]) for b in blocks) + 1
    for z in heights:
        coordinates.extend((r, z) for r in radii)
    elements = [[node_start+2*j, node_start+2*j+1, node_start+2*j+3, node_start+2*j+2]
                for j in range(len(heights)-1)]
    blocks.append((name, elements))
    for side, local_side, offset in [("left", 4, 0), ("right", 2, 1)]:
        boundaries.append((name+"_"+side, list(range(element_start, element_start+len(elements))),
                           [local_side]*len(elements), [node_start+2*j+offset for j in range(len(heights))]))
    boundaries.append((name+"_bottom", [element_start], [1], [node_start, node_start+1]))
    boundaries.append((name+"_top", [element_start+len(elements)-1], [3],
                       [node_start+2*(len(heights)-1), node_start+2*(len(heights)-1)+1]))

with netcdf_file(str(root.parent/"meshes/b8_rz_sliding.e"), "w") as f:
    f.api_version = np.float32(7.22)
    f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8)
    f.file_size = np.int32(1)
    f.title = "B8 deformable coaxial rings with finite sliding"
    for name, size in [("time_step", None), ("len_name", 33), ("num_dim", 2),
                       ("num_nodes", len(coordinates)), ("num_elem", 8), ("num_el_blk", 2),
                       ("num_side_sets", len(boundaries)), ("num_node_sets", len(boundaries))]:
        f.createDimension(name, size)

    def variable(name, dims, data, kind="i"):
        v = f.createVariable(name, kind, dims)
        v[:] = data
        return v

    def names(name, dimension, strings):
        data = np.zeros((len(strings), 33), dtype="S1")
        for i, string in enumerate(strings):
            data[i, :len(string)] = np.frombuffer(string.encode("ascii"), dtype="S1")
        variable(name, (dimension, "len_name"), data, "c")

    for name, dim, count in [("eb", "num_el_blk", 2), ("ss", "num_side_sets", 8), ("ns", "num_node_sets", 8)]:
        v = variable(name+"_prop1", (dim,), list(range(1, count+1)))
        v.name = "ID"
        variable(name+"_status", (dim,), [1]*count)
    variable("coordx", ("num_nodes",), [p[0] for p in coordinates], "d")
    variable("coordy", ("num_nodes",), [p[1] for p in coordinates], "d")
    names("eb_names", "num_el_blk", [b[0] for b in blocks])
    names("ss_names", "num_side_sets", [b[0] for b in boundaries])
    names("ns_names", "num_node_sets", [b[0] for b in boundaries])
    for i, (_, elements) in enumerate(blocks, 1):
        f.createDimension("num_el_in_blk%d" % i, len(elements))
        f.createDimension("num_nod_per_el%d" % i, 4)
        v = variable("connect%d" % i, ("num_el_in_blk%d" % i, "num_nod_per_el%d" % i), elements)
        v.elem_type = "QUAD4"
    for i, (_, elements, edges, nodes) in enumerate(boundaries, 1):
        f.createDimension("num_side_ss%d" % i, len(edges))
        f.createDimension("num_nod_ns%d" % i, len(nodes))
        variable("elem_ss%d" % i, ("num_side_ss%d" % i,), elements)
        variable("side_ss%d" % i, ("num_side_ss%d" % i,), edges)
        variable("node_ns%d" % i, ("num_nod_ns%d" % i,), nodes)

lines = ["** Geometry only; same coordinates and connectivity as b8_rz_sliding.e", "*Node"]
for i, p in enumerate(coordinates, 1):
    lines.append("%d, %.17g, %.17g" % (i, p[0], p[1]))
element = 0
for name, elements in blocks:
    lines.append("*Element, type=CAX4T, elset="+name.upper())
    for nodes in elements:
        element += 1
        lines.append(", ".join(str(n) for n in [element]+nodes))
for name, _, _, nodes in boundaries:
    lines += ["*Nset, nset="+name.upper(), ", ".join(str(n) for n in nodes)]
lines += ["*Nset, nset=ALL, generate", "1, %d, 1" % len(coordinates),
          "*Elset, elset=ALL_ELEMENTS", "INNER, OUTER", "*Surface, type=ELEMENT, name=S_INNER", "INNER, S2",
          "*Surface, type=ELEMENT, name=S_OUTER", "OUTER, S4"]
(root/"b8_rz_sliding_mesh.inc").write_text("\n".join(lines)+"\n")
