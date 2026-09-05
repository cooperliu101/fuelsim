"""Write geometry only: matching Exodus and Abaqus two-ring contact meshes."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

root = Path(__file__).resolve().parent
coordinates = [(0.001, 0), (0.002, 0), (0.002, 0.001), (0.001, 0.001),
               (0.001, 0.001001), (0.002, 0.001001), (0.002, 0.002001), (0.001, 0.002001)]
blocks = [("lower", [1, 2, 3, 4]), ("upper", [5, 6, 7, 8])]
sides = [("bottom", [1], [1], [1, 2]), ("lower_contact", [1], [3], [3, 4]),
         ("upper_contact", [2], [1], [5, 6]), ("top", [2], [3], [7, 8]),
         ("lower_radial", [1, 1], [2, 4], [1, 2, 3, 4]),
         ("upper_radial", [2, 2], [2, 4], [5, 6, 7, 8])]

with netcdf_file(str(root.parent / "meshes/b7_rz_contact.e"), "w") as f:
    f.api_version = np.float32(7.22)
    f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8)
    f.file_size = np.int32(1)
    f.title = "B7 axisymmetric two-ring thermal mechanical contact"
    for name, size in [("time_step", None), ("len_name", 33), ("num_dim", 2), ("num_nodes", 8),
                       ("num_elem", 2), ("num_el_blk", 2), ("num_side_sets", 6), ("num_node_sets", 6)]:
        f.createDimension(name, size)

    def variable(name, dimensions, data, kind="i"):
        v = f.createVariable(name, kind, dimensions)
        v[:] = data
        return v

    def names(name, dimension, strings):
        data = np.zeros((len(strings), 33), dtype="S1")
        for i, string in enumerate(strings):
            data[i, :len(string)] = np.frombuffer(string.encode("ascii"), dtype="S1")
        variable(name, (dimension, "len_name"), data, "c")

    for name, dimension, count in [("eb", "num_el_blk", 2), ("ss", "num_side_sets", 6), ("ns", "num_node_sets", 6)]:
        ids = variable(name + "_prop1", (dimension,), list(range(1, count + 1)))
        ids.name = "ID"
        variable(name + "_status", (dimension,), [1] * count)
    variable("coordx", ("num_nodes",), [p[0] for p in coordinates], "d")
    variable("coordy", ("num_nodes",), [p[1] for p in coordinates], "d")
    names("eb_names", "num_el_blk", [b[0] for b in blocks])
    names("ss_names", "num_side_sets", [s[0] for s in sides])
    names("ns_names", "num_node_sets", [s[0] for s in sides])
    for i, (_, nodes) in enumerate(blocks, 1):
        f.createDimension("num_el_in_blk%d" % i, 1)
        f.createDimension("num_nod_per_el%d" % i, 4)
        v = variable("connect%d" % i, ("num_el_in_blk%d" % i, "num_nod_per_el%d" % i), [nodes])
        v.elem_type = "QUAD4"
    for i, (_, elements, edges, nodes) in enumerate(sides, 1):
        f.createDimension("num_nod_ns%d" % i, len(nodes))
        variable("node_ns%d" % i, ("num_nod_ns%d" % i,), nodes)
        f.createDimension("num_side_ss%d" % i, len(edges))
        variable("elem_ss%d" % i, ("num_side_ss%d" % i,), elements)
        variable("side_ss%d" % i, ("num_side_ss%d" % i,), edges)

lines = ["** Same geometry and labels as b7_rz_contact.e", "*Node"]
for i, p in enumerate(coordinates, 1):
    lines.append("%d, %.17g, %.17g" % (i, p[0], p[1]))
for i, (name, nodes) in enumerate(blocks, 1):
    lines += ["*Element, type=CAX4T, elset=" + name.upper(), ", ".join(str(n) for n in [i] + nodes)]
for name, _, _, nodes in sides:
    lines += ["*Nset, nset=" + name.upper(), ", ".join(str(n) for n in nodes)]
lines += ["*Nset, nset=RADIAL_FIXED, generate", "1, 8, 1", "*Nset, nset=ALL, generate", "1, 8, 1", "*Surface, type=ELEMENT, name=S_LOWER", "LOWER, S3",
          "*Surface, type=ELEMENT, name=S_UPPER", "UPPER, S1"]
(root / "b7_rz_contact_mesh.inc").write_text("\n".join(lines) + "\n")
