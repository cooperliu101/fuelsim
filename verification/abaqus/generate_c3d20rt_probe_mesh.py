"""Write only the Exodus mesh for the tracked C3D20RT operator input cards."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

root = Path(__file__).resolve().parent
text = (root / "c3d20rt_thermal_probe.inp").read_text()
rows = text.split("*Node\n", 1)[1].split("*Element", 1)[0].strip().splitlines()
coordinates = np.array([[float(v) for v in row.split(",")[1:]] for row in rows])
node_sets = [("all", list(range(1, 21)))]
node_sets += [("xzero", [i + 1 for i, p in enumerate(coordinates) if p[0] == 0]),
              ("xhalf", [i + 1 for i, p in enumerate(coordinates) if p[0] == 0.5]),
              ("xone", [i + 1 for i, p in enumerate(coordinates) if p[0] == 1])]
node_sets += [("n%d" % i, [i]) for i in range(1, 9)]
node_sets += [("yzero", [i + 1 for i, p in enumerate(coordinates) if p[1] == 0]),
              ("zzero", [i + 1 for i, p in enumerate(coordinates) if p[2] == 0])]
with netcdf_file(str(root.parent / "meshes/c3d20rt_probe.e"), "w") as f:
    f.api_version = np.float32(7.22)
    f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8)
    f.file_size = np.int32(1)
    f.title = "C3D20RT independent thermal operator mesh"
    for name, size in [("time_step", None), ("len_name", 33), ("num_dim", 3),
                       ("num_nodes", 20), ("num_elem", 1), ("num_el_blk", 1),
                       ("num_node_sets", len(node_sets)), ("num_el_in_blk1", 1),
                       ("num_nod_per_el1", 20)]:
        f.createDimension(name, size)
    def put(name, dims, values, kind="i"):
        v = f.createVariable(name, kind, dims)
        v[:] = values
        return v
    def names(name, dim, values):
        data = np.zeros((len(values), 33), dtype="S1")
        for i, value in enumerate(values):
            data[i, :len(value)] = np.frombuffer(value.encode(), dtype="S1")
        put(name, (dim, "len_name"), data, "c")
    for i, axis in enumerate("xyz"):
        put("coord" + axis, ("num_nodes",), coordinates[:, i], "d")
    for prefix, dim, count in [("eb", "num_el_blk", 1), ("ns", "num_node_sets", len(node_sets))]:
        put(prefix + "_prop1", (dim,), list(range(1, count + 1))).name = "ID"
        put(prefix + "_status", (dim,), [1] * count)
    names("eb_names", "num_el_blk", ["solid"])
    names("ns_names", "num_node_sets", [name for name, _ in node_sets])
    # Exodus/Fuelsim orders vertical edge nodes before the upper face edges.
    connectivity = list(range(1, 13)) + list(range(17, 21)) + list(range(13, 17))
    put("connect1", ("num_el_in_blk1", "num_nod_per_el1"), [connectivity]).elem_type = "HEX20"
    for i, (_, nodes) in enumerate(node_sets, 1):
        dim = "num_nod_ns%d" % i
        f.createDimension(dim, len(nodes))
        put("node_ns%d" % i, (dim,), nodes)
