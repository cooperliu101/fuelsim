"""Write geometry only for the mixed-order CAX8T identification probes."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

root = Path(__file__).resolve().parent.parent / "meshes"
points = [(1, 0), (2.1, .1), (1.9, 1.2), (.9, 1),
          (1.55, .03), (2.02, .65), (1.4, 1.12), (.93, .5)]
with netcdf_file(str(root / "b10_rz_probe.e"), "w") as f:
    f.api_version = np.float32(7.22)
    f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8)
    f.file_size = np.int32(1)
    f.title = "CAX8T curved probe geometry"
    for name, size in [("time_step", None), ("len_name", 33), ("num_dim", 2), ("num_nodes", 8),
                       ("num_elem", 1), ("num_el_blk", 1), ("num_node_sets", 8),
                       ("num_el_in_blk1", 1), ("num_nod_per_el1", 8)]:
        f.createDimension(name, size)

    def variable(name, dims, values, kind="i"):
        v = f.createVariable(name, kind, dims)
        v[:] = values
        return v

    def names(name, dim, values):
        array = np.zeros((len(values), 33), dtype="S1")
        for i, value in enumerate(values):
            array[i, :len(value)] = np.frombuffer(value.encode("ascii"), dtype="S1")
        variable(name, (dim, "len_name"), array, "c")

    variable("coordx", ("num_nodes",), [p[0] for p in points], "d")
    variable("coordy", ("num_nodes",), [p[1] for p in points], "d")
    variable("connect1", ("num_el_in_blk1", "num_nod_per_el1"), [list(range(1, 9))]).elem_type = "QUAD8"
    variable("eb_prop1", ("num_el_blk",), [1]).name = "ID"
    variable("eb_status", ("num_el_blk",), [1])
    names("eb_names", "num_el_blk", ["solid"])
    variable("ns_prop1", ("num_node_sets",), list(range(1, 9))).name = "ID"
    variable("ns_status", ("num_node_sets",), [1]*8)
    names("ns_names", "num_node_sets", ["n%d" % n for n in range(1, 9)])
    for n in range(1, 9):
        dim = "num_nod_ns%d" % n
        f.createDimension(dim, 1)
        variable("node_ns%d" % n, (dim,), [n])
