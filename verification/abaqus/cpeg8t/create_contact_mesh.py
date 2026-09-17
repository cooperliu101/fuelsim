"""Generate the two independent QUAD8 blocks and their named Exodus boundaries."""
from pathlib import Path
import netCDF4
import numpy as np


def names(dataset, variable, dimension, values):
    result = dataset.createVariable(variable, "S1", (dimension, "len_name"))
    result[:] = np.zeros((len(values), 33), dtype="S1")
    for row, value in enumerate(values):
        result[row, :len(value)] = np.array(list(value), dtype="S1")


with netCDF4.Dataset(Path(__file__).with_name("contact.e"), "w", format="NETCDF3_64BIT_OFFSET") as d:
    d.title = "CPEG8T two-block thermal mechanical contact"
    d.api_version = np.float32(8)
    d.version = np.float32(8)
    d.floating_point_word_size = np.int32(8)
    d.file_size = np.int32(1)
    dimensions = {"len_name": 33, "num_dim": 2, "num_nodes": 16, "num_elem": 2,
                  "num_el_blk": 2, "num_node_sets": 3, "num_side_sets": 2, "time_step": None}
    for key, value in dimensions.items():
        d.createDimension(key, value)
    x = [-.01, .01, .01, -.01, 0, .01, 0, -.01]
    y = [-.01, -.01, 0, 0, -.01, -.005, 0, -.005]
    d.createVariable("coordx", "f8", ("num_nodes",))[:] = x + x
    d.createVariable("coordy", "f8", ("num_nodes",))[:] = y + [v + .0101 for v in y]
    names(d, "coor_names", "num_dim", ["x", "y"])
    for prefix, dim, labels in [("eb", "num_el_blk", ["lower", "upper"]),
                                ("ns", "num_node_sets", ["field", "bottom", "top"]),
                                ("ss", "num_side_sets", ["primary", "secondary"])]:
        ids = d.createVariable(prefix + "_prop1", "i4", (dim,))
        ids.setncattr("name", "ID")
        ids[:] = range(1, len(labels) + 1)
        d.createVariable(prefix + "_status", "i4", (dim,))[:] = [1] * len(labels)
        names(d, prefix + "_names", dim, labels)
    for i in (1, 2):
        ed, nd = "num_el_in_blk%d" % i, "num_nod_per_el%d" % i
        for key, value in [(ed, 1), (nd, 8)]:
            d.createDimension(key, value)
        connectivity = d.createVariable("connect%d" % i, "i4", (ed, nd))
        connectivity.elem_type = "QUAD8"
        connectivity[:] = [list(range(1 + 8 * (i - 1), 9 + 8 * (i - 1)))]
    for i, nodes in enumerate([list(range(1, 17)), [1, 2, 5], [11, 12, 15]], 1):
        dim = "num_nod_ns%d" % i
        d.createDimension(dim, len(nodes))
        d.createVariable("node_ns%d" % i, "i4", (dim,))[:] = nodes
    for i, element, side in [(1, 1, 3), (2, 2, 1)]:
        dim = "num_side_ss%d" % i
        d.createDimension(dim, 1)
        d.createVariable("elem_ss%d" % i, "i4", (dim,))[:] = [element]
        d.createVariable("side_ss%d" % i, "i4", (dim,))[:] = [side]
