"""Write only the Exodus mesh; production input cards are explicit static files."""
from pathlib import Path

import netCDF4
import numpy as np


def names(dataset, variable, dimension, values):
    result = dataset.createVariable(variable, "S1", (dimension, "len_name"))
    result[:] = np.zeros((len(values), 33), dtype="S1")
    for row, value in enumerate(values):
        result[row, :len(value)] = np.array(list(value), dtype="S1")


path = Path(__file__).with_name("rectangle.e")
with netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as dataset:
    dataset.title = "Standard QUAD8 rectangle"
    dataset.api_version = np.float32(8.0)
    dataset.version = np.float32(8.0)
    dataset.floating_point_word_size = np.int32(8)
    dataset.file_size = np.int32(1)
    for name, size in {"len_name": 33, "num_dim": 2, "num_nodes": 8, "num_elem": 1,
                       "num_el_blk": 1, "num_el_in_blk1": 1, "num_nod_per_el1": 8,
                       "num_node_sets": 6, "time_step": None}.items():
        dataset.createDimension(name, size)
    x = [-.01, .01, .01, -.01, 0., .01, 0., -.01]
    y = [-.005, -.005, .005, .005, -.005, 0., .005, 0.]
    dataset.createVariable("coordx", "f8", ("num_nodes",))[:] = x
    dataset.createVariable("coordy", "f8", ("num_nodes",))[:] = y
    names(dataset, "coor_names", "num_dim", ["x", "y"])
    dataset.createVariable("eb_status", "i4", ("num_el_blk",))[:] = [1]
    ids = dataset.createVariable("eb_prop1", "i4", ("num_el_blk",))
    ids.setncattr("name", "ID")
    ids[:] = [1]
    names(dataset, "eb_names", "num_el_blk", ["solid"])
    connectivity = dataset.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
    connectivity.elem_type = "QUAD8"
    connectivity[:] = [list(range(1, 9))]
    ids = dataset.createVariable("ns_prop1", "i4", ("num_node_sets",))
    ids.setncattr("name", "ID")
    ids[:] = list(range(1, 7))
    dataset.createVariable("ns_status", "i4", ("num_node_sets",))[:] = [1] * 6
    sets = [list(range(1, 9)), [1, 2, 3, 4], [1], [2], [3], [4]]
    names(dataset, "ns_names", "num_node_sets", ["field", "corners", "corner1", "corner2", "corner3", "corner4"])
    for index, nodes in enumerate(sets, 1):
        dimension = "num_nod_ns%d" % index
        dataset.createDimension(dimension, len(nodes))
        dataset.createVariable("node_ns%d" % index, "i4", (dimension,))[:] = nodes
