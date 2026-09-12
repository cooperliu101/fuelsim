"""Rebuild the committed BAR2 geometry fixtures; never create or edit input cards."""
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def names(dataset, variable, dimension, values):
    data = np.zeros((len(values), 64), dtype="S1")
    for row, value in enumerate(values):
        encoded = value.encode("ascii")
        data[row, :len(encoded)] = np.frombuffer(encoded, dtype="S1")
    dataset.createVariable(variable, "S1", (dimension, "len_name"))[:] = data


def write(path, nodes, blocks, node_sets, side_sets):
    with Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as data:
        data.api_version = np.float32(8.11)
        data.version = np.float32(8.11)
        data.floating_point_word_size = np.int32(8)
        data.file_size = np.int32(1)
        data.title = "Fuelsim BAR2 radial slices with explicit axial control nodes"
        for name, count in (("num_dim", 2), ("num_nodes", len(nodes)), ("len_name", 64),
                            ("num_elem", sum(len(b[1]) for b in blocks)), ("num_el_blk", len(blocks)),
                            ("num_node_sets", len(node_sets)), ("num_side_sets", len(side_sets))):
            data.createDimension(name, count)
        data.createDimension("time_step", None)
        for column, name in enumerate(("coordx", "coordy")):
            data.createVariable(name, "f8", ("num_nodes",))[:] = np.asarray(nodes)[:, column]
        names(data, "coor_names", "num_dim", ["r", "z"])
        for prefix, dimension, items in (("eb", "num_el_blk", blocks),
                                          ("ns", "num_node_sets", node_sets),
                                          ("ss", "num_side_sets", side_sets)):
            prop = data.createVariable(prefix + "_prop1", "i4", (dimension,))
            prop.setncattr("name", "ID")
            prop[:] = np.arange(1, len(items) + 1)
            data.createVariable(prefix + "_status", "i4", (dimension,))[:] = 1
            names(data, prefix + "_names", dimension, [item[0] for item in items])
        for block, (_, elements) in enumerate(blocks, 1):
            count, width, attributes = f"num_el_in_blk{block}", f"num_nod_per_el{block}", f"num_att_in_blk{block}"
            data.createDimension(count, len(elements))
            data.createDimension(width, 2)
            data.createDimension(attributes, 2)
            connectivity = data.createVariable(f"connect{block}", "i4", (count, width))
            connectivity.elem_type = "BAR2"
            connectivity[:] = [element[:2] for element in elements]
            data.createVariable(f"attrib{block}", "f8", (count, attributes))[:] = [e[2:] for e in elements]
            names(data, f"attrib_name{block}", attributes, ["axial_lower_node", "axial_upper_node"])
        for index, (_, nodeset) in enumerate(node_sets, 1):
            dimension = f"num_nod_ns{index}"
            data.createDimension(dimension, len(nodeset))
            data.createVariable(f"node_ns{index}", "i4", (dimension,))[:] = nodeset
        for index, (_, sides) in enumerate(side_sets, 1):
            dimension = f"num_side_ss{index}"
            data.createDimension(dimension, len(sides))
            data.createVariable(f"elem_ss{index}", "i4", (dimension,))[:] = [s[0] for s in sides]
            data.createVariable(f"side_ss{index}", "i4", (dimension,))[:] = [s[1] for s in sides]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    write(root / "gps_uniform.e", [(0.004, .005), (.005, .005), (0, 0), (0, .01)],
          [("solid", [(1, 2, 3, 4)])],
          [("all_radial", [1, 2]), ("bottom", [3]), ("top", [4])],
          [("inner", [(1, 1)]), ("outer", [(1, 2)])])
    write(root / "gps_two_slice.e",
          [(.004, .005), (.005, .005), (.00501, .005), (.006, .005),
           (.004, .02), (.005, .02), (.00501, .02), (.006, .02),
           (0, 0), (0, .01), (0, .03), (0, 0), (0, .01), (0, .03)],
          [("inner", [(1, 2, 9, 10), (5, 6, 10, 11)]),
           ("outer", [(3, 4, 12, 13), (7, 8, 13, 14)])],
          [("inner_radial", [1, 2, 5, 6]), ("outer_radial", [3, 4, 7, 8]),
           ("inner_axial", [9, 10, 11]), ("outer_axial", [12, 13, 14]),
           ("inner_bottom", [9]), ("inner_top", [11]), ("outer_bottom", [12]), ("outer_top", [14]),
           ("inner_layer1", [1, 2]), ("inner_layer2", [5, 6]),
           ("outer_layer1", [3, 4]), ("outer_layer2", [7, 8])],
          [("inner_inside", [(1, 1), (2, 1)]), ("inner_interface", [(1, 2), (2, 2)]),
           ("outer_interface", [(3, 1), (4, 1)]), ("outer_outside", [(3, 2), (4, 2)]),
           ("inner_interface1", [(1, 2)]), ("inner_interface2", [(2, 2)]),
           ("outer_interface1", [(3, 1)]), ("outer_interface2", [(4, 1)])])
