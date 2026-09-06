"""Write geometry only for the independent CAX4RT thermal operator checks."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

points = [(2, 0), (3, 0), (3, 1), (2, 1),
          (1, 0), (3, 0), (3, 1), (1, 1),
          (1, 0), (2, 0), (2, 2), (1, 2),
          (1, 0), (2, 0), (2.5, 1), (1.5, 1),
          (1, 0), (2, 0), (2, 1), (1, 1),
          (1, 3), (2, 3), (2.5, 4), (1.5, 4)]
connectivity = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12],
                [13, 14, 15, 16], [17, 18, 19, 20], [22, 23, 24, 21]]
sets = [list(range(1, 25))] + [list(range(n, 25, 4)) for n in range(1, 5)]
target = Path(__file__).resolve().parent.parent / 'meshes' / 'b915_rz_thermal.e'
with netcdf_file(str(target), 'w') as f:
    f.api_version = np.float32(7.22)
    f.version = np.float32(7.22)
    f.floating_point_word_size = np.int32(8)
    f.file_size = np.int32(1)
    f.title = 'b915 CAX4RT thermal geometry'
    for dim, size in [('time_step', None), ('len_name', 33), ('num_dim', 2), ('num_nodes', 24),
                      ('num_elem', 6), ('num_el_blk', 1), ('num_node_sets', 5),
                      ('num_el_in_blk1', 6), ('num_nod_per_el1', 4)]:
        f.createDimension(dim, size)

    def var(name, dims, values, kind='i'):
        result = f.createVariable(name, kind, dims)
        result[:] = values
        return result

    def names(name, dim, values):
        data = np.zeros((len(values), 33), dtype='S1')
        for i, value in enumerate(values):
            data[i, :len(value)] = np.frombuffer(value.encode('ascii'), dtype='S1')
        var(name, (dim, 'len_name'), data, 'c')

    var('coordx', ('num_nodes',), [p[0] for p in points], 'd')
    var('coordy', ('num_nodes',), [p[1] for p in points], 'd')
    var('connect1', ('num_el_in_blk1', 'num_nod_per_el1'), connectivity).elem_type = 'QUAD4'
    var('eb_prop1', ('num_el_blk',), [1]).name = 'ID'
    var('eb_status', ('num_el_blk',), [1])
    names('eb_names', 'num_el_blk', ['solid'])
    var('ns_prop1', ('num_node_sets',), [1, 2, 3, 4, 5]).name = 'ID'
    var('ns_status', ('num_node_sets',), [1, 1, 1, 1, 1])
    names('ns_names', 'num_node_sets', ['all', 't1', 't2', 't3', 't4'])
    for i, nodes in enumerate(sets, 1):
        dim = 'num_nod_ns%d' % i
        f.createDimension(dim, len(nodes))
        var('node_ns%d' % i, (dim,), nodes)
