"""Write only the geometry of the rectangular and distorted CAX4RT probes."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file
root = Path(__file__).resolve().parent.parent / 'meshes'
for name, points in [('b9_rz_probe', [(1, 0), (2.1, .1), (1.9, 1.2), (.9, 1)]),
                     ('b9_rz_rectangle', [(.001, 0), (.002, 0), (.002, .001), (.001, .001)])]:
    with netcdf_file(str(root / (name + '.e')), 'w') as f:
        f.api_version = np.float32(7.22)
        f.version = np.float32(7.22)
        f.floating_point_word_size = np.int32(8)
        f.file_size = np.int32(1)
        f.title = name
        for dim, size in [('time_step', None), ('len_name', 33), ('num_dim', 2), ('num_nodes', 4),
                          ('num_elem', 1), ('num_el_blk', 1), ('num_node_sets', 4),
                          ('num_el_in_blk1', 1), ('num_nod_per_el1', 4)]:
            f.createDimension(dim, size)
        def var(name, dims, values, kind='i'):
            v = f.createVariable(name, kind, dims)
            v[:] = values
            return v
        def names(name, dim, values):
            values_array = np.zeros((len(values), 33), dtype='S1')
            for i, value in enumerate(values):
                values_array[i, :len(value)] = np.frombuffer(value.encode('ascii'), dtype='S1')
            var(name, (dim, 'len_name'), values_array, 'c')
        var('coordx', ('num_nodes',), [p[0] for p in points], 'd')
        var('coordy', ('num_nodes',), [p[1] for p in points], 'd')
        var('connect1', ('num_el_in_blk1', 'num_nod_per_el1'), [[1, 2, 3, 4]]).elem_type = 'QUAD4'
        var('eb_prop1', ('num_el_blk',), [1]).name = 'ID'
        var('eb_status', ('num_el_blk',), [1])
        names('eb_names', 'num_el_blk', ['solid'])
        var('ns_prop1', ('num_node_sets',), [1, 2, 3, 4]).name = 'ID'
        var('ns_status', ('num_node_sets',), [1, 1, 1, 1])
        names('ns_names', 'num_node_sets', ['n1', 'n2', 'n3', 'n4'])
        for n in range(1, 5):
            dim = 'num_nod_ns%d' % n
            f.createDimension(dim, 1)
            var('node_ns%d' % n, (dim,), [n])
