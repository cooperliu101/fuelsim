"""Geometry-only meshes for the five thermal verification studies.

No material, boundary value, solver setting, or production card is generated here.
"""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent
CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
           (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 5),
         (2, 6), (3, 7), (4, 5), (5, 6), (6, 7), (7, 4)]


def write(name, kind, radial_count=1, primary_count=1, curved=False, shared=False, radial_spacing='uniform'):
    axis = kind.startswith('dcax')
    if radial_spacing not in ('uniform', 'geometric') or (radial_spacing == 'geometric' and not axis):
        raise ValueError('Geometric radial spacing is only defined for the cylinder meshes')
    radial_edges = np.geomspace(.01, .03, radial_count+1) if radial_spacing == 'geometric' else None
    quadratic = kind in ('dcax8', 'dc3d20')
    dim = 2 if axis else 3
    local = CORNERS[:4] if axis else CORNERS
    if quadratic:
        local = local + [tuple((CORNERS[a][d] + CORNERS[b][d])/2 for d in range(3))
                         for a, b in (EDGES[:4] if axis else EDGES)]
    nodes, lookup, blocks = [], {}, []
    sides = {'left': [], 'right': []}
    interface = not axis
    if interface and not shared:
        sides.update(secondary=[], primary=[])
    hot = {}
    element = 0
    for block in range(2 if interface else 1):
        cells = []
        ny = primary_count if block else 1
        for ix in range(radial_count):
            for iy in range(ny):
                conn = []
                for u, v, w in local:
                    x = (.01 + .02*(ix+u)/radial_count) if axis else .01*(block+(ix+u)/radial_count)
                    if radial_edges is not None:
                        x = (1-u)*radial_edges[ix]+u*radial_edges[ix+1]
                    y, z = .01*(iy+v)/ny, .01*w
                    if curved:
                        x += .001*((y/.01-.5)**2 + .6*(z/.01-.5)**2)
                    key = (0 if shared else block, round(x, 14), round(y, 14), round(z, 14))
                    if key not in lookup:
                        lookup[key] = len(nodes)+1
                        nodes.append((x, y, z))
                    node = lookup[key]
                    conn.append(node)
                    if block == 0 and ix == 0 and u == 0:
                        # Geometric groups on the hot face; values live in the static cards.
                        group = int(round(2*(y/.01+z/.01)))
                        hot.setdefault('hot'+str(group), set()).add(node)
                element += 1
                cells.append(conn)
                if block == 0 and ix == 0:
                    sides['left'].append((element, 4))
                if (not interface or block == 1) and ix == radial_count-1:
                    sides['right'].append((element, 2))
                if interface and not shared:
                    if block == 0 and ix == radial_count-1:
                        sides['secondary'].append((element, 2))
                    if block == 1 and ix == 0:
                        sides['primary'].append((element, 4))
        blocks.append(cells)
    groups = {'all': set(range(1, len(nodes)+1)), **hot}
    with netcdf_file(str(ROOT/(name+'.e')), 'w') as f:
        f.api_version = np.float32(7.22)
        f.version = np.float32(7.22)
        f.floating_point_word_size = np.int32(8)
        f.file_size = np.int32(1)
        f.title = name
        dims = {'time_step': None, 'len_name': 33, 'num_dim': dim,
                'num_nodes': len(nodes), 'num_elem': element, 'num_el_blk': len(blocks),
                'num_side_sets': len(sides), 'num_node_sets': len(groups)}
        for b, cells in enumerate(blocks, 1):
            dims['num_el_in_blk'+str(b)] = len(cells)
            dims['num_nod_per_el'+str(b)] = len(local)
        for i, entries in enumerate(sides.values(), 1):
            dims['num_side_ss'+str(i)] = len(entries)
        for i, entries in enumerate(groups.values(), 1):
            dims['num_nod_ns'+str(i)] = len(entries)
        for key, value in dims.items():
            f.createDimension(key, value)

        def var(key, dimensions, values, dtype='i'):
            value = f.createVariable(key, dtype, dimensions)
            value[:] = values
            return value

        def names(key, dimension, values):
            array = np.zeros((len(values), 33), dtype='S1')
            for i, value in enumerate(values):
                array[i, :len(value)] = np.frombuffer(value.encode(), dtype='S1')
            var(key, (dimension, 'len_name'), array, 'c')

        for d in range(dim):
            var('coord'+'xyz'[d], ('num_nodes',), [p[d] for p in nodes], 'd')
        for b, cells in enumerate(blocks, 1):
            var('connect'+str(b), ('num_el_in_blk'+str(b), 'num_nod_per_el'+str(b)), cells).elem_type = ('QUAD' if axis else 'HEX')+str(len(local))
        for prefix, dimension, values in [('eb', 'num_el_blk', ['solid'] if axis else ['hot', 'cold']),
                                          ('ss', 'num_side_sets', list(sides)),
                                          ('ns', 'num_node_sets', list(groups))]:
            var(prefix+'_prop1', (dimension,), list(range(1, len(values)+1))).name = 'ID'
            var(prefix+'_status', (dimension,), [1]*len(values))
            names(prefix+'_names', dimension, values)
        for i, entries in enumerate(sides.values(), 1):
            var('elem_ss'+str(i), ('num_side_ss'+str(i),), [e for e, s in entries])
            var('side_ss'+str(i), ('num_side_ss'+str(i),), [s for e, s in entries])
        for i, entries in enumerate(groups.values(), 1):
            var('node_ns'+str(i), ('num_nod_ns'+str(i),), sorted(entries))
    with (ROOT/(name+'_mesh.inp')).open('w') as f:
        f.write('*NODE, NSET=ALL\n')
        for i, point in enumerate(nodes, 1):
            f.write(str(i)+', '+', '.join(format(v, '.17g') for v in point[:dim])+'\n')
        element = 0
        for b, cells in enumerate(blocks):
            f.write('*ELEMENT, TYPE='+kind.upper()+', ELSET='+('SOLID' if axis else ['HOT', 'COLD'][b])+'\n')
            for conn in cells:
                element += 1
                c = conn[:12]+conn[16:20]+conn[12:16] if len(conn) == 20 else conn
                f.write(str(element)+', '+', '.join(map(str, c[:15]))+(',\n' if len(c)>15 else '\n'))
                if len(c)>15:
                    f.write(', '.join(map(str, c[15:]))+'\n')
        if interface:
            f.write('*ELSET, ELSET=SOLID\nHOT, COLD\n')
        for key, entries in hot.items():
            f.write('*NSET, NSET='+key.upper()+'\n'+', '.join(map(str, sorted(entries)))+'\n')
        all_cells = [c for cells in blocks for c in cells]
        face = [1, 2] if axis else [1, 2, 6, 5]
        if quadratic:
            face += [5] if axis else [9, 14, 17, 13]
        right = sorted({all_cells[e-1][j] for e, s in sides['right'] for j in face})
        f.write('*NSET, NSET=RIGHT\n'+', '.join(map(str, right))+'\n')
        if interface and not shared:
            for key, label in [('secondary', 'S4'), ('primary', 'S6')]:
                f.write('*SURFACE, NAME='+key.upper()+', TYPE=ELEMENT\n')
                for e, s in sides[key]:
                    f.write(str(e)+', '+label+'\n')


if __name__ == '__main__':
    for kind in ['dcax4', 'dcax8']:
        for count in [4, 8, 16]:
            write('study_'+kind+'_n'+str(count), kind, radial_count=count)
    write('study_dc3d20_nonmatching', 'dc3d20', primary_count=2)
    write('study_dc3d20_curved', 'dc3d20', curved=True)
    write('study_dc3d20_curved_nonmatching', 'dc3d20', primary_count=2, curved=True)
    write('study_dc3d8_shared', 'dc3d8', shared=True)
    write('study_dc3d8_interface', 'dc3d8')
    for count in [64, 128]:
        write('study_dcax4_n'+str(count), 'dcax4', radial_count=count)
    for count in [16, 64]:
        write('study_dcax4_geometric_n'+str(count), 'dcax4', radial_count=count, radial_spacing='geometric')
