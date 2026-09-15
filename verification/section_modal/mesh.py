"""Generate only the committed extrusion mesh, never create or edit input cards."""
from pathlib import Path
import numpy as np
from netCDF4 import Dataset


def names(data, variable, dimension, values):
    out = np.zeros((len(values), 64), dtype='S1')
    for i, value in enumerate(values):
        out[i, :len(value)] = np.frombuffer(value.encode(), dtype='S1')
    data.createVariable(variable, 'S1', (dimension, 'len_name'))[:] = out


def write(path, nx=2, ny=4, nz=8, graded=False, local_interfaces=False, coarsen_interior=False,
          width_lines=False, point_support=False, axial_planes=None, solid_interfaces=None):
    # Exodus HEX20 order: corner nodes, lower edges, vertical edges, upper edges.
    offsets = [(0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2),
               (1,0,0),(2,1,0),(1,2,0),(0,1,0),(0,0,1),(2,0,1),(2,2,1),(0,2,1),
               (1,0,2),(2,1,2),(1,2,2),(0,1,2)]
    ids, nodes, elements = {}, [], []
    s = np.linspace(0.0, 1.0, nz + 1)
    axial = 0.2 * (s - 0.85 * np.sin(2 * np.pi * s) / (2 * np.pi))
    axial[0], axial[-1] = 0.0, 0.2
    if nz % 2 == 0:
        axial[nz // 2] = 0.1
    if coarsen_interior:
        if not graded or nz < 32:
            raise ValueError('Interior coarsening requires a resolved graded axial mesh')
        # Preserve the original first and last eight cells, where the clamped
        # corner and applied end traction need fine axial resolution. Quintic
        # modal amplitudes use every second plane in the remaining smooth span.
        keep = np.array([i <= 8 or i >= nz - 8 or i % 2 == 0 for i in range(nz + 1)])
        axial = axial[keep]
        nz = len(axial) - 1
    if axial_planes is not None:
        axial = np.asarray(axial_planes, dtype=float)
        if (axial.ndim != 1 or len(axial) < 2 or not np.isfinite(axial).all()
                or axial[0] != 0 or axial[-1] != 0.2 or np.any(np.diff(axial) <= 0)):
            raise ValueError('Explicit extrusion planes must increase from zero to 0.2 m')
        nz, graded = len(axial) - 1, True
    sides = {key: [] for key in ('end', 'end_plus', 'end_minus', 'broad', 'half_broad')}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = []
                for a,b,c in offsets:
                    key = (2*i+a, 2*j+b, 2*k+c)
                    if key not in ids:
                        ids[key] = len(nodes)+1
                        z = 0.2 * key[2] / (2 * nz)
                        if graded:
                            layer, midpoint = divmod(key[2], 2)
                            z = axial[layer] if not midpoint else (axial[layer] + axial[layer + 1]) / 2
                        nodes.append((-0.001+0.002*key[0]/(2*nx), -0.01+0.02*key[1]/(2*ny), z))
                    cell.append(ids[key])
                elements.append(cell)
                e = len(elements)
                if k == nz-1:
                    sides['end'].append((e,6))
                    sides['end_plus' if i >= nx//2 else 'end_minus'].append((e,6))
                if i == nx-1:
                    sides['broad'].append((e,2))
                    if j >= ny//2: sides['half_broad'].append((e,2))
    sets = {'root': [i+1 for i,p in enumerate(nodes) if p[2] == 0], 'all_nodes': list(range(1,len(nodes)+1))}
    if local_interfaces:
        end_planes = axial if graded else np.linspace(0.0, 0.2, nz + 1)
        for name, target in [('lower_interface', 0.03), ('upper_interface', 0.17),
                             ('lower_fine', 0.006), ('upper_fine', 0.194)]:
            plane = end_planes[np.argmin(abs(end_planes - target))]
            sets[name] = [i + 1 for i, p in enumerate(nodes) if p[2] == plane]
    if solid_interfaces is not None:
        for name, plane in zip(('solid_lower', 'solid_upper'), solid_interfaces):
            if plane not in axial:
                raise ValueError('A solid interface must be an existing extrusion plane')
            sets[name] = [i + 1 for i, p in enumerate(nodes) if p[2] == plane]
    if width_lines:
        if ny % 4:
            raise ValueError('The width study requires four mesh-aligned strips')
        for j in range(5):
            y = nodes[ids[(0, j * ny // 2, 0)] - 1][1]
            sets['width_' + str(j)] = [i + 1 for i, p in enumerate(nodes) if p[2] == 0 and p[1] == y]
    if point_support:
        for name, key in [('support_a', (nx, ny, 0)), ('support_b', (2 * nx, ny, 0)),
                          ('support_c', (nx, 2 * ny, 0))]:
            if key not in ids:
                raise ValueError('The support point must coincide with an existing mesh node')
            sets[name] = [ids[key]]
    with Dataset(path, 'w', format='NETCDF3_64BIT_OFFSET') as data:
        data.api_version = np.float32(8.11); data.version = np.float32(8.11)
        data.floating_point_word_size = np.int32(8); data.file_size = np.int32(1)
        data.title = 'Fixed plate extrusion for identical solid and section-modal inputs'
        for name,size in [('num_dim',3),('num_nodes',len(nodes)),('num_elem',len(elements)),('num_el_blk',1),
                          ('len_name',64),('num_node_sets',len(sets)),('num_side_sets',len(sides))]:
            data.createDimension(name,size)
        data.createDimension('time_step',None)
        for c,name in enumerate(['coordx','coordy','coordz']):
            data.createVariable(name,'f8',('num_nodes',))[:] = np.array(nodes)[:,c]
        names(data,'coor_names','num_dim',['x','y','z'])
        for prefix,dimension,items in [('eb','num_el_blk',{'plate':[]}),('ns','num_node_sets',sets),('ss','num_side_sets',sides)]:
            data.createVariable(prefix+'_status','i4',(dimension,))[:] = 1
            prop = data.createVariable(prefix+'_prop1','i4',(dimension,)); prop.setncattr('name','ID'); prop[:] = range(1,len(items)+1)
            names(data,prefix+'_names',dimension,list(items))
        data.createDimension('num_el_in_blk1',len(elements)); data.createDimension('num_nod_per_el1',20)
        conn=data.createVariable('connect1','i4',('num_el_in_blk1','num_nod_per_el1')); conn.elem_type='HEX20'; conn[:]=elements
        for i,(_,values) in enumerate(sets.items(),1):
            dim=f'num_nod_ns{i}';data.createDimension(dim,len(values));data.createVariable(f'node_ns{i}','i4',(dim,))[:]=values
        for i,(_,values) in enumerate(sides.items(),1):
            dim=f'num_side_ss{i}';data.createDimension(dim,len(values))
            data.createVariable(f'elem_ss{i}','i4',(dim,))[:]=[v[0] for v in values]
            data.createVariable(f'side_ss{i}','i4',(dim,))[:]=[v[1] for v in values]

def write_solid_ends(path, end_layers, interior_elements=32, tip_layers=None):
    tip_layers = end_layers if tip_layers is None else tip_layers
    if end_layers <= 0 or tip_layers <= 0 or end_layers + tip_layers >= 1024 or interior_elements <= 0:
        raise ValueError('Solid end fixture must retain a nonempty modal interior')
    reference = np.linspace(0, 0.2, 1025)
    lo, hi = reference[end_layers], reference[-tip_layers-1]
    planes = np.concatenate([reference[:end_layers+1], np.linspace(lo, hi, interior_elements+1)[1:-1],
                             reference[-tip_layers-1:]])
    write(path, nx=4, ny=8, point_support=True, axial_planes=planes, solid_interfaces=(lo, hi))


if __name__ == '__main__':
    write(Path(__file__).with_name('plate.e'))
    write(Path(__file__).with_name('plate_refined.e'), nx=4, ny=8, nz=128, graded=True)
    write(Path(__file__).with_name('plate_local.e'), nx=4, ny=8, nz=128, graded=True,
          local_interfaces=True, coarsen_interior=True)
    write(Path(__file__).with_name('plate_reference.e'), nx=4, ny=8, nz=1024)
    write(Path(__file__).with_name('plate_reference_check.e'), nx=4, ny=8, nz=512)
    write_solid_ends(Path(__file__).with_name('plate_solid_ends.e'), 32)
    write_solid_ends(Path(__file__).with_name('plate_solid_ends_128_32.e'), 128, tip_layers=32)
