"""Generate only the committed extrusion mesh, never create or edit input cards."""
from pathlib import Path
import numpy as np
from netCDF4 import Dataset


def names(data, variable, dimension, values):
    out = np.zeros((len(values), 64), dtype='S1')
    for i, value in enumerate(values):
        out[i, :len(value)] = np.frombuffer(value.encode(), dtype='S1')
    data.createVariable(variable, 'S1', (dimension, 'len_name'))[:] = out


def write(path, nx=2, ny=4, nz=8):
    # Exodus HEX20 order: corner nodes, lower edges, vertical edges, upper edges.
    offsets = [(0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2),
               (1,0,0),(2,1,0),(1,2,0),(0,1,0),(0,0,1),(2,0,1),(2,2,1),(0,2,1),
               (1,0,2),(2,1,2),(1,2,2),(0,1,2)]
    ids, nodes, elements = {}, [], []
    sides = {key: [] for key in ('end', 'end_plus', 'end_minus', 'broad', 'half_broad')}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cell = []
                for a,b,c in offsets:
                    key = (2*i+a, 2*j+b, 2*k+c)
                    if key not in ids:
                        ids[key] = len(nodes)+1
                        nodes.append((-0.001+0.002*key[0]/(2*nx), -0.01+0.02*key[1]/(2*ny), 0.2*key[2]/(2*nz)))
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

if __name__ == '__main__':
    write(Path(__file__).with_name('plate.e'))
    write(Path(__file__).with_name('plate_refined.e'), nx=4, ny=8, nz=32)
