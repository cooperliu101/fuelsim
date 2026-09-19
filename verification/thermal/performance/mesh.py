"""Generate only identical structured DC3D8 geometry for the two solvers."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent


def write(interface):
    name = 'interface' if interface else 'bulk'
    n = 24
    corners = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
               (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    nodes, blocks, sides = [], {}, {'left': [], 'right': []}
    groups = {'left': [], 'right': []}
    if interface:
        sides.update(secondary=[], primary=[])
    eid = 0
    for b in range(2 if interface else 1):
        nx = n//2 if interface else n
        offset = b*nx
        lookup = {}
        for i in range(nx+1):
            for j in range(n+1):
                for k in range(n+1):
                    nodes.append((.02*(offset+i)/n, .02*j/n, .02*k/n))
                    lookup[i,j,k] = len(nodes)
                    if offset+i == 0:
                        groups['left'].append(len(nodes))
                    if offset+i == n:
                        groups['right'].append(len(nodes))
        cells = []
        for i in range(nx):
            for j in range(n):
                for k in range(n):
                    cells.append([lookup[i+a,j+c,k+d] for a,c,d in corners])
                    eid += 1
                    if offset+i == 0:
                        sides['left'].append((eid,4))
                    if offset+i == n-1:
                        sides['right'].append((eid,2))
                    if interface and b == 0 and i == nx-1:
                        sides['secondary'].append((eid,2))
                    if interface and b == 1 and i == 0:
                        sides['primary'].append((eid,4))
        blocks[('hot','cold')[b] if interface else 'solid'] = cells
    with netcdf_file(str(ROOT/(name+'.e')), 'w') as f:
        f.api_version = np.float32(7.22)
        f.version = np.float32(7.22)
        f.floating_point_word_size = np.int32(8)
        f.file_size = np.int32(1)
        f.title = name
        dims = {'time_step':None, 'len_name':33, 'num_dim':3, 'num_nodes':len(nodes),
                'num_elem':eid, 'num_el_blk':len(blocks), 'num_side_sets':len(sides), 'num_node_sets':2}
        for i,cells in enumerate(blocks.values(),1):
            dims.update({f'num_el_in_blk{i}':len(cells), f'num_nod_per_el{i}':8})
        for i,entries in enumerate(sides.values(),1):
            dims[f'num_side_ss{i}'] = len(entries)
        for i,entries in enumerate(groups.values(),1):
            dims[f'num_nod_ns{i}'] = len(entries)
        for key,value in dims.items():
            f.createDimension(key,value)
        def var(key,dim,values,kind='i'):
            v=f.createVariable(key,kind,dim);v[:]=values
            return v
        for d in range(3):
            var('coord'+'xyz'[d],('num_nodes',),[p[d] for p in nodes],'d')
        for i,cells in enumerate(blocks.values(),1):
            var(f'connect{i}',(f'num_el_in_blk{i}',f'num_nod_per_el{i}'),cells).elem_type='HEX8'
        for prefix,dim,keys in [('eb','num_el_blk',blocks),('ss','num_side_sets',sides),('ns','num_node_sets',groups)]:
            var(prefix+'_prop1',(dim,),range(1,len(keys)+1)).name='ID'
            var(prefix+'_status',(dim,),[1]*len(keys))
            names=np.zeros((len(keys),33),dtype='S1')
            for i,key in enumerate(keys):
                names[i,:len(key)]=np.frombuffer(key.encode(),dtype='S1')
            var(prefix+'_names',(dim,'len_name'),names,'c')
        for i,entries in enumerate(sides.values(),1):
            var(f'elem_ss{i}',(f'num_side_ss{i}',),[e for e,s in entries])
            var(f'side_ss{i}',(f'num_side_ss{i}',),[s for e,s in entries])
        for i,entries in enumerate(groups.values(),1):
            var(f'node_ns{i}',(f'num_nod_ns{i}',),entries)
    with (ROOT/(name+'_mesh.inp')).open('w') as f:
        f.write('*NODE, NSET=ALL\n')
        for i,point in enumerate(nodes,1):
            f.write(str(i)+', '+', '.join(format(x,'.17g') for x in point)+'\n')
        eid=0
        for block,cells in blocks.items():
            f.write('*ELEMENT, TYPE=DC3D8, ELSET='+block.upper()+'\n')
            for cell in cells:
                eid+=1;f.write(str(eid)+', '+', '.join(map(str,cell))+'\n')
        for key,entries in groups.items():
            f.write('*NSET, NSET='+key.upper()+'\n')
            for start in range(0,len(entries),16):
                f.write(', '.join(map(str,entries[start:start+16]))+'\n')
        if interface:
            for key,label in [('secondary','S4'),('primary','S6')]:
                f.write('*SURFACE, NAME='+key.upper()+', TYPE=ELEMENT\n')
                for e,s in sides[key]:
                    f.write(str(e)+', '+label+'\n')
    print(name, 'elements', eid, 'temperature_unknowns', len(nodes))


if __name__ == '__main__':
    write(False)
    write(True)
