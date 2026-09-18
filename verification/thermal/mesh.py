"""Generate only the four fixed meshes; production input cards are maintained explicitly."""
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent

def mesh(kind):
    axis = kind.startswith('dcax')
    quadratic = kind in ('dcax8', 'dc3d20')
    corners = [(0,0,0),(1,0,0),(1,1,0),(0,1,0)]
    if not axis:
        corners += [(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    edges = [(0,1),(1,2),(2,3),(3,0)]
    if not axis:
        edges += [(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]
    local = list(corners)
    if quadratic:
        local += [tuple((corners[a][d]+corners[b][d])/2 for d in range(3)) for a,b in edges]
    points, elements, lookup = [], [], {}
    for ix in range(4):
        conn = []
        for x,y,z in local:
            key = (round((ix+x)*0.005+(0.01 if axis else 0),10), round(y*0.01,10), round(z*0.01,10))
            if key not in lookup:
                lookup[key] = len(points)+1
                points.append(key)
            conn.append(lookup[key])
        elements.append(conn)
    with netcdf_file(str(ROOT/(kind+'.e')), 'w') as f:
        f.api_version=np.float32(7.22)
        f.version=np.float32(7.22)
        f.floating_point_word_size=np.int32(8)
        f.file_size=np.int32(1)
        f.title=kind+' fixed thermal mesh'
        for name,size in [('time_step',None),('len_name',33),('num_dim',2 if axis else 3),('num_nodes',len(points)),('num_elem',4),('num_el_blk',1),('num_el_in_blk1',4),('num_nod_per_el1',len(local)),('num_side_sets',2),('num_side_ss1',1),('num_side_ss2',1)]:
            f.createDimension(name,size)
        def var(name,dims,values,kind='i'):
            v=f.createVariable(name,kind,dims);v[:]=values;return v
        def names(name,dim,values):
            out=np.zeros((len(values),33),dtype='S1')
            for i,value in enumerate(values):out[i,:len(value)]=np.frombuffer(value.encode(),dtype='S1')
            var(name,(dim,'len_name'),out,'c')
        for d in range(2 if axis else 3):var('coord'+'xyz'[d],('num_nodes',),[p[d] for p in points],'d')
        var('connect1',('num_el_in_blk1','num_nod_per_el1'),elements).elem_type=('QUAD' if axis else 'HEX')+str(len(local))
        var('eb_prop1',('num_el_blk',),[1]).name='ID'
        var('eb_status',('num_el_blk',),[1]);names('eb_names','num_el_blk',['solid'])
        var('ss_prop1',('num_side_sets',),[1,2]).name='ID'
        var('ss_status',('num_side_sets',),[1,1]);names('ss_names','num_side_sets',['left','right'])
        var('elem_ss1',('num_side_ss1',),[1]);var('side_ss1',('num_side_ss1',),[4])
        var('elem_ss2',('num_side_ss2',),[4]);var('side_ss2',('num_side_ss2',),[2])
    # Mesh-only Abaqus include; material, boundaries, and execution are explicit .inp files.
    with (ROOT/(kind+'_mesh.inp')).open('w') as f:
        f.write('*NODE, NSET=ALL\n')
        for i,p in enumerate(points,1):f.write(str(i)+', '+', '.join(format(v,'.12g') for v in p[:2 if axis else 3])+'\n')
        f.write('*ELEMENT, TYPE='+kind.upper()+', ELSET=SOLID\n')
        for i,c in enumerate(elements,1):
            if len(c) == 20:
                c = c[:12]+c[16:20]+c[12:16]
            if len(c)>15:
                f.write(str(i)+', '+', '.join(map(str,c[:15]))+',\n'+', '.join(map(str,c[15:]))+'\n')
            else:f.write(str(i)+', '+', '.join(map(str,c))+'\n')
        f.write('*NSET, NSET=LEFT\n'+', '.join(str(i) for i,p in enumerate(points,1) if p[0]==points[0][0])+'\n')
if __name__ == '__main__':
    for kind in ['dcax4','dcax8','dc3d8','dc3d20']:mesh(kind)

def contact_mesh():
    points=[]
    for r0,r1 in [(0.01,0.02),(0.02,0.03)]:
        points += [(r0,0),(r1,0),(r1,0.01),(r0,0.01),((r0+r1)/2,0),(r1,0.005),((r0+r1)/2,0.01),(r0,0.005)]
    with netcdf_file(str(ROOT/'dcax8_contact.e'),'w') as f:
        f.api_version=np.float32(7.22);f.version=np.float32(7.22)
        f.floating_point_word_size=np.int32(8);f.file_size=np.int32(1);f.title='two independent thermal rings'
        for name,size in [('time_step',None),('len_name',33),('num_dim',2),('num_nodes',16),('num_elem',2),('num_el_blk',2),('num_side_sets',4),('num_el_in_blk1',1),('num_el_in_blk2',1),('num_nod_per_el1',8),('num_nod_per_el2',8)]:f.createDimension(name,size)
        def var(name,dims,values,kind='i'):
            v=f.createVariable(name,kind,dims);v[:]=values;return v
        def names(name,dim,values):
            out=np.zeros((len(values),33),dtype='S1')
            for i,value in enumerate(values):out[i,:len(value)]=np.frombuffer(value.encode(),dtype='S1')
            var(name,(dim,'len_name'),out,'c')
        var('coordx',('num_nodes',),[p[0] for p in points],'d');var('coordy',('num_nodes',),[p[1] for p in points],'d')
        for b in [1,2]:var('connect'+str(b),('num_el_in_blk'+str(b),'num_nod_per_el'+str(b)),[list(range(1+8*(b-1),9+8*(b-1)))]).elem_type='QUAD8'
        var('eb_prop1',('num_el_blk',),[1,2]).name='ID';var('eb_status',('num_el_blk',),[1,1]);names('eb_names','num_el_blk',['inner','outer'])
        var('ss_prop1',('num_side_sets',),[1,2,3,4]).name='ID';var('ss_status',('num_side_sets',),[1]*4);names('ss_names','num_side_sets',['left','right','secondary','primary'])
        for b,(element,side) in enumerate([(1,4),(2,2),(1,2),(2,4)],1):
            dim='num_side_ss'+str(b);f.createDimension(dim,1);var('elem_ss'+str(b),(dim,),[element]);var('side_ss'+str(b),(dim,),[side])
    with (ROOT/'dcax8_contact_mesh.inp').open('w') as f:
        f.write('*NODE, NSET=ALL\n')
        for n,p in enumerate(points,1):f.write(str(n)+', '+', '.join(map(str,p))+'\n')
        f.write('*ELEMENT, TYPE=DCAX8, ELSET=SOLID\n1,1,2,3,4,5,6,7,8\n2,9,10,11,12,13,14,15,16\n')
        f.write('*NSET, NSET=LEFT\n1,4,8\n*NSET, NSET=RIGHT\n10,11,14\n')
if __name__ == '__main__':contact_mesh()

def capacity_mesh(kind):
    axis = kind == 'dcax4'
    points = ([(.01,0),(.025,.002),(.023,.013),(.008,.01)] if axis else
              [(0,0,0),(.012,0,0),(.011,.01,.001),(0,.009,0),
               (.001,0,.01),(.011,.001,.012),(.013,.011,.01),(0,.01,.009)])
    count=len(points)
    name=kind+'_capacity'
    with netcdf_file(str(ROOT/(name+'.e')),'w') as f:
        f.api_version=np.float32(7.22);f.version=np.float32(7.22)
        f.floating_point_word_size=np.int32(8);f.file_size=np.int32(1)
        f.title='distorted thermal capacity probe'
        for key,size in [('time_step',None),('len_name',33),('num_dim',len(points[0])),('num_nodes',count),('num_elem',1),('num_el_blk',1),('num_el_in_blk1',1),('num_nod_per_el1',count),('num_node_sets',1),('num_nod_ns1',count)]:f.createDimension(key,size)
        def var(key,dims,values,typ='i'):
            v=f.createVariable(key,typ,dims);v[:]=values;return v
        def names(key,dim,value):
            a=np.zeros((1,33),dtype='S1');a[0,:len(value)]=np.frombuffer(value.encode(),dtype='S1');var(key,(dim,'len_name'),a,'c')
        for d in range(len(points[0])):var('coord'+'xyz'[d],('num_nodes',),[p[d] for p in points],'d')
        var('connect1',('num_el_in_blk1','num_nod_per_el1'),[list(range(1,count+1))]).elem_type='QUAD4' if axis else 'HEX8'
        var('eb_prop1',('num_el_blk',),[1]).name='ID';var('eb_status',('num_el_blk',),[1]);names('eb_names','num_el_blk','solid')
        var('ns_prop1',('num_node_sets',),[1]).name='ID';var('ns_status',('num_node_sets',),[1]);names('ns_names','num_node_sets','all')
        var('node_ns1',('num_nod_ns1',),list(range(1,count+1)))
    with (ROOT/(name+'_mesh.inp')).open('w') as f:
        f.write('*NODE, NSET=ALL\n')
        for n,p in enumerate(points,1):f.write(str(n)+', '+', '.join(map(str,p))+'\n')
        f.write('*ELEMENT, TYPE='+kind.upper()+', ELSET=SOLID\n1, '+', '.join(map(str,range(1,count+1)))+'\n')

if __name__ == '__main__':
    capacity_mesh('dcax4')
    capacity_mesh('dc3d8')

def distorted_mesh(kind):
    """Curve the quadratic meshes without changing source nodes or connectivity."""
    import shutil
    axis = kind == 'dcax8'
    target = ROOT/(kind+'_distorted.e')
    shutil.copyfile(ROOT/(kind+'.e'), target)
    with netcdf_file(str(target), 'a', mmap=False) as f:
        x = f.variables['coordx'][:].copy()
        y = f.variables['coordy'][:].copy()
        z = np.zeros_like(x) if axis else f.variables['coordz'][:].copy()
        u = (x-(.01 if axis else 0))/.02
        v, w = y/.01, z/.01
        x += .00035*np.sin(np.pi*u)*(v-.3)
        y += .0006*u*(.3+v)
        if not axis:
            x += .0002*np.sin(np.pi*u)*(w-.2)
            y += .0002*w*u
            z += .0004*u*(.2+w)+.0002*v*u
            f.variables['coordz'][:] = z
        f.variables['coordx'][:] = x
        f.variables['coordy'][:] = y
    lines = (ROOT/(kind+'_mesh.inp')).read_text().splitlines()
    with (ROOT/(kind+'_distorted_mesh.inp')).open('w') as stream:
        in_nodes = False
        for line in lines:
            if line.startswith('*'):
                in_nodes = line.upper().startswith('*NODE,')
                stream.write(line+'\n')
            elif in_nodes:
                node = int(line.split(',')[0])-1
                position = [x[node], y[node]] if axis else [x[node], y[node], z[node]]
                stream.write(str(node+1)+', '+', '.join(format(c,'.17g') for c in position)+'\n')
            else:
                stream.write(line+'\n')

if __name__ == '__main__':
    distorted_mesh('dcax8')
    distorted_mesh('dc3d20')

def contact_meshes():
    """Independent-node thermal contact meshes: matched DCAX4, nonmatching DCAX4,
    DC3D8, curved DC3D20, and an open-clearance DCAX8 gas gap."""
    quad4_sides = {1:(0,1),2:(1,2),3:(2,3),4:(3,0)}
    quad8_sides = {1:(0,4,1),2:(1,5,2),3:(2,6,3),4:(3,7,0)}
    hex8_sides = {1:(0,1,5,4),2:(1,2,6,5),3:(2,3,7,6),4:(0,4,7,3),5:(0,3,2,1),6:(4,5,6,7)}
    hex20_sides = {1:(0,1,5,4,8,13,16,12),2:(1,2,6,5,9,14,17,13),3:(2,3,7,6,10,15,18,14),
                   4:(3,0,4,7,11,12,19,15),5:(0,3,2,1,11,10,9,8),6:(4,5,6,7,16,17,18,19)}

    def write(name, axis, exodus_type, abaqus_type, nodes, blocks, block_names, sides, table):
        dim = 2 if axis else 3
        dims = {'time_step':None,'len_name':33,'num_dim':dim,'num_nodes':len(nodes),'num_elem':sum(len(b) for b in blocks),
                'num_el_blk':len(blocks),'num_side_sets':len(sides)}
        for b,block in enumerate(blocks,1):
            dims['num_el_in_blk'+str(b)] = len(block)
            dims['num_nod_per_el'+str(b)] = len(block[0])
        for k,(_,entries) in enumerate(sides,1):
            dims['num_side_ss'+str(k)] = len(entries)
        with netcdf_file(str(ROOT/(name+'.e')),'w') as f:
            f.api_version=np.float32(7.22);f.version=np.float32(7.22)
            f.floating_point_word_size=np.int32(8);f.file_size=np.int32(1);f.title=name+' thermal contact'
            for key,size in dims.items():f.createDimension(key,size)
            def var(key,dims_,values,kind='i'):
                v=f.createVariable(key,kind,dims_);v[:]=values;return v
            def names(key,dim_,values):
                out=np.zeros((len(values),33),dtype='S1')
                for i,value in enumerate(values):out[i,:len(value)]=np.frombuffer(value.encode(),dtype='S1')
                var(key,(dim_,'len_name'),out,'c')
            for d in range(dim):var('coord'+'xyz'[d],('num_nodes',),[p[d] for p in nodes],'d')
            for b,block in enumerate(blocks,1):
                var('connect'+str(b),('num_el_in_blk'+str(b),'num_nod_per_el'+str(b)),block).elem_type=exodus_type
            var('eb_prop1',('num_el_blk',),list(range(1,len(blocks)+1))).name='ID'
            var('eb_status',('num_el_blk',),[1]*len(blocks));names('eb_names','num_el_blk',list(block_names))
            var('ss_prop1',('num_side_sets',),list(range(1,len(sides)+1))).name='ID'
            var('ss_status',('num_side_sets',),[1]*len(sides));names('ss_names','num_side_sets',[s[0] for s in sides])
            for k,(_,entries) in enumerate(sides,1):
                var('elem_ss'+str(k),('num_side_ss'+str(k),),[e for e,_ in entries])
                var('side_ss'+str(k),('num_side_ss'+str(k),),[s for _,s in entries])
        with (ROOT/(name+'_mesh.inp')).open('w') as f:
            f.write('*NODE, NSET=ALL\n')
            for n,p in enumerate(nodes,1):
                f.write(str(n)+', '+', '.join(format(v,'.17g') for v in p[:dim])+'\n')
            f.write('*ELEMENT, TYPE='+abaqus_type+', ELSET=SOLID\n')
            first = 1
            for block in blocks:
                for conn in block:
                    c = list(conn)
                    if len(c) == 20:c = c[:12]+c[16:20]+c[12:16]
                    if len(c)>15:
                        f.write(str(first)+', '+', '.join(map(str,c[:15]))+',\n'+', '.join(map(str,c[15:]))+'\n')
                    else:f.write(str(first)+', '+', '.join(map(str,c))+'\n')
                    first += 1
            for setname,filename in [('left','LEFT'),('right','RIGHT')]:
                elements = [conn for group in blocks for conn in group]
                picked = sorted({elements[entry[0]-1][node] for side_name,entries in sides if side_name == setname
                                 for entry in entries for node in table[entry[1]]})
                f.write('*NSET, NSET='+filename+'\n'+', '.join(map(str,picked))+'\n')

    def ring8(r0,r1):
        return [(r0,0),(r1,0),(r1,0.01),(r0,0.01),((r0+r1)/2,0),(r1,0.005),((r0+r1)/2,0.01),(r0,0.005)]
    interface_sides = [('left',[(1,4)]),('right',[(2,2)]),('secondary',[(1,2)]),('primary',[(2,4)])]

    write('dcax4_contact',True,'QUAD4','DCAX4',
          [(0.01,0),(0.02,0),(0.02,0.01),(0.01,0.01),(0.02,0),(0.03,0),(0.03,0.01),(0.02,0.01)],
          [[[1,2,3,4]],[[5,6,7,8]]],('inner','outer'),interface_sides,quad4_sides)

    write('dcax4_contact_nonmatching',True,'QUAD4','DCAX4',
          [(0.01,0),(0.02,0),(0.02,0.005),(0.01,0.005),(0.02,0.01),(0.01,0.01),
           (0.03,0),(0.03,0.01),(0.02,0.01),(0.02,0)],
          [[[1,2,3,4],[4,3,5,6]],[[10,7,8,9]]],('inner','outer'),
          [('left',[(1,4),(2,4)]),('right',[(3,2)]),('secondary',[(1,2),(2,2)]),('primary',[(3,4)])],quad4_sides)

    cube = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    def block(ox,scale):return [tuple((c[0]*scale[0]+ox,c[1]*scale[1],c[2]*scale[2])) for c in cube]
    write('dc3d8_contact',False,'HEX8','DC3D8',block(0,(.01,.01,.01))+block(.01,(.01,.01,.01)),
          [[list(range(1,9))],[list(range(9,17))]],('hot','cold'),interface_sides,hex8_sides)

    hex_edges = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]
    def curved_block(ox):
        corners = block(ox,(.01,.01,.01))
        local = corners+[tuple((corners[a][d]+corners[b][d])/2 for d in range(3)) for a,b in hex_edges]
        def bend(p):
            u,v,w = p[0]/0.02,p[1]/0.01,p[2]/0.01
            return (p[0]+.00035*np.sin(np.pi*u)*(v-.3)+.0002*np.sin(np.pi*u)*(w-.2),
                    p[1]+.0006*u*(.3+v)+.0002*w*u,
                    p[2]+.0004*u*(.2+w)+.0002*v*u)
        return [bend(p) for p in local]
    write('dc3d20_contact',False,'HEX20','DC3D20',curved_block(0)+curved_block(.01),
          [[list(range(1,21))],[list(range(21,41))]],('hot','cold'),interface_sides,hex20_sides)

    write('dcax8_contact_gap',True,'QUAD8','DCAX8',ring8(0.01,0.0195)+ring8(0.02,0.03),
          [[list(range(1,9))],[list(range(9,17))]],('inner','outer'),interface_sides,quad8_sides)

if __name__ == '__main__':
    contact_meshes()
