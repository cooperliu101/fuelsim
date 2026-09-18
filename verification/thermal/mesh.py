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
