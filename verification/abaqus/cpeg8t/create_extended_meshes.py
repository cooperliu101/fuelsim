"""Generate standard Exodus meshes and native decks, never production input cards."""
from pathlib import Path

import netCDF4
import numpy as np

ROOT = Path(__file__).resolve().parent


def names(d, key, dim, values):
    v = d.createVariable(key, 'S1', (dim, 'len_name'))
    v[:] = np.zeros((len(values), 33), dtype='S1')
    for i, value in enumerate(values):
        v[i, :len(value)] = np.array(list(value), dtype='S1')


def mesh(path, coordinates, blocks, nodesets, sidesets):
    with netCDF4.Dataset(path, 'w', format='NETCDF3_64BIT_OFFSET') as d:
        d.title = 'CPEG8T independent validation mesh'
        d.api_version = np.float32(8)
        d.version = np.float32(8)
        d.floating_point_word_size = np.int32(8)
        d.file_size = np.int32(1)
        dims = dict(len_name=33, num_dim=2, num_nodes=len(coordinates),
                    num_elem=sum(len(b) for b in blocks.values()), num_el_blk=len(blocks),
                    time_step=None)
        if nodesets:
            dims['num_node_sets'] = len(nodesets)
        if sidesets:
            dims['num_side_sets'] = len(sidesets)
        for key, size in dims.items():
            d.createDimension(key, size)
        for c, axis in enumerate('xy'):
            d.createVariable('coord'+axis, 'f8', ('num_nodes',))[:] = np.array(coordinates)[:, c]
        names(d, 'coor_names', 'num_dim', ['x', 'y'])
        for prefix, dim, entries in [('eb', 'num_el_blk', blocks), ('ns', 'num_node_sets', nodesets),
                                     ('ss', 'num_side_sets', sidesets)]:
            if not entries:
                continue
            ids = d.createVariable(prefix+'_prop1', 'i4', (dim,))
            ids.setncattr('name', 'ID')
            ids[:] = np.arange(1, len(entries)+1)
            d.createVariable(prefix+'_status', 'i4', (dim,))[:] = 1
            names(d, prefix+'_names', dim, list(entries))
        for i, elements in enumerate(blocks.values(), 1):
            ed, nd = f'num_el_in_blk{i}', f'num_nod_per_el{i}'
            d.createDimension(ed, len(elements))
            d.createDimension(nd, 8)
            v = d.createVariable(f'connect{i}', 'i4', (ed, nd))
            v.elem_type = 'QUAD8'
            v[:] = elements
        for i, nodes in enumerate(nodesets.values(), 1):
            dim = f'num_nod_ns{i}'
            d.createDimension(dim, len(nodes))
            d.createVariable(f'node_ns{i}', 'i4', (dim,))[:] = nodes
        for i, sides in enumerate(sidesets.values(), 1):
            dim = f'num_side_ss{i}'
            d.createDimension(dim, len(sides))
            d.createVariable(f'elem_ss{i}', 'i4', (dim,))[:] = [s[0] for s in sides]
            d.createVariable(f'side_ss{i}', 'i4', (dim,))[:] = [s[1] for s in sides]


def inelastic():
    rectangle = [(-.01,-.005),(.01,-.005),(.01,.005),(-.01,.005),
                 (0,-.005),(.01,0),(0,.005),(-.01,0)]
    labels = ['plastic', 'creep', 'coupled']
    coords = [(x+.03*b, y) for b in range(3) for x,y in rectangle]
    blocks = {name: [list(range(8*b+1, 8*b+9))] for b,name in enumerate(labels)}
    corners = [8*b+i for b in range(3) for i in range(1,5)]
    mesh(ROOT/'inelastic.e', coords, blocks, {'field':list(range(1,25)), 'corners':corners}, {})
    lines = ['*Heading', '** Three independent mechanisms; uniform axial strain history.',
             '*Preprint,echo=NO,model=NO,history=NO,contact=NO', '*Node']
    lines += [f'{i},{x:.16g},{y:.16g}' for i,(x,y) in enumerate(coords,1)]
    lines += [f'{101+b},{.03*b},0' for b in range(3)]
    for b,label in enumerate(labels):
        lines += [f'*Element,type=CPEG8T,elset={label}',
                  ','.join(map(str,[b+1]+blocks[label][0]))]
        lines += [f'*Material,name={label}', '*Elastic','1e6,.25']
        if label != 'creep':
            lines += ['*Plastic','1000,0','2000,.1']
        if label != 'plastic':
            lines += ['*Creep,law=TIME','1e-7,1,0']
        lines += ['*Conductivity','10','*Density','1000','*Specific Heat','100',
                  f'*Solid Section,elset={label},material={label},ref node={101+b}', '.1,0,0']
    lines += ['*Nset,nset=FIELD,generate','1,24,1','*Nset,nset=CONTROLS','101,102,103',
              '*Initial Conditions,type=TEMPERATURE','FIELD,300',
              '*Amplitude,name=EXTENSION,time=TOTAL TIME','0,0,1,.001,2,.001,3,.0002',
              '*Step,name=HISTORY,nlgeom=NO,inc=100', '*Coupled Temperature-Displacement',
              '.25,3,.25,.25', '*Controls,parameters=FIELD,field=DISPLACEMENT',
              '1e-12,1e-12,,,,1e-12', '*Boundary','FIELD,1,2,0','FIELD,11,11,300',
              'CONTROLS,4,5,0','*Boundary,amplitude=EXTENSION','CONTROLS,3,3,1',
              '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL',
              '*Element Output','S,EE,PE,PEEQ,CE,CEEQ,COORD,IVOL',
              '*Output,history,frequency=1','*Node Output,nset=CONTROLS',
              'U3,UR1,UR2,RF3,RM1,RM2','*End Step']
    deck = '\n'.join(lines)+'\n'
    (ROOT/'inelastic_history.inp').write_text(deck)
    (ROOT/'inelastic_finite.inp').write_text(deck.replace('nlgeom=NO', 'nlgeom=YES'))


def distorted():
    # Bilinear trapezoid with straight sides; its Jacobian varies spatially.
    coords = [(-.01,-.005),(.01,-.005),(.014,.005),(-.008,.005),
              (0,-.005),(.012,0),(.003,.005),(-.009,0)]
    mesh(ROOT/'distorted.e', coords, {'solid':[list(range(1,9))]},
         {'field':list(range(1,9)),'bottom':[1,2,5],'top':[3,4,7],'middle':[6,8]},
         {'bottom':[(1,1)],'right':[(1,2)],'top':[(1,3)],'left':[(1,4)]})
    # Exact trapezoid centroid from the polygon area formula.
    polygon = np.array(coords[:4])
    cross = polygon[:,0]*np.roll(polygon[:,1],-1)-np.roll(polygon[:,0],-1)*polygon[:,1]
    centroid = np.sum((polygon+np.roll(polygon,-1,axis=0))*cross[:,None],axis=0)/(3*sum(cross))
    u3 = .001-.002*centroid[1]-.003*centroid[0]
    lines = ['*Heading','** Distorted CPEG8T: shear and temperature-dependent properties.',
             '*Preprint,echo=NO,model=NO,history=NO,contact=NO','*Node']
    lines += [f'{i},{x:.16g},{y:.16g}' for i,(x,y) in enumerate(coords,1)]
    lines += ['100,0,0','*Element,type=CPEG8T,elset=SOLID','1,1,2,3,4,5,6,7,8',
              '*Nset,nset=FIELD,generate','1,8,1','*Nset,nset=BOTTOM','1,2,5',
              '*Nset,nset=TOP','3,4,7','*Nset,nset=MIDDLE','6,8','*Nset,nset=CONTROLS','100',
              '*Material,name=SOLID','*Elastic','1e6,.25,300','1.1e6,.25,400',
              '*Expansion,zero=300','1e-5','*Conductivity','10,300','11,400',
              '*Density','1000','*Specific Heat','100,300','110,400',
              '*Solid Section,elset=SOLID,material=SOLID,ref node=100','.1,0,0',
              '*Initial Conditions,type=TEMPERATURE','FIELD,300',
              '*Step,name=PROBE,nlgeom=NO,inc=1','*Coupled Temperature-Displacement,steady state',
              '1,1,1,1','*Boundary','BOTTOM,1,1,-.00005','TOP,1,1,.00005','MIDDLE,1,1,0',
              'FIELD,2,2,0',f'100,3,3,{u3:.16g}','100,4,4,.002','100,5,5,-.003',
              'BOTTOM,11,11,300','TOP,11,11,400',
              '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL',
              '*Element Output','S,EE,COORD,IVOL,HFL','*Output,history,frequency=1',
              '*Node Output,nset=CONTROLS','U3,UR1,UR2,RF3,RM1,RM2','*End Step']
    (ROOT/'shear_distorted.inp').write_text('\n'.join(lines)+'\n')


def sliding():
    coords=[(-.02,-.01),(0,-.01),(0,0),(-.02,0),(-.01,-.01),(0,-.005),(-.01,0),(-.02,-.005),
            (.02,-.01),(.02,0),(.01,-.01),(.02,-.005),(.01,0),
            (-.01,.0001),(0,.0001),(0,.0101),(-.01,.0101),(-.005,.0001),(0,.0051),
            (-.005,.0101),(-.01,.0051)]
    blocks={'lower':[[1,2,3,4,5,6,7,8],[2,9,10,3,11,12,13,6]],'upper':[list(range(14,22))]}
    mesh(ROOT/'sliding.e',coords,blocks,{'lower':list(range(1,14)),'upper':list(range(14,22))},
         {'primary':[(1,3),(2,3)],'secondary':[(3,1)]})
    lines=['*Heading','** Fully prescribed nonmatching sliding contact with changing section thickness.',
           '*Preprint,echo=NO,model=NO,history=NO,contact=NO','*Node']
    lines += [f'{i},{x:.16g},{y:.16g}' for i,(x,y) in enumerate(coords,1)]
    lines += ['101,0,0','102,0,0']
    for label,els in blocks.items():
        lines += [f'*Element,type=CPEG8T,elset={label}']
        lines += [','.join(map(str,[i+1 if label=='lower' else 3]+e)) for i,e in enumerate(els)]
    lines += ['*Nset,nset=LOWERNODES,generate','1,13,1','*Nset,nset=UPPERNODES,generate','14,21,1',
              '*Nset,nset=FIELD,generate','1,21,1','*Nset,nset=CONTROLS','101,102',
              '*Material,name=SOLID','*Elastic','1e6,.25','*Conductivity','10','*Density','1000',
              '*Specific Heat','100','*Solid Section,elset=lower,material=SOLID,ref node=101','.1,0,0',
              '*Solid Section,elset=upper,material=SOLID,ref node=102','.1,0,0',
              '*Surface,name=PRIMARY,type=ELEMENT','lower,S3','*Surface,name=SECONDARY,type=ELEMENT','upper,S1',
              '*Surface Interaction,name=INTERFACE','.1','*Surface Behavior,pressure-overclosure=LINEAR','1e9',
              '*Gap Conductance','1000,0','1000,1',
              '*Contact Pair,interaction=INTERFACE,type=SURFACE TO SURFACE','SECONDARY,PRIMARY',
              '*Initial Conditions,type=TEMPERATURE','FIELD,300',
              '*Amplitude,name=RAMP,time=TOTAL TIME','0,0,3,1',
              '*Amplitude,name=VERTICAL,time=TOTAL TIME','0,0,.5,-.0002,1,-.0002,1.5,0','2,0,2.5,-.0002,3,-.0002',
              '*Amplitude,name=TEMP,time=TOTAL TIME','0,300,3,400',
              '*Step,name=CYCLE,nlgeom=YES,inc=100','*Coupled Temperature-Displacement','.25,3,.25,.25',
              '*Boundary','LOWERNODES,1,2,0','LOWERNODES,11,11,300','101,3,5,0',
              '*Boundary,amplitude=RAMP','UPPERNODES,1,1,.015','102,3,3,.010048','102,4,4,.02','102,5,5,-.03',
              '*Boundary,amplitude=VERTICAL','UPPERNODES,2,2,1',
              '*Boundary,amplitude=TEMP','UPPERNODES,11,11,1',
              '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL','*Element Output','S,EE,COORD,IVOL,HFL',
              '*Contact Output','CSTRESS,CDISP,CFORCE','*Output,history,frequency=1',
              '*Node Output,nset=CONTROLS','U3,UR1,UR2,RF3,RM1,RM2','*End Step']
    (ROOT/'contact_cycle.inp').write_text('\n'.join(lines)+'\n')


def thermal_boundary():
    coords=[(-.01,-.005),(.01,-.005),(.01,.005),(-.01,.005),
            (0,-.005),(.01,0),(0,.005),(-.01,0)]
    mesh(ROOT/'boundary.e',coords,{'solid':[list(range(1,9))]},
         {'field':list(range(1,9)),'bottom':[1,2,5],'top':[3,4,7],
          **{f'node{i}':[i] for i in range(1,9)}},
         {'bottom':[(1,1)],'right':[(1,2)],'top':[(1,3)],'left':[(1,4)]})
    lines=['*Heading','** Heat flux, convection and volume heating on CPEG8T.',
           '*Preprint,echo=NO,model=NO,history=NO,contact=NO','*Node']
    lines += [f'{i},{x:.16g},{y:.16g}' for i,(x,y) in enumerate(coords,1)]
    lines += ['100,0,0','*Element,type=CPEG8T,elset=SOLID','1,1,2,3,4,5,6,7,8',
              '*Nset,nset=FIELD,generate','1,8,1','*Nset,nset=BOTTOM','1,2,5','*Nset,nset=CONTROLS','100',
              '*Material,name=SOLID','*Elastic','1e6,.25','*Expansion,zero=300','1e-5',
              '*Conductivity','10','*Density','1000','*Specific Heat','100',
              '*Solid Section,elset=SOLID,material=SOLID,ref node=100','.1,0,0',
              '*Initial Conditions,type=TEMPERATURE','FIELD,300',
              '*Step,name=PROBE,nlgeom=NO,inc=1','*Coupled Temperature-Displacement,steady state','1,1,1,1',
              '*Boundary','FIELD,1,2,0','CONTROLS,3,5,0','BOTTOM,11,11,300',
              '*Dflux','SOLID,BF,1e5','SOLID,S4,100','SOLID,S2,50','*Film','SOLID,F3,350,1000',
              '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL','*Element Output','S,EE,COORD,IVOL,HFL',
              '*Output,history,frequency=1','*Node Output,nset=CONTROLS','U3,UR1,UR2,RF3,RM1,RM2','*End Step']
    (ROOT/'thermal_boundary.inp').write_text('\n'.join(lines)+'\n')


def convergence_meshes():
    for count in [2,4,8]:
        coordinates=[]
        index={}
        elements=[]
        for j in range(count):
            for i in range(count):
                nodes=[]
                for di,dj in [(0,0),(2,0),(2,2),(0,2),(1,0),(2,1),(1,2),(0,1)]:
                    key=(2*i+di,2*j+dj)
                    if key not in index:
                        index[key]=len(coordinates)+1
                        coordinates.append((key[0]/(2*count),key[1]/(2*count)))
                    nodes.append(index[key])
                elements.append(nodes)
        mesh(ROOT/f'convergence_{count}.e',coordinates,{'solid':elements},
             {'field':list(range(1,len(coordinates)+1)),
              'bottom':[i for (x,y),i in index.items() if y==0]},
             {'top':[(count*(count-1)+i+1,3) for i in range(count)]})


def curved():
    coordinates = [(-.01, -.01), (.01, -.01), (.01, .001), (-.01, .001),
                   (0., -.01), (.01, -.0045), (0., 0.), (-.01, -.0045),
                   (-.008, .00074), (.008, .00074), (.008, .01074), (-.008, .01074),
                   (0., .0001), (.008, .00574), (0., .0101), (-.008, .00574)]
    mesh(ROOT/'curved.e', coordinates,
         {'lower': [list(range(1, 9))], 'upper': [list(range(9, 17))]},
         {'lower': list(range(1, 9)), 'upper': list(range(9, 17))},
         {'primary': [(1, 3)], 'secondary': [(2, 1)]})


if __name__ == '__main__':
    inelastic()
    distorted()
    sliding()
    thermal_boundary()
    convergence_meshes()
    curved()
