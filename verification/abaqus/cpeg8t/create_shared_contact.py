"""Create a connected quadratic secondary surface; no input cards are generated."""
from pathlib import Path
from create_extended_meshes import mesh

ROOT = Path(__file__).resolve().parent


def create():
    coordinates, blocks, nodesets = [], {}, {}
    for label, bounds in [('lower',[-.02,0,.02]),('upper',[-.01,0,.01])]:
        node_map, elements = {}, []
        for left,right in zip(bounds[:-1],bounds[1:]):
            points = [(left,0),(right,0),(right,1),(left,1),((left+right)/2,0),
                      (right,.5),((left+right)/2,1),(left,.5)]
            element = []
            for x,level in points:
                if (x,level) not in node_map:
                    node_map[x,level] = len(coordinates)+1
                    y = -.01+.01*level if label=='lower' else .0001+.5*x*x+.01*level
                    coordinates.append((x,y))
                element.append(node_map[x,level])
            elements.append(element)
        blocks[label] = elements
        nodesets[label] = list(node_map.values())
    mesh(ROOT/'shared_contact.e',coordinates,blocks,nodesets,
         {'primary':[(1,3),(2,3)],'secondary':[(3,1),(4,1)]})
    lines = ['*Heading','** Curved two-edge secondary surface shares a nonuniform pressure constraint.',
             '*Preprint,echo=NO,model=NO,history=NO,contact=NO','*Node']
    lines += [f'{n},{x:.17g},{y:.17g}' for n,(x,y) in enumerate(coordinates,1)]
    lines += ['101,0,0','102,0,0']
    number = 0
    for label,elements in blocks.items():
        lines.append(f'*Element,type=CPEG8T,elset={label}')
        for element in elements:
            number += 1
            lines.append(','.join(map(str,[number]+element)))
        lines += [f'*Nset,nset={label}',','.join(map(str,nodesets[label]))]
    lines += ['*Nset,nset=FIELD,generate',f'1,{len(coordinates)},1','*Nset,nset=CONTROLS',
              '101,102','*Material,name=SOLID','*Elastic','1e6,.25','*Conductivity','10',
              '*Density','1000','*Specific Heat','100',
              '*Solid Section,elset=lower,material=SOLID,ref node=101','.1,0,0',
              '*Solid Section,elset=upper,material=SOLID,ref node=102','.1,0,0',
              '*Surface,name=PRIMARY,type=ELEMENT','lower,S3',
              '*Surface,name=SECONDARY,type=ELEMENT','upper,S1',
              '*Surface Interaction,name=INTERFACE','.1','*Surface Behavior,pressure-overclosure=LINEAR','1e9',
              '*Contact Pair,interaction=INTERFACE,type=SURFACE TO SURFACE','SECONDARY,PRIMARY',
              '*Initial Conditions,type=TEMPERATURE','FIELD,300',
              '*Step,name=CONTACT,nlgeom=NO,inc=10','*Coupled Temperature-Displacement,steady state',
              '1,1,1,1','*Boundary','FIELD,1,1,0','FIELD,11,11,300','lower,2,2,0',
              'upper,2,2,-.0003','CONTROLS,3,5,0',
              '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL',
              '*Element Output','S,EE,COORD','*Contact Output','CDISP,CFORCE',
              '*End Step']
    (ROOT/'shared_contact.inp').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    create()
