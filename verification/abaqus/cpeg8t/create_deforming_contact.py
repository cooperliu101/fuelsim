"""Generate mesh and native reference deck, never a production input card."""
from pathlib import Path

from create_extended_meshes import mesh

ROOT = Path(__file__).resolve().parent


def create():
    coordinates, blocks, nodesets = [], {}, {'bottom': [], 'top': []}
    for label, bounds in [('lower', [-.03, 0, .03]), ('upper', [-.018, -.006, .006])]:
        node_map, elements = {}, []
        for left, right in zip(bounds[:-1], bounds[1:]):
            points = [(left, 0), (right, 0), (right, 1), (left, 1),
                      ((left+right)/2, 0), (right, .5), ((left+right)/2, 1), (left, .5)]
            element = []
            for x, level in points:
                if (x, level) not in node_map:
                    node = len(coordinates)+1
                    node_map[x, level] = node
                    y = -.01+.01*level if label == 'lower' else .0001+2*(x+.006)**2+.01*level
                    coordinates.append((x, y))
                    if label == 'lower' and level == 0:
                        nodesets['bottom'].append(node)
                    if label == 'upper' and level == 1:
                        nodesets['top'].append(node)
                element.append(node_map[x, level])
            elements.append(element)
        blocks[label] = elements
        nodesets[label] = list(node_map.values())
    nodesets['field'] = list(range(1, len(coordinates)+1))
    mesh(ROOT/'deforming_contact.e', coordinates, blocks, nodesets,
         {'primary': [(1, 3), (2, 3)], 'secondary': [(3, 1), (4, 1)]})
    lines = ['*Heading', '** Deformable curved contact crossing a primary edge boundary.',
             '*Preprint,echo=NO,model=NO,history=NO,contact=NO', '*Node']
    lines += [f'{n},{x:.17g},{y:.17g}' for n, (x, y) in enumerate(coordinates, 1)]
    lines += ['101,0,0', '102,0,0']
    number = 0
    for label, elements in blocks.items():
        lines.append(f'*Element,type=CPEG8T,elset={label}')
        for element in elements:
            number += 1
            lines.append(','.join(map(str, [number]+element)))
    for label, nodes in nodesets.items():
        lines += [f'*Nset,nset={label}', ','.join(map(str, nodes))]
    lines += ['*Nset,nset=CONTROLS', '101,102', '*Material,name=SOLID', '*Elastic', '1e6,.25',
              '*Conductivity', '10', '*Density', '1000', '*Specific Heat', '100',
              '*Solid Section,elset=lower,material=SOLID,ref node=101', '.1,0,0',
              '*Solid Section,elset=upper,material=SOLID,ref node=102', '.1,0,0',
              '*Surface,name=PRIMARY,type=ELEMENT', 'lower,S3',
              '*Surface,name=SECONDARY,type=ELEMENT', 'upper,S1',
              '*Surface Interaction,name=INTERFACE', '.1',
              '*Surface Behavior,pressure-overclosure=LINEAR', '1e9', '*Gap Conductance', '1000,0', '1000,1',
              '*Contact Pair,interaction=INTERFACE,type=SURFACE TO SURFACE', 'SECONDARY,PRIMARY',
              '*Initial Conditions,type=TEMPERATURE', 'field,300',
              '*Amplitude,name=SLIDE,time=TOTAL TIME', '0,0,.25,0,1,.012',
              '*Amplitude,name=PRESS,time=TOTAL TIME', '0,0,.25,-.0005,1,-.001',
              '*Amplitude,name=HOT,time=TOTAL TIME', '0,300,1,400',
              '*Step,name=HISTORY,nlgeom=YES,inc=1000', '*Coupled Temperature-Displacement,steady state',
              '.0625,1,.0625,.0625', '*Controls,parameters=FIELD,field=DISPLACEMENT',
              '1e-12,1e-12,,,,1e-12', '*Boundary', 'bottom,1,2,0', 'bottom,11,11,300', 'CONTROLS,3,5,0',
              '*Boundary,amplitude=SLIDE', 'top,1,1,1', '*Boundary,amplitude=PRESS', 'top,2,2,1',
              '*Boundary,amplitude=HOT', 'top,11,11,1', '*Output,field,frequency=1',
              '*Node Output', 'NT,U,RF,RFL', '*Element Output', 'S,EE,COORD',
              '*Contact Output', 'CDISP,CFORCE', '*Output,history,frequency=1',
              '*Node Output,nset=CONTROLS', 'U3,UR1,UR2,RF3,RM1,RM2', '*End Step']
    (ROOT/'deforming_contact.inp').write_text('\n'.join(lines)+'\n')
    fixed = [line.replace('bottom,1,2,0', 'lower,1,2,0') for line in lines]
    (ROOT/'deforming_fixed_primary.inp').write_text('\n'.join(fixed)+'\n')


if __name__ == '__main__':
    create()
