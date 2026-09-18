"""Generate uniform contact mesh and native decks; input cards are independent."""
from pathlib import Path
import netCDF4
from create_extended_meshes import mesh

ROOT = Path(__file__).resolve().parent


def create():
    with netCDF4.Dataset(ROOT/'contact.e') as d:
        coordinates = list(zip(d['coordx'][:], d['coordy'][:]))
    mesh(ROOT/'contact_law.e', coordinates,
         {'lower': [list(range(1, 9))], 'upper': [list(range(9, 17))]},
         {'field': list(range(1, 17)), 'bottom': [1, 2, 5], 'top': [11, 12, 15],
          'lower': list(range(1, 9)), 'upper': list(range(9, 17))},
         {'primary': [(1, 3)], 'secondary': [(2, 1)]})
    base = (ROOT/'contact.inp').read_text().split('*Step,')[0]
    base = base.replace('*Boundary\nFIELD,1,1,0\nBOTTOM,2,2,0\nCONTROLS,3,5,0\n', '')
    base += '*Nset,nset=LOWER_NODES,generate\n1,8,1\n*Nset,nset=UPPER_NODES,generate\n9,16,1\n'
    for case in ['thermal_pressure', 'thermal_clearance']:
        pressure = case == 'thermal_pressure'
        table = ('*Gap Conductance,pressure\n1000,0,300\n3000,1e6,300\n1200,0,400\n3200,1e6,400'
                 if pressure else '*Gap Conductance\n1000,0,300\n500,.0005,300\n1200,0,400\n700,.0005,400')
        deck = base.replace('*Gap Conductance\n1000,0\n1000,1', table)
        deck += '\n'.join(['*Amplitude,name=MOTION,time=TOTAL TIME',
                           '0,0,.5,-.0003,1,0' if pressure else '0,0,.5,.0003,1,0',
                           '*Amplitude,name=HOT,time=TOTAL TIME', '0,300,1,400',
                           '*Step,name=HISTORY,nlgeom=NO,inc=100',
                           '*Coupled Temperature-Displacement,steady state', '.125,1,.125,.125',
                           '*Boundary', 'FIELD,1,1,0', 'LOWER_NODES,2,2,0', 'BOTTOM,11,11,300',
                           'CONTROLS,3,5,0', '*Boundary,amplitude=MOTION', 'UPPER_NODES,2,2,1',
                           '*Boundary,amplitude=HOT', 'TOP,11,11,1', '*Output,field,frequency=1',
                           '*Node Output', 'NT,U,RF,RFL', '*Element Output', 'S,EE,COORD',
                           '*Contact Output', 'CDISP,CFORCE', '*End Step'])+'\n'
        (ROOT/(case+'.inp')).write_text(deck)
    coupled = base.replace('*Elastic\n1e6,.25', '*Elastic\n1e6,.25\n*Plastic\n500,0\n1500,.1\n'
                           '*Creep,law=TIME\n1e-7,1,0\n*Expansion,zero=300\n1e-5')
    coupled += '\n'.join(['*Amplitude,name=PRESS,time=TOTAL TIME',
                          '0,0,.25,-.0003,.5,-.0005,.75,0', '1,-.0007',
                          '*Amplitude,name=HOT,time=TOTAL TIME', '0,300,1,400',
                          '*Amplitude,name=EXTEND,time=TOTAL TIME', '0,0,1,.003',
                          '*Step,name=HISTORY,nlgeom=NO,inc=1000',
                          '*Coupled Temperature-Displacement', '.0625,1,.0625,.0625',
                          '*Controls,parameters=FIELD,field=DISPLACEMENT', '1e-12,1e-12,,,,1e-12',
                          '*Controls,parameters=FIELD,field=TEMPERATURE', '1e-12,1e-12,,,,1e-12',
                          '*Boundary', 'FIELD,1,1,0', 'BOTTOM,2,2,0', 'BOTTOM,11,11,300',
                          'CONTROLS,4,5,0', '*Boundary,amplitude=EXTEND', 'CONTROLS,3,3,1',
                          '*Boundary,amplitude=PRESS', 'TOP,2,2,1', '*Boundary,amplitude=HOT',
                          'TOP,11,11,1', '*Output,field,frequency=1', '*Node Output', 'NT,U,RF,RFL',
                          '*Element Output', 'S,EE,PE,CE,PEEQ,CEEQ,COORD',
                          '*Contact Output', 'CDISP,CFORCE', '*Output,history,frequency=1',
                          '*Node Output,nset=CONTROLS', 'U3,UR1,UR2,RF3,RM1,RM2', '*End Step'])+'\n'
    (ROOT/'coupled_contact.inp').write_text(coupled)
    # With both rotations fixed to zero this changes only the moment origin,
    # not the imposed thickness or any physical boundary condition. Direct
    # centroid moments avoid subtracting single-precision history quantities.
    (ROOT/'coupled_contact_centered.inp').write_text(
        coupled.replace('17,0,0\n18,0,0','17,0,-.005\n18,0,.0051'))


if __name__ == '__main__':
    create()
