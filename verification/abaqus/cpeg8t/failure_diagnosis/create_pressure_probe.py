"""Prescribed uniform opening separates pressure from thermal activation."""
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GAPS=[1e-4,1e-5,1e-6,1e-7,1e-8,0.,-1e-8,-1e-4,1e-4]


def create():
    prefix=(ROOT.parent/'thermal_pressure.inp').read_text().split('*Amplitude')[0]
    lines=[prefix.rstrip()]
    for i,gap in enumerate(GAPS):
        lines+= [f'*Step,name=GAP_{i},nlgeom=NO,inc=10',
                 '*Coupled Temperature-Displacement,steady state','1,1,1,1',
                 '*Boundary,op=NEW','FIELD,1,1,0','LOWER_NODES,2,2,0',
                 f'UPPER_NODES,2,2,{gap-.0001:.17g}','CONTROLS,3,5,0',
                 'LOWER_NODES,11,11,300','UPPER_NODES,11,11,400',
                 '*Output,field,frequency=1','*Node Output','NT,U,RF,RFL',
                 '*Contact Output','CSTRESS,CDISP,CFORCE','*End Step']
    (ROOT/'pressure_opening.inp').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':create()
