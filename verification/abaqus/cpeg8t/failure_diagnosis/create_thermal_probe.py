"""Native frozen-geometry thermal matrices; does not generate production cards."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent
CORNERS = [1,2,3,4,9,10,14,15,16,17,22,23]


def create():
    native = json.loads((SOURCE/'deforming_fixed_primary_fields.json').read_text())['steps']['HISTORY']
    prefix = (SOURCE/'deforming_fixed_primary.inp').read_text().split('*Amplitude')[0]
    for enabled in (True,False):
        lines = [prefix if enabled else prefix.replace('*Gap Conductance\n1000,0\n1000,1\n',
                                                      '*Gap Conductance\n0,0\n0,1\n')]
        for at in (1,13):
            frame = native[at]
            displacements = {v['nodeLabel']:v['data'] for v in frame['fields']['U']['values']}
            for column in [-1]+list(range(len(CORNERS))):
                lines += [f'*Step,name=F{at}_C{column+1},nlgeom=YES,inc=10',
                          '*Coupled Temperature-Displacement,steady state','1,1,1,1',
                          '*Boundary,op=NEW','CONTROLS,3,5,0']
                for node in range(1,27):
                    for component in range(2):
                        lines.append(f'{node},{component+1},{component+1},{displacements[node][component]:.17g}')
                for i,node in enumerate(CORNERS):
                    lines.append(f'{node},11,11,{301 if i==column else 300}')
                lines += ['*Output,field,frequency=1','*Node Output','NT,U,RF,RFL','*End Step']
        name = 'thermal_matrix_on' if enabled else 'thermal_matrix_off'
        (ROOT/(name+'.inp')).write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    create()
