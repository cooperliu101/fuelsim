"""Generate only the mesh and native reference deck; never write production cards."""
from pathlib import Path
import shutil
import netCDF4

ROOT = Path(__file__).resolve().parent


def create():
    shutil.copyfile(ROOT/'contact_law.e', ROOT/'contact_small_gap.e')
    with netCDF4.Dataset(ROOT/'contact_small_gap.e', 'r+') as mesh:
        mesh['coordy'][8:16] -= .000099
        coordinates = list(zip(mesh['coordx'][:], mesh['coordy'][:]))
    if len(coordinates) != 16:
        raise ValueError('The reference mesh must contain exactly two independent QUAD8 elements')
    source = (ROOT/'thermal_pressure.inp').read_text()
    start = source.index('*Node\n')+len('*Node\n')
    end = source.index('*Element,', start)
    nodes = ''.join(f'{i+1},{x:.17g},{y:.17g}\n' for i, (x,y) in enumerate(coordinates))
    source = source[:start]+nodes+'17,0,0\n18,0,0\n'+source[end:]
    motion = '0,0,.5,-.0003,1,0'
    if source.count(motion) != 1:
        raise ValueError('Expected one original pressure-contact motion amplitude')
    source = source.replace(motion, '0,0,.25,0,.5,-.000011,.75,0\n1,0')
    (ROOT/'thermal_pressure_small_gap.inp').write_text(source)


if __name__ == '__main__':
    create()
