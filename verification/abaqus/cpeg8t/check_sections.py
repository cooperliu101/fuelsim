"""Check two independent sections against small-strain closed-form elasticity."""
import argparse
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--executable', type=Path, required=True)
parser.add_argument('--work', type=Path, required=True)
args = parser.parse_args()
source = Path(__file__).resolve().parent
args.work.mkdir(parents=True, exist_ok=True)
for name in ('sections.fsi', 'contact.e'):
    shutil.copyfile(source / name, args.work / name)
    assert (source / name).read_bytes() == (args.work / name).read_bytes()
process = subprocess.run([str(args.executable.resolve()), '-i', 'sections.fsi'],
                         cwd=args.work, capture_output=True, text=True)
(args.work / 'production.log').write_text(process.stdout + process.stderr)
if process.returncode or 'completed=true' not in process.stdout:
    raise RuntimeError(process.stdout + process.stderr)


def check(actual, expected, label, atol=1e-10):
    if not np.allclose(actual, expected, rtol=1e-10, atol=atol):
        raise AssertionError(f'{label}: {actual} != {expected}')


with netCDF4.Dataset(args.work / 'sections_results.e') as output:
    times = np.asarray(output['time_whole'][:])
    check(times, [0, .5, 1], 'Initial state and both time steps must be exported')
    global_names = list(netCDF4.chartostring(output['name_glo_var'][:]))
    element_names = list(netCDF4.chartostring(output['name_elem_var'][:]))
    nodal_names = list(netCDF4.chartostring(output['name_nod_var'][:]))
    for frame, time in enumerate(times):
        globals_ = dict(zip(global_names, output['vals_glo_var'][frame]))
        for block, section, height, center, controls in [
            (1, 'lower', .1, -.005, np.array([.001, .02, -.03]) * time),
            (2, 'upper', .2, .0051, np.array([0, -.04, .05]) * max(0, 2*time-1)),
        ]:
            prefix = 'section_' + section + '_'
            for name, value in zip(('u3', 'rotation_x', 'rotation_y'), controls):
                check(globals_[prefix + name], value, prefix + name, 1e-13)
            for name, value in [('origin_x', 0), ('origin_y', center), ('initial_thickness', height)]:
                check(globals_[prefix + name], value, prefix + name, 1e-13)
            u3, rx, ry = controls
            # The constrained in-plane strain is zero; C3333=1.2e6 Pa.
            for name, value in [('axial_force', 1.2e6*.02*.01*u3/height),
                                ('moment_x', 1.2e6*.02*.01**3/12*rx/height),
                                ('moment_y', 1.2e6*.01*.02**3/12*ry/height)]:
                check(globals_[prefix + name], value, prefix + name)
            # CPEG8T ordering is eta outside, xi inside.
            gauss = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
            for j, eta in enumerate(gauss):
                for i, xi in enumerate(gauss):
                    strain = (u3 + rx*.005*eta - ry*.01*xi) / height
                    for component, modulus in [('xx', .4e6), ('yy', .4e6), ('zz', 1.2e6), ('xy', 0)]:
                        name = f'stress_{component}_q{3*j+i}'
                        variable = element_names.index(name)+1
                        check(output[f'vals_elem_var{variable}eb{block}'][frame, 0], modulus*strain,
                              section+' '+name)
        for name, expected in [('temperature', 300), ('displacement_x', 0), ('displacement_y', 0)]:
            values = np.asarray(output[f'vals_nod_var{nodal_names.index(name)+1}'][frame])
            valid = np.isfinite(values)
            expected_nodes = np.ones(16, dtype=bool)
            if name == 'temperature':
                expected_nodes[4:8] = False
                expected_nodes[12:16] = False
            if not np.array_equal(valid, expected_nodes):
                raise AssertionError('Incorrect field node coverage: ' + name)
            check(values[np.isfinite(values)], expected, name)
print('Independent section values, time functions, free extension, reactions and all 18 material points passed')
