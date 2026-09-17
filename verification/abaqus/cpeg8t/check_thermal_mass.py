"""Independent uniform-temperature energy balance and global mechanical balance."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np


def check(path, case='thermal_mass', adaptive=False):
    with netCDF4.Dataset(path/(case+'_results.e')) as d:
        times=np.asarray(d['time_whole'][:])
        nodal=list(netCDF4.chartostring(d['name_nod_var'][:]))
        element=list(netCDF4.chartostring(d['name_elem_var'][:]))
        g=dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][:].T))
        x,y=np.asarray(d['coordx'][:]),np.asarray(d['coordy'][:])
        cross=x[:4]*np.roll(y[:4],-1)-np.roll(x[:4],-1)*y[:4]
        area=.5*np.sum(cross)
        mass=1000*.1*area
        temperature=300.
        axial_strain=0.
        records=[]
        for frame,time in enumerate(times):
            height=.1+.01*time
            if frame:
                divisions=2 if adaptive else 1
                sub_times=np.linspace(times[frame-1],time,divisions+1)
                for t0,t1 in zip(sub_times[:-1],sub_times[1:]):
                    dt=t1-t0
                    h0,h1=.1+.01*t0,.1+.01*t1
                    axial_strain+=2*(h1-h0)/(h1+h0)
                    # The independently identified CPEG8T source uses initial
                    # thickness; the in-plane area is unchanged in this case.
                    # The stable positive quadratic root gives the increment directly.
                    cp_old=100+.1*(temperature-300)
                    rhs=1e6*dt/1000
                    increment=2*rhs/(cp_old+np.sqrt(cp_old**2+.4*rhs))
                    temperature+=increment
            actual=np.asarray(d[f'vals_nod_var{nodal.index("temperature")+1}'][frame])
            np.testing.assert_array_equal(np.isfinite(actual),[True]*4+[False]*4)
            np.testing.assert_allclose(actual[:4],temperature,rtol=1e-11,atol=1e-10)
            for name in ['displacement_x','displacement_y']:
                np.testing.assert_allclose(d[f'vals_nod_var{nodal.index(name)+1}'][frame],0,rtol=0,atol=1e-13)
            np.testing.assert_allclose(g['section_solid_u3'][frame],.01*time,rtol=0,atol=1e-13)
            for component,modulus in [('xx',4e5),('yy',4e5),('zz',1.2e6),('xy',0)]:
                for q in range(9):
                    var=element.index(f'stress_{component}_q{q}')+1
                    np.testing.assert_allclose(d[f'vals_elem_var{var}eb1'][frame,0],modulus*axial_strain,
                                               rtol=1e-10,atol=1e-7)
            np.testing.assert_allclose(g['section_solid_axial_force'][frame],1.2e6*axial_strain*area,
                                       rtol=1e-10,atol=1e-9)
            if frame:
                rx=np.sum(d[f'vals_nod_var{nodal.index("reaction_x")+1}'][frame])
                ry=np.sum(d[f'vals_nod_var{nodal.index("reaction_y")+1}'][frame])
                right_length=np.hypot(x[2]-x[1],y[2]-y[1])
                # Top pressure is downward; bottom y traction is upward.
                expected_x=-mass*2-50*right_length*height
                expected_y=mass*9.81+100*(x[2]-x[3])*height-30*(x[1]-x[0])*.1
                np.testing.assert_allclose([rx,ry],[expected_x,expected_y],rtol=1e-10,atol=1e-9)
            records.append(dict(time=float(time),temperature=float(temperature),thickness=float(height),
                                axial_strain=float(axial_strain)))
        with (path/(case+'_history.csv')).open() as stream:
            history=list(csv.DictReader(stream))
        if len(history)!=len(times):
            raise AssertionError('Engineering history omits accepted steps')
        for frame,row in enumerate(history[1:],1):
            t=times[frame]
            expected=1e6*area*.1
            for key in ['generated_heat_rate','stored_heat_rate']:
                np.testing.assert_allclose(float(row[key]),expected,rtol=1e-10,atol=1e-9)
    (path/'thermal_mass_comparison.json').write_text(json.dumps(dict(status='passed',frames=records),indent=2)+'\n')
    print('Uniform thermal history, fixed initial mass, temperature-dependent heat capacity, all stresses, '
          'section force, gravity and current/reference surface load balance passed')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--adaptive',action='store_true')
    args=p.parse_args()
    source=Path(__file__).resolve().parent
    args.work.mkdir(parents=True,exist_ok=True)
    case='thermal_mass_adaptive' if args.adaptive else 'thermal_mass'
    for name in [case+'.fsi','distorted.e']:
        shutil.copyfile(source/name,args.work/name)
        assert (source/name).read_bytes()==(args.work/name).read_bytes()
    r=subprocess.run([str(args.executable.resolve()),'-i',case+'.fsi'],cwd=args.work,capture_output=True,text=True)
    (args.work/'production.log').write_text(r.stdout+r.stderr)
    if r.returncode or 'completed=true' not in r.stdout:
        raise RuntimeError(r.stdout+r.stderr)
    if args.adaptive and not any(line.startswith('progress.cutbacks=') and int(line.split('=')[1])>0
                                 for line in r.stdout.splitlines()):
        raise AssertionError('Expected initial rejected time steps were not exercised')
    check(args.work,case,args.adaptive)
