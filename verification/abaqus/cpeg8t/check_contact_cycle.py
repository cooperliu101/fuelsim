"""Independent flat-interface geometry, penalty force and thermal balance over a cycle."""
import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from check_restart import equivalent
from compare import compare, verify_references


def check(directory, source):
    records=[]
    distribution_actual,distribution_exact,distribution_zero=[],[],[]
    with netCDF4.Dataset(directory/'contact_cycle_results.e') as d:
        times=np.asarray(d['time_whole'][:])
        np.testing.assert_allclose(times,np.arange(13)*.25,rtol=0,atol=1e-14)
        names=list(netCDF4.chartostring(d['name_nod_var'][:]))
        g=dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][:].T))
        owners=[]
        for frame,time in enumerate(times):
            ramp=time/3
            ux=.015*ramp
            uy=np.interp(time,[0,.5,1,1.5,2,2.5,3],[0,-.0002,-.0002,0,0,-.0002,-.0002])
            gap=.0001+uy
            pressure=1e9*max(0,-gap)
            total_area=0.
            owner=[]
            # Two corner-neighborhood thermal constraints share the entire edge.
            # Here the rigid straight edge remains fully inside the primary surface.
            total_area = .01 * .1
            owner = [p for p,(left,right) in enumerate([(-.02,0.),(0.,.02)])
                     if min(ux,right) > max(ux-.01,left)+1e-14]
            for q in range(2):
                area = total_area / 2
                for name, expected, atol in [
                        ('gap', gap, 1e-13), ('pressure', 0., 1e-7),
                        ('area', area, 1e-13), ('force', 0., 1e-8),
                        ('heat_rate', 1000*(100*ramp)*area, 1e-8)]:
                    np.testing.assert_allclose(g[f'contact_0_q{q+3}_{name}'][frame],expected,
                                               rtol=1e-11,atol=atol,err_msg=f't={time}, q={q}, {name}')
            # Mechanical entries precede the separate thermal quadrature entries.
            for q in range(3):
                np.testing.assert_allclose(g[f'contact_0_q{q}_gap'][frame],gap,rtol=0,atol=1e-13)
                np.testing.assert_allclose(g[f'contact_0_q{q}_pressure'][frame],pressure,rtol=1e-11,atol=1e-7)
                np.testing.assert_allclose(g[f'contact_0_q{q}_force'][frame],
                    pressure*g[f'contact_0_q{q}_area'][frame],rtol=1e-11,atol=1e-8)
            np.testing.assert_allclose(sum(g[f'contact_0_q{q}_area'][frame] for q in range(3)),
                                       total_area,rtol=1e-11,atol=1e-13)
            owners.append(owner)
            for name,expected in [('u3',.01*ramp),('rotation_x',.02*ramp),('rotation_y',-.03*ramp)]:
                np.testing.assert_allclose(g['section_upper_'+name][frame],expected,rtol=0,atol=1e-13)
                np.testing.assert_allclose(g['section_lower_'+name][frame],0,rtol=0,atol=1e-13)
            for name,upper in [('displacement_x',ux),('displacement_y',uy)]:
                a=np.asarray(d[f'vals_nod_var{names.index(name)+1}'][frame])
                np.testing.assert_allclose(a,np.r_[np.zeros(13),np.full(8,upper)],rtol=0,atol=1e-13)
            temp=np.asarray(d[f'vals_nod_var{names.index("temperature")+1}'][frame])
            thermal_nodes=np.array([1,2,3,4,9,10,14,15,16,17])-1
            np.testing.assert_array_equal(np.flatnonzero(np.isfinite(temp)),thermal_nodes)
            np.testing.assert_allclose(temp[thermal_nodes],np.r_[np.full(6,300),np.full(4,300+100*ramp)],
                                       rtol=0,atol=1e-10)
            if frame:
                heat=np.asarray(d[f'vals_nod_var{names.index("heat_reaction")+1}'][frame])
                heat_rate=1000*100*ramp*total_area
                storage=1000*100*.1*.01*.01*(100/3)
                np.testing.assert_allclose([np.nansum(heat[:13]),np.nansum(heat[13:])],
                                           [-heat_rate,heat_rate+storage],rtol=1e-10,atol=1e-8)
                force=np.asarray(d[f'vals_nod_var{names.index("reaction_y")+1}'][frame])
                np.testing.assert_allclose(np.sum(force[:13]),pressure*total_area,rtol=1e-10,atol=1e-8)
                # Exact integration of each primary quadratic shape over its
                # overlap with the translating secondary edge. Rational bounds
                # identify analytical zeros without a magnitude threshold.
                exact=[Fraction(0) for _ in range(13)]
                if frame in (2,3,4,10,11,12):
                    for center,nodes in [(-8,[3,2,6]),(8,[2,9,12])]:
                        left,right=max(center-8,frame-8),min(center+8,frame)
                        if right<=left:
                            continue
                        a,b=Fraction(left-center,8),Fraction(right-center,8)
                        integrals=[(b**3-a**3)/6-(b*b-a*a)/4,
                                   (b**3-a**3)/6+(b*b-a*a)/4,
                                   b-a-(b**3-a**3)/3]
                        for node,value in zip(nodes,integrals):
                            exact[node]+=100*value
                distribution_actual.append(force[:13])
                distribution_exact.append([float(v) for v in exact])
                distribution_zero.append([v==0 for v in exact])
            records.append(dict(time=float(time),gap=float(gap),pressure=float(pressure),area=float(total_area),
                                primary_elements=[p+1 for p in owner]))
        if owners[0]!=[0] or owners[-1]!=[1] or not any(len(set(o))==2 for o in owners):
            raise AssertionError('The cycle did not cross the primary element boundary')
    with (directory/'contact_cycle_history.csv').open() as stream:
        history=list(csv.DictReader(stream))
    for row in history[1:]:
        np.testing.assert_allclose(float(row['stored_heat_rate']),100/3,rtol=1e-10,atol=1e-9)
        np.testing.assert_allclose(float(row['interface_heat_imbalance']),0,rtol=0,atol=1e-9)
    distribution_report,pointwise_report,pointwise_failures={},{},[]
    try:
        compare(distribution_actual,distribution_exact,distribution_zero,1e-8,
                'primary_nodal_contact_force',pointwise_report)
    except AssertionError as error:
        pointwise_failures.append(str(error))
    native=json.loads((source/'contact_cycle_fields.json').read_text())['steps']['CYCLE']
    native_forces=[]
    for frame in native[1:]:
        values={v['nodeLabel']:v['data'][1] for v in frame['fields']['RF']['values']}
        native_forces.append([values[n] for n in range(1,14)])
    native_failures=[]
    try:
        compare(distribution_actual,native_forces,distribution_zero,1e-8,
                'native_primary_nodal_contact_force',distribution_report)
    except AssertionError as error:
        native_failures.append(str(error))
    result=dict(balance_status='passed',frames=records,
                nodal_distribution_status='failed' if native_failures else 'passed',
                nodal_distribution=distribution_report,failures=native_failures,
                pointwise_integral_diagnostic=dict(metrics=pointwise_report,differences=pointwise_failures),
                native_distribution_status='failed' if native_failures else 'passed',
                native_distribution_failures=native_failures)
    (directory/'cycle_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Contact closure, release, reclosure, changing thickness, nonmatching element crossing, '
          'all contact quadrature quantities and transient heat balance passed')
    print(f'Native averaged primary nodal force distribution: {len(native_failures)} failed metrics')
    return native_failures


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--mpiexec',type=Path)
    p.add_argument('--require-distribution',action='store_true')
    args=p.parse_args()
    source=Path(__file__).resolve().parent
    verify_references(source, 'extended_reference.sha256')
    for mode in (['serial','parallel'] if args.mpiexec else ['serial']):
        work=args.work/mode
        work.mkdir(parents=True,exist_ok=True)
        for name in ['contact_cycle.fsi','sliding.e']:
            shutil.copyfile(source/name,work/name)
            assert (source/name).read_bytes()==(work/name).read_bytes()
        command=[str(args.executable.resolve()),'-i','contact_cycle.fsi']
        if mode=='parallel':
            command=[str(args.mpiexec),'-n','4']+command
        r=subprocess.run(command,cwd=work,capture_output=True,text=True)
        (work/'production.log').write_text(r.stdout+r.stderr)
        if r.returncode or 'completed=true' not in r.stdout:
            raise RuntimeError(r.stdout+r.stderr)
        failures=check(work,source)
        if args.require_distribution and failures:
            raise AssertionError('\n'.join(failures))
    if args.mpiexec:
        equivalent(args.work/'parallel/contact_cycle_results.e',args.work/'serial/contact_cycle_results.e')
        print('Four-process fields and histories agree with serial throughout contact reassignment')
