"""Check exported production histories against closed-form uniaxial-strain updates.

For n=1 the Norton equation and the active hardening consistency equation are
linear, so no reference nonlinear solver or production library is used here.
The native comparison is reported independently; it must never inherit a pass
from the analytical check when Abaqus uses a different time integration rule.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from compare import compare, verify_references


def analytical(times, plastic, creep, finite=False, initial_plastic=False, initial_yield=1e3):
    shear, bulk, hardening = 4e5, 2e6/3, 1e4
    signed_stress, ep, ec, previous = 0., 0., 0., 0.
    plastic_tensor, creep_tensor = np.zeros(3), np.zeros(3)
    rows = []
    previous_extension, accumulated_strain = 0., 0.
    for frame, time in enumerate(times):
        extension = np.interp(time, [0, 1, 2, 3],
                              [0, .01, .02, .03] if initial_plastic else [0, .01, .01, .002])
        accumulated_strain += 2*(extension-previous_extension)/(2+extension+previous_extension)
        strain = accumulated_strain if finite else extension
        previous_extension = extension
        dt = time-times[frame-1] if frame else 0.
        trial = signed_stress + 2*shear*(strain-previous)
        direction = np.sign(trial)
        c = dt*1e-7 if creep else 0.
        q = abs(trial)/(1+3*shear*c)
        dp = 0.
        if plastic and q > initial_yield+hardening*ep:
            dp = (abs(trial)-(1+3*shear*c)*(initial_yield+hardening*ep)) / (
                3*shear+(1+3*shear*c)*hardening)
            q = initial_yield+hardening*(ep+dp)
        dc = c*q
        ep += dp
        ec += dc
        plastic_tensor += direction*dp*np.array([-.5,-.5,1.])
        creep_tensor += direction*dc*np.array([-.5,-.5,1.])
        signed_stress = direction*q
        stress = bulk*strain+signed_stress*np.array([-1/3,-1/3,2/3])
        elastic = np.array([0.,0.,strain])-plastic_tensor-creep_tensor
        rows.append(dict(stress=stress, elastic=elastic, plastic=plastic_tensor.copy(),
                         creep=creep_tensor.copy(), equivalent_plastic_strain=ep,
                         equivalent_creep_strain=ec))
        previous = strain
    return rows


def check(directory, source, finite=False, initial_plastic=False):
    report, native_report, failures = {}, {}, []
    case='inelastic_initial_plastic' if initial_plastic else ('inelastic_finite' if finite else 'inelastic_history')
    native = json.loads((source/(case+'_fields.json')).read_text())
    frames = native['steps']['HISTORY']
    with netCDF4.Dataset(directory/(case+'_results.e')) as d:
        times = np.asarray(d['time_whole'][:])
        np.testing.assert_allclose(times, np.arange(13)*.25, rtol=0, atol=1e-14)
        np.testing.assert_allclose(times, [f['time'] for f in frames], rtol=0, atol=1e-14)
        names = list(netCDF4.chartostring(d['name_elem_var'][:]))
        nodal_names = list(netCDF4.chartostring(d['name_nod_var'][:]))
        globals_ = dict(zip(netCDF4.chartostring(d['name_glo_var'][:]), d['vals_glo_var'][:].T))
        for block, label, plastic, creep in [(1,'plastic',True,False), (2,'creep',False,True),
                                              (3,'coupled',True,True)]:
            plastic = plastic or initial_plastic
            expected = analytical(times, plastic, creep, finite, initial_plastic,
                                  500. if initial_plastic and block == 2 else 1000.)
            prefix = f'section_{label}_'
            for name, reference in [('u3',np.interp(times,[0,1,2,3],
                                                   [0,.001,.002,.003] if initial_plastic else [0,.001,.001,.0002])),
                                    ('rotation_x',np.zeros(13)),('rotation_y',np.zeros(13)),
                                    ('axial_force',np.array([r['stress'][2] for r in expected])*.0002),
                                    ('moment_x',np.zeros(13)),('moment_y',np.zeros(13))]:
                np.testing.assert_allclose(globals_[prefix+name],reference,rtol=1e-10,atol=1e-10,
                                           err_msg=prefix+name)
            for quantity, native_key in [('stress','S'),('elastic','EE'),('plastic','PE'),('creep','CE'),
                                        ('equivalent_plastic_strain','PEEQ'),
                                        ('equivalent_creep_strain','CEEQ')]:
                scalar = quantity.startswith('equivalent')
                components = [''] if scalar else (['xx','yy','zz','xy'] if quantity=='stress'
                                                   else ['xx','yy','zz','xy','yz','xz'])
                for component, suffix in enumerate(components):
                    base = quantity+('_'+suffix if suffix else '')
                    ref = np.array([r[quantity] if scalar else (r[quantity][component] if component<3 else 0.)
                                    for r in expected])
                    actual = np.array([d[f'vals_elem_var{names.index(base+"_q"+str(q))+1}eb{block}'][:,0]
                                       for q in range(9)]).T
                    np.testing.assert_allclose(actual,np.repeat(ref[:,None],9,axis=1),rtol=1e-10,atol=1e-10,
                                               err_msg=label+' '+base)
                    report[label+' '+base] = dict(count=int(actual.size),
                                                 maximum_absolute=float(np.max(abs(actual-ref[:,None]))))
                    native_values = []
                    for frame in frames:
                        rows = {v['integrationPoint']:v['data'] for v in frame['fields'][native_key]['values']
                                if v['elementLabel']==block}
                        # Absent mechanisms have no native field entries. The analytical
                        # check above still requires all production entries to be zero.
                        if not rows:
                            native_values.append([0.]*9)
                        else:
                            native_values.append([rows[q+1] if scalar else
                                                  (rows[q+1][component] * (.5 if component>=3 else 1.)
                                                   if component<4 else 0.) for q in range(9)])
                    # Mechanism activation and exact zero components come from the
                    # analytical solution, never from a reference magnitude cutoff.
                    zeros = np.repeat((ref==0)[:,None],9,axis=1)
                    try:
                        compare(actual,np.array(native_values),zeros,1e-7 if quantity=='stress' else 1e-10,
                                label+' '+base,native_report)
                    except AssertionError as error:
                        failures.append(str(error))
            if plastic and not any(r['equivalent_plastic_strain']>0 for r in expected):
                raise AssertionError('Plasticity never activated')
            if creep and not any(r['equivalent_creep_strain']>0 for r in expected):
                raise AssertionError('Creep never activated')
            if plastic and creep:
                increments = np.diff([[r['equivalent_plastic_strain'],r['equivalent_creep_strain']]
                                      for r in expected],axis=0)
                if not np.any(np.all(increments>0,axis=1)):
                    raise AssertionError('No simultaneous mechanism increments')
            if initial_plastic:
                native_plastic = np.array([[v['data'] for v in frame['fields']['PEEQ']['values']
                                            if v['elementLabel']==block] for frame in frames])
                if native_plastic.shape != (13,9) or not np.all(np.diff(native_plastic,axis=0)>0):
                    raise AssertionError('Native plasticity must activate at every point from the first increment')
                if creep:
                    native_creep = np.array([[v['data'] for v in frame['fields']['CEEQ']['values']
                                              if v['elementLabel']==block] for frame in frames])
                    if native_creep.shape != (13,9) or not np.all(np.diff(native_creep,axis=0)>0):
                        raise AssertionError('Native creep must activate alongside plasticity at every increment')
            history_region=next(name for name in native['history']['HISTORY']
                                if name.startswith('Node ') and name.endswith('.'+str(100+block)))
            native_history=native['history']['HISTORY'][history_region]
            for name,key in [('u3','U3'),('rotation_x','UR1'),('rotation_y','UR2'),
                             ('axial_force','RF3'),('moment_x','RM1'),('moment_y','RM2')]:
                reference=np.array([row[1] for row in native_history[key]])
                zeros=np.ones(len(times),dtype=bool) if key in ('UR1','UR2','RM1','RM2') else times==0
                try:
                    compare(globals_[prefix+name],reference,zeros,1e-9,label+' '+name,native_report)
                except AssertionError as error:
                    failures.append(str(error))
        for name,key,component,zero_nodes in [('reaction_x','RF',0,(5,7)),('reaction_y','RF',1,(6,8))]:
            actual=np.asarray(d[f'vals_nod_var{nodal_names.index(name)+1}'][:])
            reference=[]
            for frame in frames:
                values={v['nodeLabel']:v['data'][component] for v in frame['fields'][key]['values']}
                reference.append([values[n] for n in range(1,25)])
            zeros=np.array([[time==0 or (n-1)%8+1 in zero_nodes for n in range(1,25)] for time in times])
            try:
                compare(actual,reference,zeros,1e-8,name,native_report)
            except AssertionError as error:
                failures.append(str(error))
        for name, value in [('temperature',300),('displacement_x',0),('displacement_y',0)]:
            actual = np.asarray(d[f'vals_nod_var{nodal_names.index(name)+1}'][:])
            valid = np.isfinite(actual)
            expected_count = 12 if name=='temperature' else 24
            if not np.all(np.sum(valid,axis=1)==expected_count):
                raise AssertionError('Incorrect exported field coverage')
            np.testing.assert_allclose(actual[valid],value,rtol=0,atol=1e-12)
    result = dict(analytical_status='passed',analytical=report,
                  native_status='failed' if failures else 'passed', native=native_report,
                  native_failures=failures)
    (directory/'inelastic_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('All 13 times, 27 material points, mechanism tensors and section forces match closed-form backward Euler')
    print(f'Native Abaqus comparison: {len(failures)} failed field metrics')
    return failures


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--require-native',action='store_true')
    parser.add_argument('--finite',action='store_true')
    parser.add_argument('--initial-plastic',action='store_true')
    args=parser.parse_args()
    if args.finite and args.initial_plastic:
        parser.error('--initial-plastic is a small-strain case')
    source=Path(__file__).resolve().parent
    verify_references(source, 'extended_reference.sha256')
    args.work.mkdir(parents=True,exist_ok=True)
    case='inelastic_initial_plastic' if args.initial_plastic else ('inelastic_finite' if args.finite else 'inelastic_history')
    for name in [case+'.fsi','inelastic.e']:
        shutil.copyfile(source/name,args.work/name)
        assert (source/name).read_bytes()==(args.work/name).read_bytes()
    run=subprocess.run([str(args.executable.resolve()),'-i',case+'.fsi'],
                       cwd=args.work,capture_output=True,text=True)
    (args.work/'production.log').write_text(run.stdout+run.stderr)
    if run.returncode or 'completed=true' not in run.stdout:
        raise RuntimeError(run.stdout+run.stderr)
    failed=check(args.work,source,args.finite,args.initial_plastic)
    if args.require_native and failed:
        raise AssertionError('\n'.join(failed))
