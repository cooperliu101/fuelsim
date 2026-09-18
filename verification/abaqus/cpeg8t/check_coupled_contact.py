"""Audit coupled production history, contact transitions and native material fields."""
import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np

from check_contact_supplement import run
from compare import compare, verify_references


def check(source, work):
    case = 'coupled_contact'
    native = json.loads((source/(case+'_centered_fields.json')).read_text())
    frames = native['steps']['HISTORY']
    report, failures = {}, []
    def metric(a, r, z, absolute, label):
        try:
            compare(a, r, z, absolute, label, report)
        except AssertionError as error:
            failures.append(str(error))
    with netCDF4.Dataset(work/(case+'_results.e')) as d:
        nn = list(netCDF4.chartostring(d['name_nod_var'][:]))
        en = list(netCDF4.chartostring(d['name_elem_var'][:]))
        gn = list(netCDF4.chartostring(d['name_glo_var'][:]))
        times = np.asarray(d['time_whole'][:])
        if len(times) != 17 or len(frames) != 17:
            raise AssertionError('Both solvers must export all 17 prescribed times')
        active = []
        for at, time in enumerate(times):
            frame = next(f for f in frames if abs(f['time']-time)<1e-12)['fields']
            if time == 0:
                for name in ['temperature','displacement_x','displacement_y','reaction_x','reaction_y',
                             'heat_reaction']:
                    a = np.asarray(d[f'vals_nod_var{nn.index(name)+1}'][at])
                    a = a[np.isfinite(a)]
                    metric(a,np.full_like(a,300 if name == 'temperature' else 0),name != 'temperature',
                           1e-12 if name.startswith('displacement') else 1e-8,'0/initial/'+name)
                for b in (1,2):
                    for i,name in enumerate(en,1):
                        if name.startswith(('stress_','elastic_','plastic_','creep_','equivalent_')):
                            a = np.asarray(d[f'vals_elem_var{i}eb{b}'][0])
                            metric(a,np.zeros_like(a),True,1e-7 if name.startswith('stress_') else 1e-10,
                                   f'0/initial/block{b}/{name}')
                for key in ('NT11','U','RF','RFL11','S','EE','PE','CE','PEEQ','CEEQ'):
                    a = np.asarray([v['data'] for v in frame[key]['values'] if v.get('nodeLabel',0)<=16])
                    metric(a,np.full_like(a,300 if key == 'NT11' else 0),key != 'NT11',
                           1e-7 if key == 'S' else 1e-10,'0/native/'+key)
                continue
            gap_key = next(k for k in frame if k.startswith('COPEN'))
            # An everywhere-open interface carries no normal reaction. This
            # mask comes from signed geometry, not a force magnitude cutoff.
            open_interface = all(v['data'] > 0 for v in frame[gap_key]['values'])
            for name, key, c, nodes, zero, absolute in [
                ('temperature', 'NT11', None, [1,2,3,4,9,10,11,12], lambda n: False, 1e-9),
                ('displacement_x', 'U', 0, range(1,17), lambda n: True, 1e-12),
                ('displacement_y', 'U', 1, range(1,17), lambda n: n in (1,2,5) or
                 (time == .75 and n in (11,12,15)), 1e-12),
                ('reaction_x', 'RF', 0, range(1,17), lambda n: n in (5,7,13,15), 1e-8),
                ('reaction_y', 'RF', 1, range(1,17),
                 lambda n: open_interface or n not in (1,2,5,11,12,15), 1e-8),
                ('heat_reaction', 'RFL11', None, [1,2,3,4,9,10,11,12],
                 lambda n: n not in (1,2,11,12), 1e-8),
            ]:
                values = {r['nodeLabel']: r['data'] for r in frame[key]['values']}
                reference = [values[n] if c is None else values[n][c] for n in nodes]
                actual = np.asarray(d[f'vals_nod_var{nn.index(name)+1}'][at])[np.asarray(list(nodes))-1]
                metric(actual, reference, [zero(n) for n in nodes], absolute, f'{time}/{name}')
            for key, prefix in [('S','stress'), ('EE','elastic'), ('PE','plastic'), ('CE','creep'),
                                ('PEEQ','equivalent_plastic_strain'), ('CEEQ','equivalent_creep_strain')]:
                scalar = key in ('PEEQ','CEEQ')
                components = [''] if scalar else (['xx','yy','zz','xy'] if key == 'S' else
                                                   ['xx','yy','zz','xy','yz','xz'])
                for c, component in enumerate(components):
                    actual, reference = [], []
                    for r in frame[key]['values']:
                        suffix = '' if scalar else '_'+component
                        v = en.index(f'{prefix}{suffix}_q{r["integrationPoint"]-1}')+1
                        actual.append(float(d[f'vals_elem_var{v}eb{r["elementLabel"]}'][at,0]))
                        reference.append(r['data'] if scalar else (0.0 if c >= 4 else
                                         r['data'][c]*(.5 if c == 3 and key != 'S' else 1)))
                    metric(actual, reference, not scalar and c >= 3,
                           1e-7 if key == 'S' else 1e-10, f'{time}/{prefix}_{component}')
            gl = dict(zip(gn,d['vals_glo_var'][at]))
            force = sum(v for k,v in gl.items() if k.startswith('contact_0_q') and k.endswith('_force'))
            active.append(bool(force>1e-8))
            native_force_key = next(k for k in frame if k.startswith('CNORMF'))
            native_force = sum(v['data'][1] for v in frame[native_force_key]['values']
                               if v['nodeLabel'] in (9,10,13))
            metric([force],[native_force],open_interface,1e-8,f'{time}/contact_force')
            for section, node in [('lower',17),('upper',18)]:
                region = next(k for k in native['history']['HISTORY'] if k.endswith('.'+str(node)))
                history = native['history']['HISTORY'][region]
                for name,key in [('u3','U3'),('rotation_x','UR1'),('rotation_y','UR2'),
                                 ('axial_force','RF3'),('moment_x','RM1'),('moment_y','RM2')]:
                    h = {k:next(v for t,v in pairs if abs(t-time)<1e-12) for k,pairs in history.items()}
                    # Native control nodes are at the same initial centroids.
                    # Never subtract single-precision moment/force histories.
                    r = h[key]
                    metric([gl['section_'+section+'_'+name]],[r],key in ('UR1','UR2','RM2'),
                           1e-9,f'{time}/section_{section}_{name}')
        increments = []
        for b in (1,2):
            for q in range(9):
                ep = np.asarray(d[f'vals_elem_var{en.index(f"equivalent_plastic_strain_q{q}")+1}eb{b}'][:,0])
                ec = np.asarray(d[f'vals_elem_var{en.index(f"equivalent_creep_strain_q{q}")+1}eb{b}'][:,0])
                increments.append(bool(np.any((np.diff(ep)>0)&(np.diff(ec)>0))))
        if not all(increments):
            failures.append('Every material point must exhibit simultaneous plasticity and creep')
        native_simultaneous = 0
        for b in (1,2):
            for q in range(1,10):
                series = []
                for f in frames:
                    series.append([next(v['data'] for v in f['fields'][k]['values']
                                        if v['elementLabel']==b and v['integrationPoint']==q)
                                   for k in ('PEEQ','CEEQ')])
                native_simultaneous += bool(np.any(np.all(np.diff(series,axis=0)>0,axis=1)))
        if native_simultaneous != 18:
            failures.append('Native reference must activate both mechanisms at all 18 points')
        if not (any(active[:8]) and not active[11] and any(active[12:])):
            failures.append('Required closure, release and recontact were not observed')
    result = dict(status='failed' if failures else 'passed', failures=failures, metrics=report,
                  simultaneous_points=sum(increments), native_simultaneous_points=native_simultaneous,
                  contact_active=active)
    (work/'coupled_contact_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    return not failures


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True)
    args = p.parse_args()
    source = Path(__file__).resolve().parent
    verify_references(source,'supplement_reference.sha256')
    if not run('coupled_contact',source,args.work,args.executable) or not check(source,args.work):
        raise SystemExit('Coupled contact validation failed: '+str(args.work))
