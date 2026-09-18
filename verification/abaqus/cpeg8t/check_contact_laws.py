"""Uniform contact law: independent series thermal resistance and native fields."""
import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np

from check_contact_supplement import run
from compare import compare, verify_references


def check(case, source, work):
    frames = json.loads((source/(case+'_fields.json')).read_text())['steps']['HISTORY']
    report, failures, failed_labels, history, opened = {}, [], [], [], set()
    def metric(a, r, z, absolute, label):
        try:
            compare(a, r, z, absolute, label, report)
        except AssertionError as error:
            failures.append(str(error))
            failed_labels.append(label)
    with netCDF4.Dataset(work/(case+'_results.e')) as d:
        names = list(netCDF4.chartostring(d['name_nod_var'][:]))
        gn = list(netCDF4.chartostring(d['name_glo_var'][:]))
        en = list(netCDF4.chartostring(d['name_elem_var'][:]))
        for at, time in enumerate(d['time_whole'][:]):
            if time == 0:
                continue
            pressure_case = case in ('thermal_pressure', 'thermal_pressure_small_gap')
            if case == 'thermal_pressure_small_gap':
                motion = np.interp(time, [0,.25,.5,.75,1], [0,0,-.000011,0,0])
                gap = 1e-6+motion
            else:
                motion = np.interp(time, [0, .5, 1], [0, -.0003 if pressure_case else .0003, 0])
                gap = .0001+motion
            pressure = max(-1e9*gap, 0)
            h = 1000+100*time+(.002*pressure if pressure_case else -1e6*gap)
            heat = 100*time/(1+1/(h*.002))
            if gap > 0:
                opened.add(str(float(time)))
            history.append(dict(time=float(time),gap_m=float(gap),pressure_Pa=float(pressure),
                                conductance_W_per_m2_K=float(h),analytic_heat_W=float(heat)))
            corners = [1, 2, 3, 4, 9, 10, 11, 12]
            temperature = [300, 300, 300+.5*heat, 300+.5*heat,
                           300+100*time-.5*heat, 300+100*time-.5*heat, 300+100*time, 300+100*time]
            reaction = [-heat/2, -heat/2, 0, 0, 0, 0, heat/2, heat/2]
            frame = next(f for f in frames if abs(f['time']-time) < 1e-12)['fields']
            for name, key, reference, zero, absolute in [
                ('temperature', 'NT11', temperature, [False]*8, 1e-9),
                ('heat_reaction', 'RFL11', reaction, [False, False, True, True, True, True, False, False], 1e-8),
            ]:
                a = np.asarray(d[f'vals_nod_var{names.index(name)+1}'][at])[np.asarray(corners)-1]
                native = {r['nodeLabel']: r['data'] for r in frame[key]['values']}
                metric(a, reference, zero, absolute, f'{time}/analytic/{name}')
                metric(a, [native[n] for n in corners], zero, absolute, f'{time}/native/{name}')
            gl = dict(zip(gn, d['vals_glo_var'][at]))
            total_heat = sum(v for k, v in gl.items() if k.startswith('contact_0_q') and k.endswith('_heat_rate'))
            total_force = sum(v for k, v in gl.items() if k.startswith('contact_0_q') and k.endswith('_force'))
            metric([total_heat], [heat], False, 1e-8, f'{time}/contact_heat')
            metric([total_force], [pressure*.002], pressure == 0, 1e-8, f'{time}/contact_force')
            for name, key, component in [('displacement_x', 'U', 0), ('displacement_y', 'U', 1),
                                         ('reaction_x', 'RF', 0), ('reaction_y', 'RF', 1)]:
                values = {r['nodeLabel']: r['data'][component] for r in frame[key]['values']}
                a = np.asarray(d[f'vals_nod_var{names.index(name)+1}'][at])
                zeros = [name.endswith('_x') or (name == 'displacement_y' and (n <= 8 or motion == 0))
                         or (name == 'reaction_y' and (pressure == 0 or n not in (3, 4, 7, 9, 10, 13)))
                         for n in range(1, 17)]
                metric(a, [values[n] for n in range(1, 17)], zeros,
                       1e-12 if name.startswith('displacement') else 1e-8, f'{time}/native/{name}')
            for component, label in enumerate(['xx', 'yy', 'zz', 'xy']):
                a, r = [], []
                for row in frame['S']['values']:
                    v = en.index(f'stress_{label}_q{row["integrationPoint"]-1}')+1
                    a.append(d[f'vals_elem_var{v}eb{row["elementLabel"]}'][at, 0])
                    r.append(row['data'][component])
                metric(a, r, True, 1e-7, f'{time}/native/stress_{label}')
    if case == 'thermal_pressure_small_gap':
        np.testing.assert_allclose([row['time'] for row in history],np.arange(1,9)*.125,rtol=0,atol=1e-14)
        assert len(opened)==5 and len(history)-len(opened)==3
    result = dict(status='failed' if failures else 'passed', failures=failures, metrics=report,
                  history=history,
                  analytical_status='failed' if any('/native/' not in k for k in failed_labels) else 'passed',
                  native_closed_status='failed' if any('/native/' in k and k.split('/')[0] not in opened
                                                       for k in failed_labels) else 'passed',
                  native_open_status='failed' if any('/native/' in k and k.split('/')[0] in opened
                                                     for k in failed_labels) else 'passed')
    (work/(case+'_comparison.json')).write_text(json.dumps(result, indent=2)+'\n')
    return not failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['thermal_pressure', 'thermal_clearance', 'thermal_pressure_small_gap'], required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    verify_references(source, 'supplement_reference.sha256')
    if args.case == 'thermal_pressure_small_gap':
        verify_references(source, 'small_gap_reference.sha256')
    if not run(args.case, source, args.work, args.executable) or not check(args.case, source, args.work):
        raise SystemExit('Contact law validation failed: '+str(args.work))
    print(args.case+': complete history passes analytical and native comparisons')
