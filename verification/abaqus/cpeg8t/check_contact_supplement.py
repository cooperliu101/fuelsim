"""Run complete static cards and audit every available native frame and field."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from compare import compare


def run(case, source, work, executable):
    work.mkdir(parents=True, exist_ok=True)
    card = source/(case+'.fsi')
    # Reading the literal mesh filename does not alter the problem definition.
    mesh = next(line.split('=', 1)[1].strip() for line in card.read_text().splitlines()
                if line.strip().startswith('file ='))
    for name in [case+'.fsi', mesh]:
        shutil.copyfile(source/name, work/name)
        assert (source/name).read_bytes() == (work/name).read_bytes()
    result = subprocess.run([str(executable.resolve()), '-i', case+'.fsi'], cwd=work,
                            capture_output=True, text=True)
    (work/(case+'.log')).write_text(result.stdout+result.stderr)
    return result.returncode == 0 and 'completed=true' in result.stdout


def check(case, source, work):
    native = json.loads((source/(case+'_fields.json')).read_text())
    frames = next(iter(native['steps'].values()))
    report, failures = {}, []
    def metric(a, r, z, absolute, label):
        try:
            compare(a, r, z, absolute, label, report)
        except AssertionError as error:
            failures.append(str(error))
        report.setdefault(label, {}).update(actual=np.asarray(a).tolist(), reference=np.asarray(r).tolist())
    with netCDF4.Dataset(source/'deforming_contact.e') as mesh:
        x = np.asarray(mesh['coordx'][:])
        sets = dict(zip(netCDF4.chartostring(mesh['ns_names'][:]),
                        [set(map(int, mesh[f'node_ns{i+1}'][:])) for i in range(len(mesh['ns_names']))]))
        corners = set()
        element_map = {}
        for b in range(1, len(mesh['eb_names'])+1):
            for e, nodes in enumerate(mesh[f'connect{b}'][:]):
                element_map[len(element_map)+1] = (b, e)
                corners.update(map(int, nodes[:4]))
        count = len(mesh['coordx'])
    with netCDF4.Dataset(work/(case+'_results.e')) as d:
        nn = list(netCDF4.chartostring(d['name_nod_var'][:]))
        en = list(netCDF4.chartostring(d['name_elem_var'][:]))
        times = np.asarray(d['time_whole'][:])
        for frame in frames:
            time = frame['time']
            if time == 0:
                continue
            at = int(np.argmin(abs(times-time)))
            if abs(times[at]-time) > 1e-12:
                failures.append(f'Missing production time {time}')
                continue
            fields = frame['fields']
            fixed_primary = case == 'deforming_fixed_primary'
            fixed = sets['lower'] if fixed_primary else sets['bottom']
            constrained = fixed | sets['top']
            for name, key, component, nodes, zeros, absolute in [
                ('temperature', 'NT11', None, sorted(corners), lambda n: False, 1e-9),
                ('displacement_x', 'U', 0, range(1, count+1),
                 lambda n: n in fixed or (time <= .25 and (n in sets['top'] or
                           (fixed_primary and n in sets['upper'] and abs(x[n-1]+.006)<1e-14))), 1e-12),
                ('displacement_y', 'U', 1, range(1, count+1), lambda n: n in fixed, 1e-12),
                ('reaction_x', 'RF', 0, range(1, count+1),
                 lambda n: n not in constrained or (fixed_primary and
                 ((n in sets['lower'] and n not in {3,4,7,10,13}) or
                  (n in sets['top'] and abs(x[n-1]+.006)<1e-14))), 1e-8),
                ('reaction_y', 'RF', 1, range(1, count+1),
                 lambda n: n not in constrained or (fixed_primary and n in sets['lower'] and
                                                   n not in {3,4,7,10,13}), 1e-8),
                ('heat_reaction', 'RFL11', None, sorted(corners),
                 lambda n: n not in sets['bottom'] | sets['top'], 1e-8),
            ]:
                rows = {v['nodeLabel']: v['data'] for v in fields[key]['values']}
                r = [rows[n] if component is None else rows[n][component] for n in nodes]
                a = np.asarray(d[f'vals_nod_var{nn.index(name)+1}'][at])[np.asarray(list(nodes))-1]
                metric(a, r, [zeros(n) for n in nodes], absolute, f'{time}/{name}')
            for key, prefix, absolute in [('S', 'stress', 1e-7), ('EE', 'elastic', 1e-10)]:
                for c, component in enumerate(['xx', 'yy', 'zz', 'xy']):
                    a, r, zeros = [], [], []
                    for value in fields[key]['values']:
                        b, e = element_map[value['elementLabel']]
                        variable = en.index(f'{prefix}_{component}_q{value["integrationPoint"]-1}')+1
                        a.append(float(d[f'vals_elem_var{variable}eb{b}'][at, e]))
                        r.append(value['data'][c]*(.5 if key == 'EE' and c == 3 else 1))
                        zeros.append((fixed_primary and b == 1) or (key == 'EE' and c == 2))
                    metric(a, r, zeros, absolute, f'{time}/{prefix}_{component}')
    result = {'status': 'failed' if failures else 'passed', 'failures': failures, 'metrics': report}
    (work/(case+'_comparison.json')).write_text(json.dumps(result, indent=2)+'\n')
    return not failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    if not run(args.case, source, args.work, args.executable):
        raise SystemExit('Production did not complete; inspect '+str(args.work/(args.case+'.log')))
    if not check(args.case, source, args.work):
        raise SystemExit('Comparison failed; inspect '+str(args.work/(args.case+'_comparison.json')))
