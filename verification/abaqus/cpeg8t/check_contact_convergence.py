"""Manual surface refinement and penalty study, using complete static cards."""
import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np

from check_contact_supplement import run


def summarize(path):
    with netCDF4.Dataset(path) as d:
        g = dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][-1]))
        forces = [float(v) for k,v in g.items() if k.startswith('contact_0_q') and k.endswith('_force')]
        pressures = [float(v) for k,v in g.items() if k.startswith('contact_0_q') and k.endswith('_pressure')]
        gaps = [float(g[k.removesuffix('_pressure')+'_gap']) for k,v in g.items()
                if k.startswith('contact_0_q') and k.endswith('_pressure') and v>0]
        nodal = list(netCDF4.chartostring(d['name_nod_var'][:]))
        force_y = np.asarray(d[f'vals_nod_var{nodal.index("reaction_y")+1}'][-1])
        if abs(float(sum(force_y))) > 1e-8:
            raise AssertionError('Global vertical force balance failed')
        lower_nodes = np.unique(np.asarray(d['connect1'][:],dtype=int))-1
        return dict(force=abs(float(sum(force_y[lower_nodes]))),normal_scalar_force=sum(forces),
                    peak_pressure=max(pressures),maximum_penetration=-min(gaps),
                    active_constraints=sum(p>0 for p in pressures),
                    force_balance=float(sum(force_y)),pressure=pressures)


def study(source, work, executable=None):
    cases = [f'contact_refine_{n}' for n in (2,4,8,16)]+['contact_penalty_low','contact_penalty_high']
    results = {}
    for case in cases:
        if executable and not run(case,source,work/case,executable):
            raise AssertionError(case+' did not complete')
        directory = work/case if executable else work
        results[case] = summarize(directory/(case+'_results.e'))
    force = [results[f'contact_refine_{n}']['force'] for n in (2,4,8,16)]
    changes = np.diff(force)
    fine_change = abs(changes[-1])/abs(force[-1])
    results['surface_refinement'] = dict(elements=[2,4,8,16],successive_force_changes=changes.tolist(),
        last_relative_force_change=float(fine_change),
        monotonically_decreasing_changes=bool(np.all(np.diff(abs(changes))<0)),
        scope='Only the contact direction is refined; no continuum error bound or uniform convergence order claimed.')
    low,normal,high = [results[c] for c in ('contact_penalty_low','contact_refine_8','contact_penalty_high')]
    if not (low['maximum_penetration']>normal['maximum_penetration']>high['maximum_penetration']):
        raise AssertionError('Increasing penalty must reduce penetration in this prescribed indentation')
    results['penalty_sensitivity'] = dict(penalties=[5e8,1e9,2e9],
        forces=[r['force'] for r in (low,normal,high)],
        penetrations=[r['maximum_penetration'] for r in (low,normal,high)])
    (work/'contact_convergence.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True)
    args = p.parse_args()
    study(Path(__file__).resolve().parent,args.work,args.executable)
