"""Black-box thermal studies: unchanged cards, production outputs, independent checks."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import netCDF4
import numpy as np
from compare import compare, metric


def run(args, name, directory=None, ranks=1):
    directory = directory or args.work/name
    directory.mkdir(parents=True, exist_ok=True)
    # A repeated test must not accidentally inspect a prior successful result.
    for pattern in [name+'_results*.e', name+'_history*.csv', name+'.checkpoint']:
        for old in directory.glob(pattern):
            old.unlink()
    card = args.source/(name+'.fsi')
    mesh = re.search(r'\[Mesh\].*?file\s*=\s*(\S+)', card.read_text(), re.S).group(1)
    for filename in [card.name, mesh]:
        shutil.copyfile(args.source/filename, directory/filename)
        if (directory/filename).read_bytes() != (args.source/filename).read_bytes():
            raise AssertionError('Input changed while staging '+filename)
    command = [str(args.fuelsim.resolve()), '-i', str((directory/card.name).resolve())]
    if ranks > 1:
        command = [args.mpiexec, '-n', str(ranks)] + command
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    log = result.stdout+result.stderr
    (directory/(name+'.log')).write_text(log)
    if result.returncode:
        raise RuntimeError(log)
    return directory, log


def output(directory, name):
    files = list(directory.glob(name+'_results*.e'))
    if len(files) != 1:
        raise AssertionError('Expected exactly one output file: '+str(files))
    with netCDF4.Dataset(files[0]) as data:
        result = {key: np.asarray(value[:]) for key, value in data.variables.items()
                  if key.startswith(('vals_', 'coord', 'connect')) or key == 'time_whole'}
        names = netCDF4.chartostring(data['name_nod_var'][:]).tolist()
        for i, key in enumerate(names, 1):
            result[key] = np.asarray(data['vals_nod_var'+str(i)][:])
        names = netCDF4.chartostring(data['name_elem_var'][:]).tolist()
        for i, key in enumerate(names, 1):
            result[key] = np.concatenate([np.asarray(data[f'vals_elem_var{i}eb{b+1}'][:])
                                          for b in range(len(data.dimensions['num_el_blk']))], axis=1)
    return result


def history(directory, name):
    files = list(directory.glob(name+'_history*.csv'))
    if len(files) != 1:
        raise AssertionError('Expected exactly one history file: '+str(files))
    with files[0].open() as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.asarray([row[key] for row in rows], float) for key in rows[0]}


def close(a, b, message, atol=1e-9, rtol=1e-11):
    if np.shape(a) != np.shape(b) or not np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=True):
        raise AssertionError(message)


def balance(h):
    net = (h['generated_heat_rate']+h['dirichlet_heat_input_rate']+
           h['surface_heat_input_rate']-h['convection_heat_rate']-h['stored_heat_rate'])
    close(net, np.zeros_like(net), 'Independent rate balance', atol=1e-8, rtol=0)
    close(h['interface_heat_imbalance'], np.zeros_like(net), 'Interface conservation', atol=1e-10, rtol=0)
    return float(np.max(np.abs(net)))


def native(args):
    results = []
    for name in ['study_dc3d20_nonmatching', 'study_dc3d20_curved', 'study_dcax4_k',
                 'study_dc3d20_k', 'study_dcax4_cp', 'study_dc3d20_cp', 'study_dc3d8_cp_capacity']:
        folder, log = run(args, name)
        zero = [] if name.endswith(('nonmatching', 'curved')) else list('xyz' if name.endswith('capacity') else 'yz')
        result = compare(folder, name, args.source, zero)
        results.append(result)
    # Compare the prescribed variable-cp ramp against the constant-cp native probe.
    # The end-point rule predicts exactly twice every nodal reaction, while the
    # enthalpy-increment rule predicts 1.5 times, independent of geometric weights.
    def reactions(name):
        with (args.source/(name+'.nodes.csv')).open() as stream:
            return np.array([float(row['reaction']) for row in csv.DictReader(stream)])
    ratio = reactions('study_dc3d8_cp_capacity')/reactions('dc3d8_capacity')
    close(ratio, np.full_like(ratio, 2), 'Native end-point specific heat rule', atol=2e-6, rtol=0)
    return {'cases': results, 'native_capacity_reaction_ratios': ratio.tolist(),
            'passed': all(r['passed'] for r in results)}


def restart(args):
    full, log = run(args, 'study_adaptive_full')
    split = args.work/'split'
    run(args, 'study_adaptive_first', split)
    run(args, 'study_adaptive_resume', split)
    count = re.search(r'^time_error_rejections\s*[=:]\s*(\d+)', log, re.M)
    if count is None or int(count.group(1)) == 0:
        raise AssertionError('The requested error-based time-step rejection did not occur')
    a = output(full, 'study_adaptive_full')
    b = output(split, 'study_adaptive_resume')
    mask = a['time_whole'] >= 1-1e-12
    for key in a:
        if key.startswith('vals_') or key == 'time_whole':
            close(a[key][mask], b[key], 'Restart changed '+key)
    h = history(full, 'study_adaptive_full')
    first = history(split, 'study_adaptive_first')
    resumed = history(split, 'study_adaptive_resume')
    for key in h:
        close(h[key], np.concatenate([first[key], resumed[key][1:]]), 'Restart history changed '+key)
    error = balance(h)
    t = a['time_whole']
    # Integral of the independently prescribed piecewise-linear source.
    integral = np.where(t <= 1, t+.5*t*t, 1.5+2*(t-1)-.75*(t-1)**2)
    exact = 200+np.sqrt(10000+2000*integral)
    maximum = float(np.max(np.abs(a['temperature']-exact[:, None])))
    if maximum > .02:
        raise AssertionError('Adaptive solution differs from analytical heating by '+str(maximum))
    # Accepted solution consists of two half steps. Source diagnostics must average
    # both end-point loads, including the load that precedes the midpoint.
    midpoint = t[1:]-.5*np.diff(t)
    source = lambda time: np.interp(time, [0, 1, 2], [1, 2, .5])
    close(h['generated_heat_rate'][1:], source(midpoint)+source(t[1:]), 'Source was not restored/averaged')
    return {'passed': True, 'rejected_steps': int(count.group(1)),
            'accepted_steps': len(t)-1, 'maximum_temperature_error_K': maximum,
            'maximum_energy_rate_defect_W': error}


def energy(args):
    results = []
    for name in ['study_dc3d8_shared', 'study_dc3d8_interface']:
        folder, log = run(args, name)
        a, h = output(folder, name), history(folder, name)
        defect = balance(h)
        t = a['temperature']
        energy_value = .1*np.mean(t[:, a['connect1'][0]-1]-300, axis=1)+.4*np.mean(t[:, a['connect2'][0]-1]-300, axis=1)
        storage = np.diff(energy_value)/np.diff(a['time_whole'])
        close(storage, h['stored_heat_rate'][1:], 'Independent nodal energy increment', atol=1e-8)
        close(h['generated_heat_rate'][1:], np.full(len(storage), 3.), 'Two material source integration')
        convection = .01*(np.mean(t[:, np.isclose(a['coordx'], .02)], axis=1)-300)
        close(convection[1:], h['convection_heat_rate'][1:], 'Independent boundary heat loss')
        native_result = compare(folder, name, args.source, ['y', 'z'])
        if not native_result['passed']:
            raise AssertionError('Multi-material native comparison: '+json.dumps(native_result))
        results.append({'case': name, 'maximum_energy_rate_defect_W': defect,
                        'maximum_storage_error_W': float(np.max(np.abs(storage-h['stored_heat_rate'][1:]))),
                        'native': native_result})
    name = 'study_dc3d8_interface'
    parallel, log = run(args, name, args.work/'mpi2', ranks=2)
    b = output(parallel, name)
    for key in a:
        close(a[key], b[key], 'MPI changed '+key)
    hb = history(parallel, name)
    for key in h:
        close(h[key], hb[key], 'MPI history changed '+key)
    return {'passed': True, 'cases': results, 'mpi_ranks': [1, 2]}


def convergence(args):
    space = {}
    for kind in ['dcax4', 'dcax8']:
        errors, flux_errors, heat_errors = [], [], []
        for n in [4, 8, 16]:
            name = 'study_'+kind+'_n'+str(n)
            folder, log = run(args, name)
            a = output(folder, name)
            exact = 400-100*np.log(a['coordx']/.01)/np.log(3)
            errors.append(float(np.max(np.abs(a['temperature'][-1]-exact))))
            q_errors = []
            for q in range(4 if kind == 'dcax4' else 9):
                exact_q = 1000/(a[f'point_x_q{q}'][-1]*np.log(3))
                q_errors.extend(np.abs(a[f'heat_flux_x_q{q}'][-1]/exact_q-1))
                close(a[f'heat_flux_y_q{q}'][-1], np.zeros(n), 'Spurious axial heat flux', atol=1e-7, rtol=0)
            flux_errors.append(float(max(q_errors)))
            heat = np.sum(a['heat_reaction'][-1, np.isclose(a['coordx'], .01)])
            heat_errors.append(float(abs(heat/(20*np.pi/np.log(3))-1)))
        rates = np.log2(np.array(errors[:-1])/errors[1:])
        if np.min(rates) < (1.8 if kind == 'dcax4' else 2.8):
            raise AssertionError('Spatial convergence order '+kind+': '+str(rates))
        flux_rates = np.log2(np.array(flux_errors[:-1])/flux_errors[1:])
        if np.min(flux_rates) < (.8 if kind == 'dcax4' else 1.7):
            raise AssertionError('Heat-flux convergence order '+kind+': '+str(flux_rates))
        space[kind] = {'maximum_temperature_errors_K': errors, 'nodal_observed_orders': rates.tolist(),
                       'maximum_relative_flux_errors': flux_errors, 'flux_observed_orders': flux_rates.tolist(),
                       'relative_total_heat_errors': heat_errors}
    temporal = []
    for level in ['coarse', 'medium', 'fine']:
        name = 'study_time_'+level
        folder, log = run(args, name)
        a, h = output(folder, name), history(folder, name)
        balance(h)
        exact = 200+np.sqrt(14000)
        error = float(np.max(np.abs(a['temperature'][-1]-exact)))
        temporal.append(error)
        if np.max(np.ptp(a['temperature'], axis=1)) > 1e-9:
            raise AssertionError('Uniform heating acquired a spatial gradient')
    rates = np.log2(np.array(temporal[:-1])/temporal[1:])
    if np.min(rates) < .9 or np.max(rates) > 1.1:
        raise AssertionError('Backward Euler temporal order: '+str(rates))
    if temporal[-1]/(np.sqrt(14000)-100) >= .005:
        raise AssertionError('Fine time-grid temperature rise misses 0.5 percent')
    return {'passed': True, 'spatial': space, 'temporal_errors_K': temporal,
            'temporal_orders': rates.tolist()}


def analytic(args):
    results = []
    for name in ['study_dcax4_n128', 'study_dcax8_n16', 'study_time_fine']:
        folder, log = run(args, name)
        a = output(folder, name)
        if name == 'study_time_fine':
            exact = 200+np.sqrt(10000+2000*a['time_whole'])
            exact = exact[:, None]+np.zeros_like(a['temperature'])
            balance(history(folder, name))
        else:
            exact = 400-100*np.log(a['coordx']/.01)/np.log(3)
            exact = np.broadcast_to(exact, a['temperature'][-1:].shape)
        actual = a['temperature'] if name == 'study_time_fine' else a['temperature'][-1:]
        rise = exact-300
        nonzero = rise != 0
        relative = float(np.max(np.abs((actual-exact)[nonzero]/rise[nonzero])))
        if relative >= .005:
            raise AssertionError('Analytical temperature-rise error '+name+': '+str(relative))
        close(actual[~nonzero], exact[~nonzero], 'Zero temperature rise', atol=1e-7, rtol=0)
        result = {'case': name, 'maximum_pointwise_temperature_rise_error': relative}
        if name == 'study_dcax4_n128':
            flux = np.concatenate([a[f'heat_flux_x_q{q}'][-1] for q in range(4)])
            radius = np.concatenate([a[f'point_x_q{q}'][-1] for q in range(4)])
            transverse = np.concatenate([a[f'heat_flux_{d}_q{q}'][-1] for d in 'yz' for q in range(4)])
            heat = np.sum(a['heat_reaction'][-1, np.isclose(a['coordx'], .01)])
            fields = {'radial_flux': metric(flux, 1000/(radius*np.log(3))),
                      'transverse_flux': metric(transverse, np.zeros_like(transverse), zero=True),
                      'total_hot_boundary_heat': metric([heat], [20*np.pi/np.log(3)])}
            if not all(field['passed'] for field in fields.values()):
                raise AssertionError('Selected 128-element cylinder misses the heat-flow gate: '+json.dumps(fields))
            result['heat_flow_fields'] = fields
        results.append(result)
    return {'passed': True, 'cases': results}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fuelsim', type=Path, required=True)
    p.add_argument('--source', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--mpiexec', default='mpiexec')
    p.add_argument('--study', choices=['native', 'restart', 'energy', 'convergence', 'analytic'], required=True)
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    result = globals()[args.study](args)
    (args.work/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise AssertionError('Original thermal comparison thresholds not met')


if __name__ == '__main__':
    main()
