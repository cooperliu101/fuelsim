"""Production verification of curved nonmatching and transient thermal interfaces."""
import argparse
import json
from pathlib import Path
import re
import shutil
import numpy as np
from compare import compare
from study import balance, close, history, output, run


def scalar(log, key):
    match = re.search(r'^'+re.escape(key)+r'=([^\n]+)', log, re.M)
    if match is None:
        raise AssertionError('Missing production diagnostic '+key)
    return float(match.group(1))


def fields_equal(first, second, mask=None):
    if first.keys() != second.keys():
        raise AssertionError('Output variable sets differ')
    maximum_temperature = 0.
    for key in first:
        # Static coordinates/connectivity are compared separately from all saved
        # physical fields and time frames. No interpolation of differing time grids.
        values = first[key]
        if mask is not None and not key.startswith(('coord', 'connect')):
            values = values[mask]
        close(values, second[key], 'Output mismatch: '+key)
        if key == 'temperature':
            maximum_temperature = float(np.max(np.abs(values-second[key])))
    return maximum_temperature


def curved(args):
    name = 'study_dc3d20_curved_nonmatching'
    folder, log = run(args, name)
    native = compare(folder, name, args.source, [])
    if not native['passed']:
        raise AssertionError('Curved nonmatching native fields: '+json.dumps(native))
    a = output(folder, name)
    reaction = a['heat_reaction'][-1]
    hot_nodes = a['connect1'][0][[0, 3, 7, 4, 11, 15, 19, 12]]-1
    hot_heat = float(np.sum(reaction[hot_nodes]))
    cold_heat = float(np.sum(reaction)-hot_heat)
    close([hot_heat+cold_heat], [0.], 'Steady global heat conservation', atol=1e-10, rtol=0)
    close([hot_heat], [scalar(log, 'contact.gap.total_heat_rate')], 'Interface heat versus hot reaction', atol=1e-10)
    if scalar(log, 'contact.gap.unprojected_contact_nodes') != 0 or scalar(log, 'contact.gap.projected_contact_nodes') != 8:
        raise AssertionError('The curved secondary surface is not fully projected')
    if len(a['connect1']) != 1 or len(a['connect2']) != 2:
        raise AssertionError('Expected one secondary and two primary volume elements')
    return {'passed': True, 'native': native, 'hot_boundary_heat_W': hot_heat,
            'cold_boundary_heat_W': cold_heat, 'heat_imbalance_W': hot_heat+cold_heat,
            'projected_secondary_neighborhoods': 8}


def power_integral(t):
    return np.where(t <= 1, t+.5*t*t, 1.5+2*(t-1)-.75*(t-1)**2)


def check_energy(a, h):
    t, temperature = a['time_whole'], a['temperature']
    close(t, h['time'], 'Output/history times differ', atol=1e-10, rtol=0)
    if not np.all(np.diff(t) > 0) or t[-1] != 2:
        raise AssertionError('Accepted time sequence is incomplete or not monotone')
    source = 3*np.diff(power_integral(t))/np.diff(t)
    close(h['generated_heat_rate'][1:], source, 'Exact integrated source after retry', atol=1e-9)
    energy_hot = .1*np.mean(temperature[:, a['connect1'][0]-1]-300, axis=1)
    energy_cold = .4*np.mean(temperature[:, a['connect2'][0]-1]-300, axis=1)
    storage = np.diff(energy_hot+energy_cold)/np.diff(t)
    close(storage, h['stored_heat_rate'][1:], 'Independent initial-mass storage', atol=1e-8)
    defect = balance(h)
    secondary = a['connect1'][0][[1, 2, 6, 5]]-1
    primary = a['connect2'][0][[0, 3, 7, 4]]-1
    jump = np.mean(temperature[:, secondary], axis=1)-np.mean(temperature[:, primary], axis=1)
    if np.max(jump) < 5 or np.min(np.ptp(temperature, axis=1)) < 20:
        raise AssertionError('The interface/nonuniform-temperature branch was not exercised')
    return {'maximum_interface_temperature_jump_K': float(np.max(jump)),
            'minimum_global_temperature_range_K': float(np.min(np.ptp(temperature, axis=1))),
            'maximum_storage_error_W': float(np.max(np.abs(storage-h['stored_heat_rate'][1:]))),
            'maximum_global_heat_imbalance_W': defect,
            'maximum_interface_heat_imbalance_W': float(np.max(np.abs(h['interface_heat_imbalance']))),
            'integrated_generated_energy_J': float(np.sum(h['generated_heat_rate'][1:]*np.diff(t)))}


def transient(args):
    fixed_name = 'study_interface_variable_source'
    fixed, log = run(args, fixed_name)
    native = compare(fixed, fixed_name, args.source, ['y', 'z'])
    if not native['passed']:
        raise AssertionError('Variable-source native comparison: '+json.dumps(native))
    a, h = output(fixed, fixed_name), history(fixed, fixed_name)
    fixed_energy = check_energy(a, h)
    t, temperatures = a['time_whole'], a['temperature']
    hot = a['connect1'][0]-1
    cold = a['connect2'][0]-1
    interface_heat = .1*(np.mean(temperatures[:, hot[[1, 2, 6, 5]]], axis=1)-
                         np.mean(temperatures[:, cold[[0, 3, 7, 4]]], axis=1))
    convection = .01*(np.mean(temperatures[:, cold[[1, 2, 6, 5]]], axis=1)-300)
    source = np.diff(power_integral(t))/np.diff(t)
    hot_storage = .1*np.diff(np.mean(temperatures[:, hot], axis=1))/np.diff(t)
    cold_storage = .4*np.diff(np.mean(temperatures[:, cold], axis=1))/np.diff(t)
    hot_defect = hot_storage-source-h['dirichlet_heat_input_rate'][1:]+interface_heat[1:]
    cold_defect = cold_storage-2*source+convection[1:]-interface_heat[1:]
    close(hot_defect, np.zeros_like(hot_defect), 'Hot-region energy balance', atol=1e-8, rtol=0)
    close(cold_defect, np.zeros_like(cold_defect), 'Cold-region energy balance', atol=1e-8, rtol=0)

    full_name = 'study_interface_adaptive_full'
    first_name = 'study_interface_adaptive_first'
    resume_name = 'study_interface_adaptive_resume'
    full, full_log = run(args, full_name)
    split = args.work/'split1'
    run(args, first_name, split)
    run(args, resume_name, split)
    full_state, resumed_state = output(full, full_name), output(split, resume_name)
    mask = full_state['time_whole'] >= 1-1e-12
    restart_difference = fields_equal(full_state, resumed_state, mask)
    full_history = history(full, full_name)
    first_history, resumed_history = history(split, first_name), history(split, resume_name)
    for key in full_history:
        close(full_history[key], np.concatenate([first_history[key], resumed_history[key][1:]]), 'Restart history '+key)
    rejected = int(scalar(full_log, 'time_error_rejections'))
    if rejected == 0 or scalar(full_log, 'maximum_accepted_time_error_estimate') > 1:
        raise AssertionError('Adaptive rejection/acceptance contract was not exercised')
    adaptive_energy = check_energy(full_state, full_history)

    mpi_folder, mpi_log = run(args, full_name, args.work/'full2', ranks=2)
    mpi_state, mpi_history = output(mpi_folder, full_name), history(mpi_folder, full_name)
    mpi_difference = fields_equal(full_state, mpi_state)
    for key in full_history:
        close(full_history[key], mpi_history[key], 'MPI history '+key)
    if scalar(mpi_log, 'time_error_rejections') != rejected:
        raise AssertionError('MPI changed the rejected-step count')
    check_energy(mpi_state, mpi_history)

    # A checkpoint written by one process is resumed by two processes.
    cross_rank = args.work/'resume2_from_rank1'
    cross_rank.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(split/(first_name+'.checkpoint'), cross_rank/(first_name+'.checkpoint'))
    run(args, resume_name, cross_rank, ranks=2)
    cross_state = output(cross_rank, resume_name)
    cross_difference = fields_equal(full_state, cross_state, mask)
    cross_history = history(cross_rank, resume_name)
    for key in resumed_history:
        close(resumed_history[key], cross_history[key], 'Cross-rank restart history '+key)
    return {'passed': True, 'native': native, 'fixed_step_energy': fixed_energy,
            'maximum_hot_region_defect_W': float(np.max(np.abs(hot_defect))),
            'maximum_cold_region_defect_W': float(np.max(np.abs(cold_defect))),
            'adaptive_energy': adaptive_energy, 'time_error_rejections': rejected,
            'accepted_steps': len(full_state['time_whole'])-1,
            'restart_maximum_temperature_difference_K': restart_difference,
            'mpi_maximum_temperature_difference_K': mpi_difference,
            'cross_rank_restart_maximum_temperature_difference_K': cross_difference,
            'mpi_ranks': [1, 2]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fuelsim', type=Path, required=True)
    p.add_argument('--source', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--mpiexec', default='mpiexec')
    p.add_argument('--study', choices=['curved', 'transient'], required=True)
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    result = curved(args) if args.study == 'curved' else transient(args)
    (args.work/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
