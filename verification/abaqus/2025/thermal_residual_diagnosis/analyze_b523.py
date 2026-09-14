"""Reconstruct boundary thermal operators from saved nodal histories; no solve."""
import csv
import importlib.util
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ABAQUS = HERE.parents[1]
spec = importlib.util.spec_from_file_location('weights', ABAQUS / 'b523_diagnosis/analyze_thermal_weights.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
PREFIX = HERE.parent / 'production/b523_hex8_c3d8t_integrated_path'
reference, elements = w.geometry(Path(str(PREFIX) + '.inp'))
nodes, points, energies = w.load_native(PREFIX)
actual = w.load_actual()

def operators(temperature, old_temperature, positions):
    fixed, native, conduction = (np.zeros(24) for _ in range(3))
    for element, indices in elements.items():
        ref = [w.metrics(reference[indices], sign / np.sqrt(3)) for sign in w.SIGNS]
        cur = [w.metrics(positions[indices], sign / np.sqrt(3)) for sign in w.SIGNS]
        ratio = sum(m[2] for m in cur) / sum(m[2] for m in ref)
        t = temperature[indices]
        rate = (t - old_temperature[indices]) / .02
        for q, index in enumerate(indices):
            capacity = ref[q][2] * 100 * (1 + .001 * (t[q] - 300)) * rate[q]
            fixed[index] += capacity
            native[index] += capacity * ratio / (cur[q][2] / ref[q][2])
            k = (15 if element <= 2 else 10) * (1 + .001 * (t[q] - 300))
            g = cur[q][1]
            conduction[indices] += ref[q][2] * ratio * k * (g @ (t @ g))
    return fixed, native, conduction

rows = []
old_native = np.full(24, 300.)
for increment in range(1, 21):
    t = np.array([nodes[increment, n]['temperature_k'] for n in range(1, 25)])
    x = reference + np.array([[nodes[increment, n]['u%d_m' % i] for i in range(1, 4)] for n in range(1, 25)])
    fixed, native, cond = operators(t, old_native, x)
    tf = actual['temperature'][increment]
    xf = reference + np.column_stack([actual['displacement_' + s][increment] for s in 'xyz'])
    ff, unused, fc = operators(tf, actual['temperature'][increment - 1], xf)
    for n in w.BOUNDARY:
        i = n - 1
        ra = nodes[increment, n]['reaction_heat_flux_w']
        rf = actual['reaction_heat_flux'][increment, i]
        rows.append(dict(time=increment * .02, node=n, native_reaction_w=ra, fuelsim_reaction_w=rf,
            total_difference_w=rf-ra, native_reconstruction_error_w=native[i]+cond[i]-ra,
            fuelsim_reconstruction_error_w=ff[i]+fc[i]-rf,
            capacity_rule_difference_w=fixed[i]-native[i],
            state_history_difference_w=ff[i]+fc[i]-fixed[i]-cond[i]))
    old_native = t
with (HERE/'b523_decomposition.tsv').open('w') as f:
    writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    writer.writeheader();writer.writerows(rows)
for key in ['native_reconstruction_error_w','fuelsim_reconstruction_error_w']:
    print(key,max(abs(row[key]) for row in rows))
print('worst sample',max(rows,key=lambda row:abs(row['total_difference_w'])))
