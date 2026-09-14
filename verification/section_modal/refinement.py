"""Manual read-only mesh-refinement diagnostic, separate from mode truncation."""
import argparse
import csv
from pathlib import Path

import numpy as np

from compare import read_modal, read_solid, relative, resultants, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coarse', type=Path)
    parser.add_argument('fine', type=Path)
    args = parser.parse_args()
    coarse = read_solid(args.coarse / 'transverse_solid.e')
    fine = read_solid(args.fine / 'transverse_refined_solid.e')
    indices = {tuple(np.round(point, 14)): i for i, point in enumerate(fine[0])}
    common = np.array([indices[tuple(np.round(point, 14))] for point in coarse[0]])
    values = dict(coarse_displacement_dofs=coarse[1].size, fine_displacement_dofs=fine[1].size,
                  displacement_relative_change=relative(coarse[1], fine[1][common]),
                  energy_relative_change=abs(coarse[6] / fine[6] - 1), coarse_energy=coarse[6], fine_energy=fine[6])
    for name, model in [('coarse', coarse), ('fine', fine)]:
        xyz, u, points, stress, volume, area, energy, lower, upper = model
        g = resultants(points, stress, area)
        z = np.array(sorted(set(np.round(points[:, :, 2].ravel(), 12))))
        # Moment convention My = -integral(x*sigma_zz dA) on the +z cut.
        expected = -10 * 0.02 * (0.2 - z) ** 2 / 2
        values[name + '_moment_equilibrium_L2'] = relative(g[:, 2], expected)
    with (args.coarse / 'transverse_12_history.csv').open() as stream:
        rows = [row for row in csv.DictReader(stream) if row['kind'] == 'node']
    modal_u = np.array([[float(row[c]) for c in ('ux', 'uy', 'uz')] for row in rows])
    values['modal12_displacement_vs_fine'] = relative(modal_u, fine[1][common])
    values['modal12_energy_vs_fine'] = abs(summary(args.coarse / 'transverse_12_summary.csv')['strain_energy'] / fine[6] - 1)
    _, _, modal_stress, _ = read_modal(args.coarse / 'transverse_12_history.csv', fine[2], args.coarse / 'plate.e')
    values['modal12_axial_stress_vs_fine'] = relative(modal_stress[:, :, 2], fine[3][:, :, 2], fine[4])
    modal_g = resultants(fine[2], modal_stress, fine[5])
    fine_g = resultants(fine[2], fine[3], fine[5])
    values['modal12_moment_vs_fine'] = relative(modal_g[:, 1:], fine_g[:, 1:])
    fine_z = np.array(sorted(set(np.round(fine[2][:, :, 2].ravel(), 12))))
    values['modal12_moment_equilibrium_L2'] = relative(modal_g[:, 2], -10 * 0.02 * (0.2 - fine_z) ** 2 / 2)
    with (args.fine / 'reference_refinement.csv').open('w') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['metric', 'value'])
        writer.writerows(values.items())
    for key, value in values.items():
        print(f'{key}={value:.12g}')
