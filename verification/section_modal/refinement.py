"""Read-only axial mesh check for the homogeneous affine HEX20 plate reference.

Stress interpolation is exact within each coarse element (degree <= 2 in each
coordinate), and retains every native integration point of the finer solid.
The separate nodal change metric explicitly uses all shared coarse nodes.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from compare import lagrange, read_solid, relative, resultants


def check(coarse_path, fine_path, output):
    coarse, fine = read_solid(coarse_path), read_solid(fine_path)
    ids = {tuple(np.round(point, 14)): i for i, point in enumerate(fine[0])}
    common = np.array([ids[tuple(np.round(point, 14))] for point in coarse[0]])
    target = fine[2].reshape(-1, 3)
    interpolated = np.empty((len(target), 6))
    coverage = np.zeros(len(target), dtype=int)
    axial = {}
    g, _ = np.polynomial.legendre.leggauss(3)
    for e, (lo, hi) in enumerate(zip(coarse[7], coarse[8])):
        interval = (lo[2], hi[2])
        if interval not in axial:
            axial[interval] = np.flatnonzero((target[:, 2] > lo[2]) & (target[:, 2] < hi[2]))
        candidates = axial[interval]
        selected = candidates[np.all((target[candidates, :2] > lo[:2]) & (target[candidates, :2] < hi[:2]), axis=1)]
        natural = 2 * (target[selected] - lo) / (hi - lo) - 1
        interpolated[selected] = np.einsum('pi,pj,pk,kjic->pc',
            lagrange(g, natural[:, 0]), lagrange(g, natural[:, 1]), lagrange(g, natural[:, 2]),
            coarse[3][e].reshape(3, 3, 3, 6), optimize=True)
        coverage[selected] += 1
    if not np.all(coverage == 1):
        raise ValueError('Refinement comparison missed or duplicated a native reference point')
    stress = interpolated.reshape(fine[3].shape)
    values = dict(coarse_displacement_dofs=coarse[1].size, fine_displacement_dofs=fine[1].size,
                  shared_displacement_nodes=len(common), reference_stress_points=len(target),
                  shared_displacement_L2=relative(coarse[1], fine[1][common]),
                  axial_stress_L2=relative(stress[:, :, 2], fine[3][:, :, 2], fine[4]),
                  stress_tensor_L2=relative(stress, fine[3], fine[4][:, :, None] * [1, 1, 1, 2, 2, 2]),
                  bending_moment_L2=relative(resultants(fine[2], stress, fine[5])[:, 1:],
                                             resultants(fine[2], fine[3], fine[5])[:, 1:]),
                  energy_relative=abs(coarse[6] / fine[6] - 1))
    for name, model in [('coarse', coarse), ('fine', fine)]:
        z = np.unique(np.round(model[2][:, :, 2], 12))
        expected = -0.2 * (0.2 - z) ** 2 / 2
        values[name + '_moment_equilibrium_L2'] = relative(resultants(model[2], model[3], model[5])[:, 2], expected)
    with output.open('w') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['metric', 'value'])
        writer.writerows(values.items())
    for key, value in values.items():
        print(f'{key}={value:.12g}')
    return values


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coarse', type=Path)
    parser.add_argument('fine', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    check(args.coarse, args.fine, args.output)
