"""Check transfer invariance on nested Quad8 primary displacement spaces."""
import argparse
import numpy as np
from analyze_h20_28_finite_transfer import extract
from diagnose_c3d20rt_contact_gap import shape


def reference_name(grid):
    if grid == 8:
        return 'h20_28_finite_transfer_probe'
    return f'h20_28_finite_transfer_{grid}_probe'


def compare(coarse_grid, fine_grid):
    coarse = extract(reference_name(coarse_grid))
    fine = extract(reference_name(fine_grid))
    lookup = {tuple(np.rint(2 * coarse_grid * p).astype(int)): i
              for i, p in enumerate(coarse[3])}
    prolongation = np.zeros((len(fine[3]), len(coarse[3])))
    for row, (x, y) in enumerate(fine[3]):
        i = min(coarse_grid - 1, int(x * coarse_grid))
        j = min(coarse_grid - 1, int(y * coarse_grid))
        values = shape(2 * (coarse_grid * x - i) - 1,
                       2 * (coarse_grid * y - j) - 1)[0]
        coordinates = [(2*i, 2*j), (2*i+2, 2*j), (2*i+2, 2*j+2),
                       (2*i, 2*j+2), (2*i+1, 2*j), (2*i+2, 2*j+1),
                       (2*i+1, 2*j+2), (2*i, 2*j+1)]
        for value, point in zip(values, coordinates):
            prolongation[row, lookup[point]] = value
    np.testing.assert_allclose(prolongation.sum(axis=1), 1, atol=1e-14)
    np.testing.assert_allclose(prolongation @ coarse[3], fine[3], atol=1e-14)
    difference = fine[2] @ prolongation - coarse[2]
    print('coarse_grid', coarse_grid, 'fine_grid', fine_grid)
    print('nested_transfer_relative',
          np.linalg.norm(difference) / np.linalg.norm(coarse[2]),
          'max', abs(difference).max())
    for row in [0, 4]:
        print(row, np.linalg.norm(difference[row]) / np.linalg.norm(coarse[2][row]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid', type=int, choices=[8, 16, 32],
                        nargs='+', default=[8, 16])
    grids = parser.parse_args().grid
    if len(grids) < 2 or grids != sorted(set(grids)):
        parser.error('provide at least two strictly increasing grid sizes')
    for coarse_grid, fine_grid in zip(grids, grids[1:]):
        compare(coarse_grid, fine_grid)
