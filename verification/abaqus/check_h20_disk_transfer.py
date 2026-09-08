"""Test circular averaging of extrapolated primary-face shape functions."""
import math
import numpy as np
from analyze_h20_28_finite_transfer import extract
from diagnose_c3d20rt_contact_gap import rules, shape


def disk_rectangle(x0, x1, y0, y1, radius):
    lo, hi = max(x0, -radius), min(x1, radius)
    if hi <= lo or y1 <= -radius or y0 >= radius:
        return 0.0
    cuts = [lo, hi]
    for y in [y0, y1]:
        if abs(y) < radius:
            x = math.sqrt(radius * radius - y * y)
            cuts.extend(v for v in [-x, x] if lo < v < hi)
    cuts = sorted(set(cuts))
    def integral(x):
        return .5 * (x * math.sqrt(max(0., radius * radius - x * x)) +
                     radius * radius * math.asin(max(-1., min(1., x / radius))))
    area = 0.0
    for a, b in zip(cuts[:-1], cuts[1:]):
        root = math.sqrt(max(0., radius * radius - (.5 * (a + b)) ** 2))
        upper, lower = min(y1, root), max(y0, -root)
        if upper <= lower:
            continue
        coefficient = int(root <= y1) + int(-root >= y0)
        constant = (y1 if y1 < root else 0.) - (y0 if y0 > -root else 0.)
        area += coefficient * (integral(b) - integral(a)) + constant * (b-a)
    return area / (math.pi * radius * radius)


def samples(index):
    result = rules()[0 if index < 4 else 1].copy()
    for row in result:
        x, y = row[:2]
        if index == 1: row[0] = 1-x
        elif index == 2: row[:2] = [1-x, 1-y]
        elif index == 3: row[1] = 1-y
        elif index == 5: row[:2] = [1-y, x]
        elif index == 6: row[1] = 1-y
        elif index == 7: row[:2] = [y, x]
    result[:, 2] /= result[:, 2].sum()
    return result


def predict(grid, coordinates):
    radius = min(
        min(np.min(samples(i)[:, :2]), np.min(1-samples(i)[:, :2])) for i in range(8))
    lookup = {tuple(np.rint(2 * grid * xy).astype(int)): i for i, xy in enumerate(coordinates)}
    result = np.zeros((8, len(coordinates)))
    output_rows = [0, 1, 3, 2, 4, 7, 5, 6]
    for local in range(8):
        for x, y, weight in samples(local):
            total = 0.0
            for i in range(max(0, int((x-radius)*grid)), min(grid, int((x+radius)*grid)+1)):
                for j in range(max(0, int((y-radius)*grid)), min(grid, int((y+radius)*grid)+1)):
                    fraction = disk_rectangle(i/grid-x, (i+1)/grid-x, j/grid-y, (j+1)/grid-y, radius)
                    if fraction <= 0: continue
                    total += fraction
                    n = shape(2*(grid*x-i)-1, 2*(grid*y-j)-1)[0]
                    keys = [(2*i,2*j),(2*i+2,2*j),(2*i+2,2*j+2),(2*i,2*j+2),
                            (2*i+1,2*j),(2*i+2,2*j+1),(2*i+1,2*j+2),(2*i,2*j+1)]
                    for key, value in zip(keys, n):
                        result[output_rows[local], lookup[key]] += weight * fraction * value
            if abs(total-1) > 1e-10: raise ValueError('Disk partition lost area')
    return radius, result


if __name__ == '__main__':
    for grid, stem in [(8, 'h20_28_finite_transfer_shallow_probe'),
                       (16, 'h20_28_finite_transfer_16_probe'),
                       (32, 'h20_28_finite_transfer_32_probe')]:
        reference = extract(stem)
        radius, predicted = predict(grid, reference[3])
        difference = predicted-reference[2]
        print('grid', grid, 'radius', radius, 'relative', np.linalg.norm(difference)/np.linalg.norm(reference[2]),
              'maximum', abs(difference).max())
        print('corner_relative', np.linalg.norm(difference[0])/np.linalg.norm(reference[2][0]),
              'edge_relative', np.linalg.norm(difference[4])/np.linalg.norm(reference[2][4]))
