"""Manual contact refinement meshes; no input cards are generated."""
from pathlib import Path
import numpy as np
from create_extended_meshes import mesh

ROOT = Path(__file__).resolve().parent


def create(count):
    coordinates, blocks = [], {}
    sets = {'field': [], 'lower': [], 'top': []}
    for label, bounds in [('lower', [-.03, 0, .03]), ('upper', np.linspace(-.012,.012,count+1))]:
        mapping, elements = {}, []
        for left, right in zip(bounds[:-1], bounds[1:]):
            element = []
            for x, level in [(left,0),(right,0),(right,1),(left,1),((left+right)/2,0),
                             (right,.5),((left+right)/2,1),(left,.5)]:
                key = (round(float(x),14), level)
                if key not in mapping:
                    n = len(coordinates)+1
                    mapping[key] = n
                    y = -.01+.01*level if label == 'lower' else .0001+2*x*x+.01*level
                    coordinates.append((x,y))
                    sets['field'].append(n)
                    if label == 'lower':
                        sets['lower'].append(n)
                    if label == 'upper' and level == 1:
                        sets['top'].append(n)
                element.append(mapping[key])
            elements.append(element)
        blocks[label] = elements
    mesh(ROOT/f'contact_refine_{count}.e',coordinates,blocks,sets,
         {'primary': [(1,3),(2,3)], 'secondary': [(i+3,1) for i in range(count)]})


if __name__ == '__main__':
    for count in (2,4,8,16):
        create(count)
