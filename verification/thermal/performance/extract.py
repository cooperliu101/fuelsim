"""Run with Abaqus Python; save all native frames, nodes and material points."""
import sys
import numpy as np
from odbAccess import openOdb

odb = openOdb(sys.argv[1], readOnly=True)


def nodal(field):
    blocks = field.bulkDataBlocks
    labels = np.concatenate([np.array(b.nodeLabels,copy=True) for b in blocks])
    values = np.concatenate([np.array(b.data,copy=True) for b in blocks]).ravel()
    order = np.argsort(labels)
    assert np.array_equal(labels[order],np.arange(1,len(labels)+1))
    return values[order]


def points(field):
    blocks = field.bulkDataBlocks
    elements = np.concatenate([np.array(b.elementLabels,copy=True) for b in blocks])
    numbers = np.concatenate([np.array(b.integrationPoints,copy=True) for b in blocks])
    values = np.concatenate([np.array(b.data,copy=True) for b in blocks])
    order = np.lexsort((numbers,elements))
    assert np.array_equal(elements[order],np.repeat(np.arange(1,len(elements)//8+1),8))
    assert np.array_equal(numbers[order],np.tile(np.arange(1,9),len(elements)//8))
    return values[order].reshape((-1,8,3))


times, temperatures, reactions, fluxes = [], [], [], []
coordinates = None
for step in odb.steps.values():
    for frame in step.frames:
        if frame.frameValue <= 0:
            continue
        nt = nodal(frame.fieldOutputs['NT11'])
        r = nodal(frame.fieldOutputs['RFL11'])
        q = points(frame.fieldOutputs['HFL'])
        xyz = points(frame.fieldOutputs['COORD'])
        if coordinates is None:
            coordinates = xyz
        else:
            assert np.array_equal(coordinates,xyz)
        times.append(frame.frameValue)
        temperatures.append(nt)
        reactions.append(r)
        fluxes.append(q)
np.savez_compressed(sys.argv[2], time=times, temperature=temperatures,
                    reaction=reactions, flux=fluxes, coordinates=coordinates)
odb.close()
