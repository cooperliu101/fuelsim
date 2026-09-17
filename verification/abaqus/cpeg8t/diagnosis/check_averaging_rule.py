"""Identify the contact sampling and primary transfer independently of Fuelsim.

The native fine scan and disconnected-cell probes were not used to tune the
original failed production cases. All positions and weights below are computed
from equal-weight polynomial quadrature and linear test functions.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def rule():
    positions, weights = [], []
    for a, b, n in [(0., 1/6, 3), (1/6, 5/6, 5), (5/6, 1., 3)]:
        positions.extend((a+b)/2 + (b-a)*(np.arange(n)-(n-1)/2)/np.sqrt(n*n-1))
        weights.extend([(b-a)/n]*n)
    x, w = np.array(positions), np.array(weights)
    xi = 2*x-1
    test = np.array([np.maximum(-xi,0),np.maximum(xi,0),1-abs(xi)])
    shape = np.array([xi*(xi-1)/2,xi*(xi+1)/2,1-xi*xi]).T
    return x, w, test, shape


def compact_native():
    result = {}
    for name in ('contact_transfer_fine','contact_primary_cells'):
        raw = (ROOT/(name+'_fields.json')).read_bytes()
        result[name+'_sha256'] = hashlib.sha256(raw).hexdigest()
        fields = json.loads(raw)['steps']
        if name == 'contact_transfer_fine':
            result['scan'] = []
            for frame in fields['SCAN'][1:]:
                rf = {v['nodeLabel']:v['data'][1] for v in frame['fields']['RF']['values']}
                result['scan'].append([frame['time']] + [rf[n] for n in (4,7,3,13,10)])
        else:
            result['cells'] = []
            for step in range(65):
                frame = fields[f'CELL{step}'][-1]
                gap = next(v for k,v in frame['fields'].items() if k.startswith('COPEN'))
                values = {v['nodeLabel']:v['data'] for v in gap['values']}
                result['cells'].append([values[n] for n in (513,514,517)])
    raw = (ROOT/'contact_operator_skew_fields.json').read_bytes()
    result['contact_operator_skew_sha256'] = hashlib.sha256(raw).hexdigest()
    fields = json.loads(raw)['steps']
    columns = []
    for node in (9,10,13):
        values = []
        for sign in ('POS','NEG'):
            frame = fields[f'N{node}_{sign}'][-1]
            field = next(v for k,v in frame['fields'].items() if k.startswith('COPEN'))
            gap = {v['nodeLabel']:v['data'] for v in field['values']}
            values.append(np.array([gap[n] for n in (9,10,13)]))
        columns.append((values[0]-values[1])/2e-10)
    result['skew_averaging'] = np.array(columns).T.tolist()
    (ROOT/'averaging_reference.json').write_text(json.dumps(result,indent=2)+'\n')


def check():
    for line in (ROOT/'averaging_reference.sha256').read_text().splitlines():
        digest,name = line.split(maxsplit=1)
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Averaging reference checksum mismatch: '+name)
    ref = json.loads((ROOT/'averaging_reference.json').read_text())
    x,w,test,shape = rule()
    area = test@w
    averaging = (test*w)@shape/area[:,None]
    stiffness = averaging.T@np.diag(2*area)@averaging
    original = json.loads((ROOT/'diagnosis.json').read_text())
    np.testing.assert_allclose(averaging,
        original['averaged_constraint_reconstruction']['gap_averaging_matrix'],rtol=0,atol=1e-9)
    np.testing.assert_allclose(stiffness,
        original['shallow_contact_operator']['native'],rtol=0,atol=1e-9)
    skew_weight = w*(1-.4*(2*x-1))
    skew_averaging = (test*skew_weight)@shape/(test@skew_weight)[:,None]
    np.testing.assert_allclose(skew_averaging,ref['skew_averaging'],rtol=0,atol=1e-9)
    radius = (1-1/np.sqrt(2))/12
    cells = np.array(ref['cells'])
    measured = -(cells[1:]-cells[0]).T/1e-7
    prediction = np.zeros((3,64))
    for cell in range(64):
        left = -.5+cell/32
        fraction = np.maximum(0,np.minimum(x+radius,left+1/32)-np.maximum(x-radius,left))/(2*radius)
        prediction[:,cell] = test@(w*fraction)/area
    np.testing.assert_allclose(prediction,measured,rtol=0,atol=1e-9)
    native = np.array(ref['scan'])
    # ODB frame times are single precision; the fixed input increment supplies
    # the exact prescribed displacement. Check that no frames were skipped.
    times = np.arange(1,2001)/2000
    np.testing.assert_allclose(native[:,0],times,rtol=0,atol=3e-8)
    predicted_force = np.zeros((2000,5))
    for frame,time in enumerate(times):
        position = -.01+.01*x+.01*time
        physical_radius = .01*radius
        for center,indices in [(-.01,[0,2,1]),(.01,[2,4,3])]:
            left, right = center-.01, center+.01
            fraction = np.maximum(0,np.minimum(position+physical_radius,right)
                                  -np.maximum(position-physical_radius,left))/(2*physical_radius)
            xi = (position-center)/.01
            primary = np.array([xi*(xi-1)/2,xi*(xi+1)/2,1-xi*xi]).T
            predicted_force[frame,indices] += 100*(w*fraction)@primary
    np.testing.assert_allclose(predicted_force,native[:,1:],rtol=0,atol=1e-8)
    report = dict(gap_averaging_matrix=averaging.tolist(),normalized_stiffness=stiffness.tolist(),
                  sampling_positions=x.tolist(),sampling_weights=w.tolist(),
                  smoothing_half_width_fraction=radius,
                  skew_averaging_maximum_absolute=float(np.max(abs(skew_averaging-ref['skew_averaging']))),
                  cell_transfer_maximum_absolute=float(np.max(abs(prediction-measured))),
                  sliding_force_maximum_absolute_N=float(np.max(abs(predicted_force-native[:,1:]))),
                  scope='Flat and nonuniformly parameterized edges; curved and connected edges use production comparisons.')
    (ROOT/'averaging_rule.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract',action='store_true')
    args = parser.parse_args()
    if args.extract:
        compact_native()
    check()
