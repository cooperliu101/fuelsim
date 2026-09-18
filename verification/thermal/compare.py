"""Compare production Exodus fields with tracked native Abaqus output, at every saved time and point."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import netCDF4
import numpy as np


def metric(value, reference, zero=False):
    value, reference = np.asarray(value), np.asarray(reference)
    difference = value-reference
    if zero:
        return {'absolute_maximum':float(np.max(np.abs(difference))),
                'passed':bool(np.max(np.abs(difference)) < 1e-7)}
    nonzero = reference != 0
    l2 = float(np.linalg.norm(difference)/np.linalg.norm(reference))
    peak = float(abs(np.max(np.abs(value))-np.max(np.abs(reference)))/np.max(np.abs(reference)))
    pointwise = float(np.max(np.abs(difference[nonzero]/reference[nonzero])))
    zeros_pass = bool(np.all(np.abs(difference[~nonzero]) < 1e-7))
    return {'relative_l2':l2, 'relative_peak':peak, 'maximum_pointwise_relative':pointwise,
            'passed':bool(max(l2,peak,pointwise) < 0.005 and zeros_pass)}


def compare(directory, name, references):
    nodes = list(csv.DictReader((references/(name+'.nodes.csv')).open()))
    points = list(csv.DictReader((references/(name+'.points.csv')).open()))
    with netCDF4.Dataset(directory/(name+'_results.e')) as data:
        times = np.asarray(data['time_whole'][:])
        nodal_names = netCDF4.chartostring(data['name_nod_var'][:]).tolist()
        element_names = netCDF4.chartostring(data['name_elem_var'][:]).tolist()
        node_variables = {n:np.asarray(data[f'vals_nod_var{i+1}'][:]) for i,n in enumerate(nodal_names)}
        element_variables = {n:np.concatenate([np.asarray(data[f'vals_elem_var{i+1}eb{b+1}'][:]) for b in range(len(data.dimensions['num_el_blk']))],axis=1) for i,n in enumerate(element_names)}
        def frame(time):
            matches = np.flatnonzero(np.isclose(times,time,rtol=0,atol=1e-10))
            if len(matches) != 1:
                raise AssertionError('Missing or duplicate output time '+str(time))
            return matches[0]
        t, tr, r, rr = [],[],[],[]
        for row in nodes:
            index, node = frame(float(row['time'])),int(row['node'])-1
            t.append(node_variables['temperature'][index,node]);tr.append(float(row['temperature']))
            r.append(node_variables['heat_reaction'][index,node]);rr.append(float(row['reaction']))
        if len(nodes) != len(set(float(row['time']) for row in nodes))*len(data.dimensions['num_nodes']):
            raise AssertionError('Native node coverage is incomplete')
        q, qr = [[],[],[]], [[],[],[]]
        count = len(element_names)//6
        used = set()
        for row in points:
            index, element = frame(float(row['time'])),int(row['element'])-1
            positions = np.array([[element_variables[f'point_{c}_q{j}'][index,element] for c in 'xyz'] for j in range(count)])
            target = np.array([float(row[c]) for c in 'xyz'])
            distances = np.linalg.norm(positions-target,axis=1)
            point = int(np.argmin(distances))
            if distances[point] > 1e-8 or (index,element,point) in used:
                raise AssertionError('Material point correspondence is not unique')
            used.add((index,element,point))
            for d,c in enumerate('xyz'):
                q[d].append(element_variables[f'heat_flux_{c}_q{point}'][index,element])
                qr[d].append(float(row['q'+c]))
        if len(used) != len(set(float(row['time']) for row in points))*len(data.dimensions['num_elem'])*count:
            raise AssertionError('Native integration point coverage is incomplete')
    fields = {'temperature':metric(t,tr), 'temperature_rise':metric(np.array(t)-300,np.array(tr)-300),
              'heat_reaction':metric(r,rr,zero=not np.any(rr)),
              'heat_flux_x':metric(q[0],qr[0],zero=name.endswith('_capacity')), 'heat_flux_y':metric(q[1],qr[1],zero=True),
              'heat_flux_z':metric(q[2],qr[2],zero=True)}
    return {'case':name,'nodes_compared':len(nodes),'points_compared':len(points),
            'passed':all(v['passed'] for v in fields.values()),'fields':fields}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fuelsim',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--diagnostic',action='store_true')
    args=parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    kind=args.case if args.case == 'dcax8_contact' or args.case.endswith('_capacity') else args.case.split('_')[0]
    for filename in [args.case+'.fsi',kind+'.e']:
        shutil.copyfile(args.source/filename,args.work/filename)
    card=args.work/(args.case+'.fsi')
    if card.read_bytes() != (args.source/card.name).read_bytes():
        raise AssertionError('Production card changed during copying')
    run=subprocess.run([str(args.fuelsim.resolve()),'-i',str(card.resolve())],capture_output=True,text=True)
    (args.work/'fuelsim.log').write_text(run.stdout+run.stderr)
    if run.returncode:
        raise RuntimeError(run.stdout+run.stderr)
    result=compare(args.work,args.case,args.source)
    (args.work/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not args.diagnostic and not result['passed']:
        raise AssertionError('Native Abaqus comparison failed; original 0.5% gates retained')
    if args.diagnostic:
        print('Diagnostic comparison only: native agreement is not asserted.')

if __name__ == '__main__':
    main()
