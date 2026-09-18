"""Identify conductivity sampling using only native temperatures, HFL and the input mesh."""
import csv
import json
from pathlib import Path
import netCDF4
import numpy as np

root = Path(__file__).resolve().parent
case = 'dcax4_nonlinear'
with (root/(case+'.nodes.csv')).open() as stream:
    temperatures = {(float(r['time']), int(r['node'])): float(r['temperature'])
                    for r in csv.DictReader(stream)}
with netCDF4.Dataset(root/'dcax4.e') as mesh:
    coordinates = np.array([mesh['coordx'][:], mesh['coordy'][:]]).T
    connectivity = np.array(mesh['connect1'][:], dtype=int)
signs = np.array([[-1,-1], [1,-1], [1,1], [-1,1]])
errors = {'gauss_temperature': [], 'paired_corner_temperature': [], 'element_mean_temperature': []}
with (root/(case+'.points.csv')).open() as stream:
    for row in csv.DictReader(stream):
        nodes = connectivity[int(row['element'])-1]
        xy = coordinates[nodes-1]
        t = np.array([temperatures[(float(row['time']), n)] for n in nodes])
        position = np.array([float(row['x']), float(row['y'])])
        natural = 2*(position-xy[0])/(xy[2]-xy[0])-1
        shape = np.prod(1+signs*natural, axis=1)/4
        derivative = np.array([[s[0]*(1+s[1]*natural[1]), s[1]*(1+s[0]*natural[0])]
                               for s in signs])/4
        gradient = t@(derivative@np.linalg.inv(xy.T@derivative))
        native_flux = np.array([float(row['qx']), float(row['qy'])])
        corner = np.argmax(signs@natural)
        samples = {'gauss_temperature': shape@t, 'paired_corner_temperature': t[corner],
                   'element_mean_temperature': t.mean()}
        for name, temperature in samples.items():
            reconstructed = -(10+.05*(temperature-300))*gradient
            errors[name].append(float(np.linalg.norm(native_flux-reconstructed)/np.linalg.norm(native_flux)))
report = {name: {'points': len(values), 'maximum_relative_flux_error': max(values)}
          for name, values in errors.items()}
print(json.dumps(report, indent=2))
if (report['paired_corner_temperature']['maximum_relative_flux_error'] >= 1e-5
        or report['gauss_temperature']['maximum_relative_flux_error'] <= .02):
    raise AssertionError('Native DCAX4 conductivity sampling is no longer identified by paired corner temperature')
