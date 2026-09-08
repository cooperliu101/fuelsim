"""Read production Exodus and native Abaqus CSVs; compare at t=1,...,20 s.

Manual temporal study only. Never generates or changes physical input cards.
Requires NumPy, Matplotlib and a NetCDF C library in the analysis environment.
"""
import csv
import ctypes as C
import ctypes.util
import gzip
import os
from pathlib import Path
import sys

library = os.environ.get('NETCDF_LIBRARY') or ctypes.util.find_library('netcdf')
if not library:
    raise RuntimeError('Set NETCDF_LIBRARY to the NetCDF C shared library path')
lib = C.CDLL(library)
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'verification/abaqus/b15_time'
LEVELS = [('100', 1.), ('050', .5), ('025', .25), ('0125', .125),
          ('00625', .0625), ('003125', .03125), ('0015625', .015625),
          ('00078125', .0078125)]


class Exodus:
    def __init__(self, path):
        self.id = C.c_int()
        if lib.nc_open(str(path).encode(), 0, C.byref(self.id)):
            raise RuntimeError('Cannot open ' + str(path))

    def get(self, name, text=False):
        var = C.c_int()
        if lib.nc_inq_varid(self.id, name.encode(), C.byref(var)):
            raise RuntimeError('Missing variable ' + name)
        ndim = C.c_int()
        lib.nc_inq_varndims(self.id, var, C.byref(ndim))
        dims = (C.c_int * ndim.value)()
        lib.nc_inq_vardimid(self.id, var, dims)
        shape = []
        for dim in dims:
            size = C.c_size_t()
            lib.nc_inq_dimlen(self.id, dim, C.byref(size))
            shape.append(size.value)
        result = np.empty(shape, dtype='S1' if text else 'float64')
        read = lib.nc_get_var_text if text else lib.nc_get_var_double
        if read(self.id, var, C.c_void_p(result.ctypes.data)):
            raise RuntimeError('Cannot read ' + name)
        return result

    def names(self, name):
        return [b''.join(x).decode().rstrip('\0 ') for x in self.get(name, True)]

    def close(self):
        lib.nc_close(self.id)


def reference_rows(directory, kind):
    path = directory / ('b13_small_cax4t_' + kind + '.csv')
    stream = path.open() if path.exists() else gzip.open(str(path) + '.gz', 'rt')
    with stream:
        for row in csv.DictReader(stream):
            time = float(row['time'])
            if time > 0 and abs(time - round(time)) < 1e-8:
                yield {k: float(v) for k, v in row.items()}


def reference(directory):
    nodes = list(reference_rows(directory, 'nodes'))
    points = list(reference_rows(directory, 'points'))
    contacts = list(reference_rows(directory, 'contact'))
    assert len(nodes) == 20 * 53 and len(points) == 20 * 136
    assert len(contacts) == 20 * 5
    fields = {'temperature': np.array([[r['temperature']] for r in nodes]),
              'displacement': np.array([[r['ur'], r['uz']] for r in nodes])}
    for prefix in ['stress', 'elastic', 'plastic', 'creep']:
        fields[prefix] = np.array([[r[prefix + '_' + c] for c in ['rr', 'zz', 'hoop', 'rz']]
                                   for r in points])
        fields[prefix][:, 3] *= np.sqrt(2.)
    for name in ['equiv_creep', 'equiv_plastic']:
        fields[name] = np.array([[r[name]] for r in points])
    for name in ['pressure', 'gap']:
        fields[name] = np.array([[r[name]] for r in contacts])
    fields['normal_force'] = np.array([[np.hypot(r['normal_r'], r['normal_z'])] for r in contacts])
    fields['normal_resultant'] = fields['normal_force'].reshape(20, 5).sum(axis=1)[:, None]
    return fields, np.array([int(r['node']) - 1 for r in contacts[:5]])


def production(path, contact_nodes):
    file = Exodus(path)
    times = file.get('time_whole')
    take = np.array([int(np.argmin(abs(times - t))) for t in range(1, 21)])
    assert np.max(abs(times[take] - np.arange(1, 21))) < 1e-8
    nodal_names, elem_names = file.names('name_nod_var'), file.names('name_elem_var')

    def nodal(name):
        return file.get('vals_nod_var' + str(nodal_names.index(name) + 1))[take]

    def element(name):
        index = elem_names.index(name) + 1
        return np.concatenate([file.get('vals_elem_var%deb%d' % (index, block))[take]
                               for block in [1, 2]], axis=1)

    fields = {'temperature': nodal('temperature').reshape(-1, 1),
              'displacement': np.stack([nodal('displacement_r'), nodal('displacement_z')], axis=-1).reshape(-1, 2)}
    for prefix in ['stress', 'elastic', 'plastic', 'creep', 'equiv_creep', 'equiv_plastic']:
        components = ['rr', 'zz', 'hoop', 'rz'] if prefix in ['stress', 'elastic', 'plastic', 'creep'] else ['']
        values = np.stack([np.stack([element(prefix + ('_' + c if c else '') + '_q' + str(q))
                                    for c in components], axis=-1) for q in [0, 1, 3, 2]], axis=2)
        fields[prefix] = values.reshape(-1, len(components))
        if len(components) == 4:
            fields[prefix][:, 3] *= np.sqrt(2.)
    for short, name in [('pressure', 'pressure'), ('gap', 'gap'), ('normal_force', 'normal_force')]:
        fields[short] = nodal('contact_' + name + '_fuel_cladding')[:, contact_nodes].reshape(-1, 1)
    fields['normal_resultant'] = fields['normal_force'].reshape(20, 5).sum(axis=1)[:, None]
    coords = np.stack([file.get('coordx'), file.get('coordy')], axis=-1)
    fixed = np.tile(coords == 0., (20, 1))
    assert np.max(abs(fields['displacement'][fixed])) < 1e-12
    fields['displacement'][fixed] = 0.
    file.close()
    return fields, fixed


def metrics(actual, ref):
    d = np.linalg.norm(actual - ref, axis=1)
    r = np.linalg.norm(ref, axis=1)
    nz = r != 0
    tensor = ref.shape[1] > 1
    peak = d.max() if tensor else abs(np.abs(actual).max() - r.max())
    return {'relative_l2_percent': 100 * np.linalg.norm(d) / np.linalg.norm(r),
            'relative_peak_percent': 100 * peak / r.max(),
            'maximum_pointwise_percent': 100 * np.max(d[nz] / r[nz]),
            'maximum_absolute_difference': d.max(),
            'zero_reference_count': int((~nz).sum()),
            'maximum_zero_reference_absolute_difference': d[~nz].max() if (~nz).any() else 0.}


def write(name, rows):
    with (STUDY / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(paired, convergence):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    for field, label in [('stress', 'Stress'), ('pressure', 'Contact pressure'),
                         ('equiv_creep', 'Equivalent creep'), ('equiv_plastic', 'Equivalent plastic')]:
        rows = [r for r in paired if r['field'] == field]
        axes[0].loglog([r['dt'] for r in rows], [r['relative_l2_percent'] for r in rows], 'o-', label=label)
    for solver, style in [('fuelsim', 'o-'), ('abaqus', 's--')]:
        for field in ['stress', 'equiv_creep']:
            rows = [r for r in convergence if r['solver'] == solver and r['field'] == field]
            axes[1].loglog([r['fine_dt'] for r in rows], [r['relative_l2_percent'] for r in rows], style,
                          label=solver + ': ' + field)
    for ax in axes:
        ax.axhline(.1, color='gray', linestyle=':', label='0.1%')
        ax.invert_xaxis()
        ax.set_xlabel('Time step (s), finer to the right')
        ax.set_ylabel('Relative L2 difference (%)')
        ax.grid(True, which='both', alpha=.2)
        ax.legend(fontsize=8)
    axes[0].set_title('Fuelsim versus Abaqus')
    axes[1].set_title('Each solver: dt versus dt/2')
    fig.suptitle('CAX4T small strain; common samples at t = 1, 2, ..., 20 s')
    fig.savefig(STUDY / 'convergence.png', dpi=180)
    plt.close(fig)


def main():
    native_root = Path(sys.argv[1]) if len(sys.argv) > 1 else STUDY
    paired, final, convergence, previous = [], [], [], None
    for tag, dt in LEVELS:
        directory = ROOT / 'verification/abaqus' if tag == '100' else native_root / ('dt' + tag)
        ref, cn = reference(directory)
        path = ROOT / ('build/blackbox/b13_small_cax4t/verification/fuelsim/transient_b13_small_cax4t_results.e'
                       if tag == '100' else 'verification/fuelsim/transient_b15_time_' + tag + '_results.e')
        if not path.exists() and tag != '100':
            path = ROOT / 'build/b15_time' / ('dt' + tag) / path.name
        actual, fixed = production(path, cn)
        assert np.max(abs(ref['displacement'][fixed])) < 1e-12
        ref['displacement'][fixed] = 0.
        for name in actual:
            paired.append(dict(dt=dt, field=name, **metrics(actual[name], ref[name])))
            count = len(actual[name]) // 20
            final.append(dict(dt=dt, field=name, **metrics(actual[name][-count:], ref[name][-count:])))
            if previous:
                for solver, data in [('fuelsim', actual), ('abaqus', ref)]:
                    convergence.append(dict(coarse_dt=previous[0], fine_dt=dt, solver=solver,
                                            field=name, **metrics(previous[1][solver][name], data[name])))
        previous = dt, {'fuelsim': actual, 'abaqus': ref}
    write('common_time_comparison.csv', paired)
    write('end_time_comparison.csv', final)
    write('self_convergence.csv', convergence)
    plot(paired, convergence)
    strict, _ = reference(native_root / 'dt0015625')
    adjusted, _ = reference(native_root / 'dt0015625_heat_tol')
    for fields in [strict, adjusted]:
        assert np.max(abs(fields['displacement'][fixed])) < 1e-12
        fields['displacement'][fixed] = 0.
    write('heat_tolerance_sensitivity.csv',
          [dict(field=name, **metrics(adjusted[name], strict[name])) for name in strict])


if __name__ == '__main__':
    main()
