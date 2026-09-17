"""Check production output for a curved secondary surface with a shared node."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from compare import compare, verify_references
from check_restart import equivalent


def check(directory, source):
    fields = json.loads((source/'shared_contact_fields.json').read_text())['steps']['CONTACT'][-1]['fields']
    report = {}
    contact_nodes = (3,4,7,10,13,14,15,18,22,24)
    with netCDF4.Dataset(directory/'shared_contact_results.e') as d:
        names = list(netCDF4.chartostring(d['name_nod_var'][:]))
        forces = []
        for component,field in enumerate(('reaction_x','reaction_y')):
            actual = np.asarray(d[f'vals_nod_var{names.index(field)+1}'][-1])
            reference = {v['nodeLabel']:v['data'][component] for v in fields['RF']['values']}
            zeros = [n not in contact_nodes or (component==0 and n in (3,15)) for n in range(1,27)]
            compare(actual,[reference[n] for n in range(1,27)],zeros,1e-8,field,report)
            np.testing.assert_allclose(sum(actual),0,rtol=0,atol=1e-8)
            forces.append(actual)
        for field,value in [('displacement_x',np.zeros(26)),
                            ('displacement_y',np.r_[np.zeros(13),np.full(13,-.0003)])]:
            actual = np.asarray(d[f'vals_nod_var{names.index(field)+1}'][-1])
            np.testing.assert_allclose(actual,value,rtol=0,atol=1e-12)
        temperature = np.asarray(d[f'vals_nod_var{names.index("temperature")+1}'][-1])
        assert np.isfinite(temperature).sum()==12
        np.testing.assert_allclose(temperature[np.isfinite(temperature)],300,rtol=0,atol=1e-12)
        heat = np.asarray(d[f'vals_nod_var{names.index("heat_reaction")+1}'][-1])
        np.testing.assert_allclose(heat[np.isfinite(heat)],0,rtol=0,atol=1e-8)
        x=np.asarray(d['coordx'][:])
        y=np.asarray(d['coordy'][:])+np.r_[np.zeros(13),np.full(13,-.0003)]
        np.testing.assert_allclose(x@forces[1]-y@forces[0],0,rtol=0,atol=1e-10)
        variables = list(netCDF4.chartostring(d['name_elem_var'][:]))
        for block in (1,2):
            for component in ('xx','yy','zz','xy'):
                for q in range(9):
                    data = np.asarray(d[f'vals_elem_var{variables.index(f"stress_{component}_q{q}")+1}eb{block}'][-1])
                    np.testing.assert_allclose(data,0,rtol=0,atol=1e-7)
        globals_ = dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][-1]))
        gap_field = next(v for k,v in fields.items() if k.startswith('COPEN'))
        native_gap = {v['nodeLabel']:v['data'] for v in gap_field['values']}
        compare([globals_[f'contact_0_q{q}_gap'] for q in range(5)],
                [native_gap[n] for n in (14,15,18,22,24)],False,1e-12,'averaged_gap',report)
    (directory/'shared_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Shared-node curved contact: all nodal reactions, constraint gaps and zero bulk stresses passed')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--mpiexec',type=Path)
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    verify_references(source,'extended_reference.sha256')
    for mode in (['serial','parallel'] if args.mpiexec else ['serial']):
        directory=args.work/mode
        directory.mkdir(parents=True,exist_ok=True)
        for name in ('shared_contact.fsi','shared_contact.e'):
            shutil.copyfile(source/name,directory/name)
            assert (source/name).read_bytes()==(directory/name).read_bytes()
        command=[str(args.executable.resolve()),'-i','shared_contact.fsi']
        if mode=='parallel':
            command=[str(args.mpiexec),'-n','2']+command
        run=subprocess.run(command,cwd=directory,text=True,capture_output=True)
        (directory/'production.log').write_text(run.stdout+run.stderr)
        if run.returncode or 'completed=true' not in run.stdout:
            raise RuntimeError(run.stdout+run.stderr)
        check(directory,source)
    if args.mpiexec:
        equivalent(args.work/'parallel/shared_contact_results.e',args.work/'serial/shared_contact_results.e')
