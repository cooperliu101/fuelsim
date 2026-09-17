"""Native comparisons for added cases; write every failed metric before failing."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import netCDF4
import numpy as np

from compare import compare, verify_references


def check(case, directory, source):
    native=json.loads((source/(case+'_fields.json')).read_text())
    step=list(native['steps'])[-1]
    fields=native['steps'][step][-1]['fields']
    contact=case=='sections_contact'
    thermal=case=='thermal_boundary'
    rotation=case=='bending_rotation'
    source_probe=case=='source_nonaffine'
    count=16 if contact else 8
    corners=[1,2,3,4,9,10,11,12] if contact else [1,2,3,4]
    report, failures = {}, []
    def metric(a,r,z,absolute,label):
        try:
            compare(a,r,z,absolute,label,report)
        except AssertionError as error:
            failures.append(str(error))
            report[label]={'status':'failed','reason':str(error),
                           'actual':np.asarray(a).tolist(),'reference':np.asarray(r).tolist()}
    with netCDF4.Dataset(directory/(case+'_results.e')) as d:
        nodal=list(netCDF4.chartostring(d['name_nod_var'][:]))
        elements=list(netCDF4.chartostring(d['name_elem_var'][:]))
        globals_=dict(zip(netCDF4.chartostring(d['name_glo_var'][:]),d['vals_glo_var'][-1]))
        for name,key,component,nodes,zeros,absolute in [
            ('temperature','NT11',None,corners,lambda n:False,1e-9),
            ('displacement_x','U',0,range(1,count+1),lambda n:(source_probe and n not in (6,8))
             or (not source_probe and not rotation and (contact or thermal or n in (6,8))),1e-12),
            ('displacement_y','U',1,range(1,count+1),lambda n:not rotation and (not contact or n in (1,2,5)),1e-12),
            ('reaction_x','RF',0,range(1,count+1),lambda n:False,1e-8),
            ('reaction_y','RF',1,range(1,count+1),lambda n:(contact and n not in (1,2,5,11,12,15))
             or (rotation and n in (5,7)) or (source_probe and n in (6,8)),1e-8),
            ('heat_reaction','RFL11',None,corners,lambda n:(contact and n not in (1,2,11,12))
             or (thermal and n in (3,4)) or rotation,1e-8),
        ]:
            data={v['nodeLabel']:v['data'] for v in fields[key]['values']}
            reference=[data[n] if component is None else data[n][component] for n in nodes]
            actual=np.asarray(d[f'vals_nod_var{nodal.index(name)+1}'][-1])[np.array(list(nodes))-1]
            metric(actual,reference,[zeros(n) for n in nodes],absolute,name)
        for component,name in enumerate(['xx','yy','zz','xy']):
            actual,reference=[],[]
            for row in fields['S']['values']:
                variable=elements.index(f'stress_{name}_q{row["integrationPoint"]-1}')+1
                actual.append(d[f'vals_elem_var{variable}eb{row["elementLabel"]}'][-1,0])
                reference.append(row['data'][component])
            # The prescribed bending map has diagonal symmetric midpoint strain;
            # the subsequent 90-degree rotation preserves zero shear stress.
            zeros=([v['integrationPoint'] in (4,5,6) for v in fields['S']['values']]
                   if source_probe and name=='xy' else (thermal or rotation) and name=='xy')
            metric(actual,reference,zeros,1e-7,'stress_'+name)
        for section,node in ([('lower',17),('upper',18)] if contact else [('solid',100)]):
            region=next(r for r in native['history'][step] if r.startswith('Node ') and r.endswith('.'+str(node)))
            h={key:values[-1][1] for key,values in native['history'][step][region].items()}
            prefix='section_'+section+'_'
            x,y=globals_[prefix+'origin_x'],globals_[prefix+'origin_y']
            for name,key in [('u3','U3'),('rotation_x','UR1'),('rotation_y','UR2'),
                             ('axial_force','RF3'),('moment_x','RM1'),('moment_y','RM2')]:
                value=h[key]
                if key=='U3': value+=y*h['UR1']-x*h['UR2']
                if key=='RM1': value-=y*h['RF3']
                if key=='RM2': value+=x*h['RF3']
                metric([globals_[prefix+name]],[value],(contact and section=='upper' and key=='RF3')
                       or (thermal and key in ('U3','UR1','UR2')) or (rotation and key in ('UR1','UR2','RM1'))
                       or (source_probe and key in ('U3','UR1','UR2','RM1')),
                       1e-9,prefix+name)
        if contact:
            force_key=next(k for k in fields if k.startswith('CNORMF'))
            force=sum(v['data'][1] for v in fields[force_key]['values'] if v['nodeLabel'] in (9,10,13))
            heat=sum(v['data'] for v in fields['RFL11']['values'] if v['nodeLabel'] in (11,12))
            total_force = sum(value for name, value in globals_.items()
                              if name.startswith('contact_0_q') and name.endswith('_force'))
            total_heat = sum(value for name, value in globals_.items()
                             if name.startswith('contact_0_q') and name.endswith('_heat_rate'))
            metric([total_force],[force],False,1e-8,'contact_force')
            metric([total_heat],[heat],False,1e-8,'contact_heat_rate')
    result={'status':'failed' if failures else 'passed','metrics':report,'failures':failures}
    (directory/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    if failures:
        raise AssertionError('\n'.join(failures))
    print(case+': all final nodes, stresses and section controls passed native comparison')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',choices=['sections_contact','shear_distorted','thermal_boundary','bending_rotation',
                                        'source_nonaffine'],required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    verify_references(source, 'extended_reference.sha256')
    args.work.mkdir(parents=True,exist_ok=True)
    mesh={'sections_contact':'contact.e','shear_distorted':'distorted.e','thermal_boundary':'boundary.e',
          'bending_rotation':'boundary.e','source_nonaffine':'boundary.e'}[args.case]
    for name in [args.case+'.fsi',mesh]:
        shutil.copyfile(source/name,args.work/name)
        assert (source/name).read_bytes()==(args.work/name).read_bytes()
    run=subprocess.run([str(args.executable.resolve()),'-i',args.case+'.fsi'],
                       cwd=args.work,capture_output=True,text=True)
    (args.work/'production.log').write_text(run.stdout+run.stderr)
    if run.returncode or 'completed=true' not in run.stdout:
        raise RuntimeError(run.stdout+run.stderr)
    check(args.case,args.work,source)
