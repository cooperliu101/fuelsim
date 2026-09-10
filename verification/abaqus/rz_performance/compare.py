"""Read-only final-state accuracy comparison; never links or runs a solver."""
import argparse
import csv
import gzip
import json
import os
from contextlib import contextmanager
from pathlib import Path
import numpy as np


@contextmanager
def result_database(path):
    # Read the existing NetCDF C library directly; no Python NetCDF package is required.
    import ctypes as c
    import ctypes.util
    library=os.environ.get('NETCDF_LIBRARY') or ctypes.util.find_library('netcdf')
    if not library:raise RuntimeError('Set NETCDF_LIBRARY to the Exodus environment libnetcdf shared library')
    lib=c.CDLL(library)
    def call(name,*args):
        code=getattr(lib,name)(*args)
        if code:raise RuntimeError('NetCDF error %d from %s'%(code,name))
    fid=c.c_int();call('nc_open',str(path).encode(),0,c.byref(fid))
    try:
        nd,nv=c.c_int(),c.c_int()
        call('nc_inq',fid,c.byref(nd),c.byref(nv),None,None)
        dimensions={};sizes=[]
        for i in range(nd.value):
            name=c.create_string_buffer(257);size=c.c_size_t()
            call('nc_inq_dim',fid,i,name,c.byref(size));dimensions[name.value.decode()]=size.value;sizes.append(size.value)
        variables={}
        for i in range(nv.value):
            name=c.create_string_buffer(257);kind=c.c_int();rank=c.c_int();ids=(c.c_int*1024)()
            call('nc_inq_var',fid,i,name,c.byref(kind),c.byref(rank),ids,None)
            shape=tuple(sizes[ids[j]] for j in range(rank.value))
            arr=np.empty(shape,dtype='S1' if kind.value==2 else np.float64)
            call('nc_get_var_text' if kind.value==2 else 'nc_get_var_double',fid,i,c.c_void_p(arr.ctypes.data))
            variables[name.value.decode()]=arr
        from types import SimpleNamespace
        yield SimpleNamespace(dimensions=dimensions,variables=variables)
    finally:call('nc_close',fid)


def names(variable):
    return [b''.join(row).decode().rstrip('\0 ') for row in variable[:]]


def metric(actual, reference, zero_tolerance):
    a,b=np.asarray(actual,dtype=float),np.asarray(reference,dtype=float)
    if a.shape != b.shape or not a.size or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Missing, nonfinite or mismatched field')
    if a.ndim == 1:a,b=a[:,None],b[:,None]
    difference=np.linalg.norm(a-b,axis=1);scale=np.linalg.norm(b,axis=1)
    nonzero=scale>0
    result={'samples':len(scale),'maximum_absolute_error':float(difference.max()),
            'zero_reference_samples':int((~nonzero).sum()),
            'zero_reference_absolute_error':float(difference[~nonzero].max(initial=0)),
            'zero_absolute_tolerance':zero_tolerance}
    if nonzero.any():
        result.update(relative_l2=float(np.linalg.norm(difference)/np.linalg.norm(scale)),
                      relative_absolute_peak=float(abs(np.linalg.norm(a,axis=1).max()-scale.max())/scale.max()),
                      relative_maximum_absolute_error=float(difference.max()/scale.max()),
                      maximum_pointwise_relative=float((difference[nonzero]/scale[nonzero]).max()))
    result['passed']=result['zero_reference_absolute_error']<=zero_tolerance and all(
        result.get(k,0)<.0001 for k in ['relative_l2','relative_absolute_peak','maximum_pointwise_relative'])
    return result


def compare(result_path, prefix, element="cax4t"):
    point_count = {"cax4rt": 1, "cax4t": 4, "cax8t": 9}[element]
    def rows(kind):
        with gzip.open(str(prefix)+'_'+kind+'.csv.gz','rt') as f:return list(csv.DictReader(f))
    nodes,points,contacts=rows('nodes'),rows('points'),rows('contact')
    with result_database(result_path) as f:
        if len(f.variables['time_whole'][:]) != 1:raise ValueError('Expected final steady output only')
        nodal={name:np.array(f.variables['vals_nod_var%d'%i][-1]) for i,name in enumerate(names(f.variables['name_nod_var']),1)}
        elem={name:np.concatenate([np.array(f.variables['vals_elem_var%deb%d'%(i,b)][-1]) for b in range(1,f.dimensions['num_el_blk']+1)]) for i,name in enumerate(names(f.variables['name_elem_var']),1)}
        glob=dict(zip(names(f.variables['name_glo_var']),map(float,f.variables['vals_glo_var'][-1])))
        connectivity=np.concatenate([np.asarray(f.variables['connect%d'%b],dtype=int) for b in range(1,f.dimensions['num_el_blk']+1)])
    if glob['load_factor'] != 1.0 or not np.all(elem['material_point_count']==point_count):
        raise ValueError('Incomplete load path or wrong element integration')
    thermal_nodes=np.unique(connectivity[:,:4])-1
    if element == 'cax8t':
        active=np.zeros(len(nodal['temperature']));active[thermal_nodes]=1
        if not np.array_equal(nodal['temperature_active'],active):raise ValueError('Temperature DOF coverage differs')
    if not (np.all(nodal['temperature']>=500) and np.all(nodal['temperature']<=2500)
            and all(500<=float(nodes[n]['temperature'])<=2500 for n in thermal_nodes)):
        raise ValueError('Temperature lies outside the sampled conductivity table')
    if [int(r['node']) for r in nodes] != list(range(1,len(nodal['temperature'])+1)):raise ValueError('Node labels/coverage differ')
    if len(points)!=point_count*len(elem['material_point_count']):raise ValueError('Material point coverage differs')
    if any(float(r['time'])!=20 for r in nodes+points+contacts):raise ValueError('Reference is not the final 20-increment state')
    metrics={}
    metrics['temperature']=metric(nodal['temperature'][thermal_nodes],[float(nodes[n]['temperature']) for n in thermal_nodes],1e-8)
    if element == 'cax8t':
        for edge in range(4):
            expected=.5*(nodal['temperature'][connectivity[:,edge]-1]+nodal['temperature'][connectivity[:,(edge+1)%4]-1])
            if not np.array_equal(nodal['temperature'][connectivity[:,edge+4]-1],expected):
                raise ValueError('Midside output must interpolate the temperature DOFs')
    a=np.stack([nodal['displacement_r'],nodal['displacement_z']],axis=1)
    b=np.array([[float(r['ur']),float(r['uz'])] for r in nodes])
    with result_database(result_path) as database:
        constrained=np.stack([database.variables['coordx']==0,database.variables['coordy']==0],axis=1)
    # These components are prescribed exactly zero by both input decks. Check
    # both solvers absolutely, so Abaqus roundoff does not create a false 100%.
    metrics['prescribed_displacement']=metric(np.concatenate([a[constrained],b[constrained]]),np.zeros(2*constrained.sum()),1e-12)
    a[constrained]=0;b[constrained]=0
    metrics['free_displacement']=metric(a,b,1e-12)
    actual,reference=[],[]
    qmap=list(range(9)) if element == "cax8t" else ([0] if element == "cax4rt" else [0,1,3,2])
    for index,r in enumerate(points):
        e,q=divmod(index,point_count)
        if int(r['element'])!=e+1 or int(r['point'])!=q+1:raise ValueError('Material labels differ')
        actual.append([elem['stress_'+c+'_q'+str(qmap[q])][e]*w for c,w in [('rr',1),('zz',1),('hoop',1),('rz',2**.5)]])
        reference.append([float(r['stress_'+c])*w for c,w in [('rr',1),('zz',1),('hoop',1),('rz',2**.5)]])
    metrics['stress_tensor']=metric(actual,reference,1e-3)
    indices=[int(r['node'])-1 for r in contacts]
    if len(set(indices))!=len(indices) or not indices:raise ValueError('Invalid contact node labels')
    expected=np.flatnonzero(np.isfinite(nodal['contact_gap_fuel_cladding']))
    if sorted(indices)!=list(expected):raise ValueError('Contact coverage differs')
    for name,field,column,tol in [('contact_pressure','pressure','pressure',1e-3),('contact_gap','gap','gap',1e-12),('contact_normal_force','normal_force','normal_r',1e-8)]:
        actual_field = 'recovered_pressure' if element == 'cax8t' and field == 'pressure' else field
        metrics[name]=metric(nodal['contact_'+actual_field+'_fuel_cladding'][indices],[(float(r['normal_r'])**2+float(r['normal_z'])**2)**.5 if column=='normal_r' else float(r[column]) for r in contacts],tol)
    metrics['contact_total_force']=metric([glob['contact_force_fuel_cladding']],[sum((float(r['normal_r'])**2+float(r['normal_z'])**2)**.5 for r in contacts)],1e-8)
    return metrics

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('result');p.add_argument('prefix');p.add_argument('--report',required=True);p.add_argument('--element',choices=('cax4t','cax4rt','cax8t'),default='cax4t');args=p.parse_args()
    m=compare(args.result,Path(args.prefix),args.element);Path(args.report).write_text(json.dumps(m,indent=2)+'\n')
    for name,r in m.items():print(name,'PASS' if r['passed'] else 'FAIL', 'relative errors (%)',*[100*r.get(k,0) for k in ['relative_l2','relative_absolute_peak','maximum_pointwise_relative']])
    raise SystemExit(0 if all(r['passed'] for r in m.values()) else 1)
