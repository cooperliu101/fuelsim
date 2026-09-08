from __future__ import print_function
import csv
import sys
from odbAccess import openOdb


def data(value):
    try:
        return value.dataDouble
    except Exception:
        return value.data


odb=openOdb(sys.argv[1],readOnly=True)
try:
    with open(sys.argv[2],'wb') as stream:
        writer=csv.writer(stream)
        writer.writerow(['time','node','t1x','t1y','t1z','t2x','t2y','t2z','slip1','slip2'])
        elapsed=0.
        for step in odb.steps.values():
            for frame_index in range(1,len(step.frames)):
                frame=step.frames[frame_index]
                fields={}
                for name in ['CTANDIR1','CTANDIR2','CSLIP1','CSLIP2']:
                    keys=[k for k in frame.fieldOutputs.keys() if k.strip()==name or k.strip().startswith(name+' ')]
                    if len(keys)!=1: raise RuntimeError('Expected one contact frame field: '+name)
                    fields[name]=dict((v.nodeLabel,data(v)) for v in frame.fieldOutputs[keys[0]].values)
                for node in sorted(fields['CTANDIR1']):
                    writer.writerow([elapsed+frame.frameValue,node]+list(fields['CTANDIR1'][node])+
                                    list(fields['CTANDIR2'][node])+[fields['CSLIP1'][node],fields['CSLIP2'][node]])
            elapsed+=step.timePeriod
finally:
    odb.close()
