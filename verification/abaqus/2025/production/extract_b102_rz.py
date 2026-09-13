from __future__ import print_function
import csv, sys
from odbAccess import openOdb

def data(value):
    try: return value.dataDouble
    except Exception: return value.data

def values(frame,name,nodal=False):
    try: field=frame.fieldOutputs[name]
    except Exception: return {}
    result={}
    for value in field.values:
        key=value.nodeLabel if nodal else (value.elementLabel,value.integrationPoint)
        result[key]=data(value)
    return result

odb=openOdb(path=sys.argv[1],readOnly=True)
try:
    with open(sys.argv[2]+"_nodes.csv","w", newline="") as nf,open(sys.argv[2]+"_points.csv","w", newline="") as pf:
        nodes,points=csv.writer(nf),csv.writer(pf)
        nodes.writerow(["time","node","temperature","ur","uz","rf_r","rf_z","reaction_heat","contact_pressure","contact_gap","contact_force_z","contact_heat_flux"])
        points.writerow(["time","element","point"]+[p+c for p in ["stress_","elastic_","plastic_","creep_"] for c in ["rr","zz","hoop","rz"]]+["equiv_plastic","equiv_creep"])
        elapsed=0
        for step in odb.steps.values():
            for frame_index in range(1,len(step.frames)):
                frame=step.frames[frame_index];time=elapsed+frame.frameValue
                nflds={n:values(frame,n,True) for n in ["NT11","U","RF","RFL11"]}
                cfields={}
                for short in ["CPRESS","COPEN","CNORMF","HFL"]:
                    matches=[n for n in frame.fieldOutputs.keys() if n.split()[0]==short and "S_LOWER/S_UPPER" in n]
                    cfields[short]=values(frame,matches[0],True) if matches else {}
                for node in sorted(nflds["U"]):
                    cpress=cfields["CPRESS"].get(node,0);copen=cfields["COPEN"].get(node,0);cnorm=cfields["CNORMF"].get(node,(0,0));hfl=cfields["HFL"].get(node,0)
                    nodes.writerow([time,node,nflds["NT11"].get(node,float("nan"))]+list(nflds["U"][node][:2])+list(nflds["RF"][node][:2])+[nflds["RFL11"].get(node,float("nan")),cpress,copen,cnorm[1],hfl])
                fld={n:values(frame,n) for n in ["S","EE","PE","CE","PEEQ","CEEQ"]}
                for key in sorted(fld["S"]):
                    row=[time,key[0],key[1]]
                    for n in ["S","EE","PE","CE"]:
                        tensor=list(fld[n].get(key,(0,0,0,0)))
                        if n!="S": tensor[3]*=.5
                        row+=tensor
                    points.writerow(row+[fld["PEEQ"].get(key,0),fld["CEEQ"].get(key,0)])
            elapsed+=step.timePeriod
finally: odb.close()
