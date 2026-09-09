#!/usr/bin/env python3
"""Serial timing of committed output-disabled Fuelsim/Abaqus inputs on one host.

Run zero warms each program/case. Runs one and two are reported measurements.
No input text is generated or modified by this runner.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'verification/abaqus/rz_performance'


def run(fuelsim, powershell, windows_source, cpu, resume=False):
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    records=json.loads((REFERENCE/'timing.json').read_text()) if resume else []
    if resume:
        import hashlib
        provenance=json.loads((REFERENCE/'provenance.json').read_text())
        for name,digest in provenance['sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise RuntimeError('Cannot resume changed evidence: '+name)
        if any(r['solver']=='fuelsim' and r['cpu']!=cpu for r in records):
            raise RuntimeError('Cannot resume with a different CPU affinity')
        if hashlib.sha256(fuelsim.read_bytes()).hexdigest()!=provenance['sha256']['build/fuelsim']:
            raise RuntimeError('Cannot resume with a different executable')
    def recorded(solver,size,repeat):
        return any(r['solver']==solver and r['size']==size and r['run']==repeat for r in records)
    for size in ('medium','large'):
        card=ROOT/'verification/fuelsim'/('steady_rz_performance_'+size+'_timing.fsi')
        for repeat in range(3):
            if recorded('fuelsim',size,repeat):continue
            log=REFERENCE/('fuelsim_'+size+'_run%d.log'%repeat)
            print('Starting Fuelsim',size,'run',repeat,flush=True)
            start=time.perf_counter()
            with log.open('w') as output:
                subprocess.run(['taskset','-c',str(cpu),str(fuelsim),'-i',str(card)],env=env,stdout=output,stderr=subprocess.STDOUT,check=True)
            seconds=time.perf_counter()-start
            text=log.read_text()
            for marker in ('completed=true','load_steps_completed=20','mpi_ranks=1'):
                if not re.search('^'+marker+'$',text,re.M):raise RuntimeError('Missing '+marker)
            values=dict(re.findall(r'^([^=\n]+)=([^\n]+)$',text,re.M))
            records.append(dict(solver='fuelsim',size=size,run=repeat,warmup=repeat==0,
                                external_wall_seconds=seconds,solver_internal_seconds=float(values['total_seconds']),
                                nonlinear_iterations=int(values['nonlinear_iterations_total']),cpu=cpu))
            (REFERENCE/'timing.json').write_text(json.dumps(records,indent=2)+'\n')
            print('Completed:',seconds,'seconds',flush=True)
        if all(recorded('abaqus',size,r) for r in range(3)):continue
        print('Starting Abaqus',size,'warmup plus two measurements',flush=True)
        completed_files=all((REFERENCE/('rz_performance_'+size+'_timing_run%d_timing.txt'%r)).exists()
                            and (REFERENCE/('rz_performance_'+size+'_timing_run%d.sta'%r)).exists()
                            and 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' in (REFERENCE/('rz_performance_'+size+'_timing_run%d.sta'%r)).read_text(encoding='latin-1') for r in range(3))
        if not (resume and completed_files):
            with (REFERENCE/('abaqus_'+size+'_timing_launcher.log')).open('w') as output:
                subprocess.run([powershell,'-NoProfile','-ExecutionPolicy','Bypass','-File',windows_source+r'\run.ps1',
                                '-SourceDirectory',windows_source,'-Size',size,'-Timing','-Runs','3'],stdout=output,stderr=subprocess.STDOUT,check=True)
        for repeat in range(3):
            if recorded('abaqus',size,repeat):continue
            stem=REFERENCE/('rz_performance_'+size+'_timing_run%d'%repeat)
            values=dict(re.findall(r'^([^=\n]+)=([^\n]+)$',Path(str(stem)+'_timing.txt').read_text(),re.M))
            dat=Path(str(stem)+'.dat').read_text(encoding='latin-1');msg=Path(str(stem)+'.msg').read_text(encoding='latin-1')
            iterations=re.findall(r'(\d+)\s+ITERATIONS INCLUDING CONTACT',msg)
            increments=re.findall(r'TOTAL OF\s+(\d+)\s+INCREMENTS',msg)
            if not increments or int(increments[-1])!=20:raise RuntimeError('Abaqus increment count differs')
            wall=re.findall(r'WALLCLOCK TIME \(SEC\)\s*=\s*([\d.]+)',dat)
            records.append(dict(solver='abaqus',size=size,run=repeat,warmup=repeat==0,
                                external_wall_seconds=float(values['external_wall_seconds']),
                                analysis_wall_seconds=float(wall[-1]),nonlinear_iterations=int(iterations[-1]),cpu=0))
        (REFERENCE/'timing.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fuelsim',type=Path,default=ROOT/'build/fuelsim')
    p.add_argument('--powershell',default='/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')
    p.add_argument('--windows-source',required=True);p.add_argument('--cpu',type=int,default=0)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();run(a.fuelsim.resolve(),a.powershell,a.windows_source,a.cpu,a.resume)
