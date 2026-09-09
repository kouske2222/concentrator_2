"""Isolated, staged C6 H-matrix scaling measurements. Never runs production Q."""
import argparse
import csv
import ctypes as C
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'K1_SUS_C6_hmatrix'
sys.path.insert(0,str(BASE))
import hmatrix as hm
# Only this measurement process is redirected to its independent library.
hm.ROOT=ROOT
MESHES={512:(8,16,16),1024:(8,32,32),2048:(16,32,32),
        4096:(16,64,64),8192:(32,64,64),16384:(32,128,128)}

def atomic(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def break_even(build,direct,fast):
    saving=direct-fast
    if saving<=0:
        return dict(real=None,non_loss_order=None,strict_profit_order=None,
                    reason='H operator is not faster')
    ratio=build/saving
    return dict(real=ratio,non_loss_order=math.ceil(ratio),
                strict_profit_order=math.floor(ratio)+1,reason=None)

def stats_empty(q,qcone):
    return dict(terminal_blocks=0,low_rank_blocks=0,direct_blocks=0,
        near_direct_blocks=0,admissible_direct_fallback_blocks=0,
        rank_sum=0,rank_max=0,factor_complex_entries=0,
        covered_direct_pairs=0,covered_low_rank_pairs=0,
        candidate_pairs=6*q*q,eligible_pairs=6*q*q-qcone*qcone-(q-qcone)**2,
        rank_histogram={key:0 for key in ['1-4','5-8','9-16','17-32','33+']})

class MeasuredH(hm.HMatrix):
    def __init__(self,kernel,settings,cache,identity,report,save):
        self.metrics=stats_empty(kernel.q,kernel.qcone)
        self.report,self.save_report=report,save
        self.csvfile=(Path(cache).parent/'blocks.csv').open('w',newline='')
        self.writer=csv.DictWriter(self.csvfile,fieldnames=[
            'd','kind','target_panels','source_panels','m','n','rank',
            'D_target','D_source','kD','box_distance','center_distance','R_over_D'])
        self.writer.writeheader()
        self.last=time.monotonic()
        self.save_seconds=0
        try:
            super().__init__(kernel,settings,cache,identity)
        finally:
            self.csvfile.flush();self.csvfile.close()
            report['block_stats']=self.metrics
            save()
    def _save(self):
        t=time.perf_counter()
        super()._save()
        self.save_seconds=time.perf_counter()-t
    def _append(self,d,a,b,u=None,v=None,norm2=0,error2=0):
        super()._append(d,a,b,u,v,norm2,error2)
        D=max(a.diameter,b.diameter)
        gap=np.maximum(0,np.maximum(a.lo[0]-b.hi[d],b.lo[d]-a.hi[0]))
        distance=float(np.linalg.norm(gap))
        center=float(np.linalg.norm((a.lo[0]+a.hi[0]-b.lo[d]-b.hi[d])/2))
        admissible=(distance>0 and D<=self.settings.eta*distance and
                    self.kernel.k*a.diameter*b.diameter<=self.settings.phase_limit*distance)
        kind='low_rank' if u is not None else ('admissible_direct_fallback' if admissible else 'near_direct')
        rank=0 if u is None else u.shape[1]
        m,n=3*len(a.ids),3*len(b.ids)
        self.writer.writerow(dict(d=d,kind=kind,target_panels=len(a.ids),source_panels=len(b.ids),
            m=m,n=n,rank=rank,D_target=a.diameter,D_source=b.diameter,kD=self.kernel.k*D,
            box_distance=distance,center_distance=center,R_over_D=center/D if D else 0))
        s=self.metrics
        s['terminal_blocks']+=1
        if u is None:
            s['direct_blocks']+=1
            s['near_direct_blocks' if not admissible else 'admissible_direct_fallback_blocks']+=1
            s['covered_direct_pairs']+=len(a.ids)*len(b.ids)
        else:
            s['low_rank_blocks']+=1
            s['rank_sum']+=rank;s['rank_max']=max(s['rank_max'],rank)
            s['factor_complex_entries']+=rank*(m+n)
            s['covered_low_rank_pairs']+=len(a.ids)*len(b.ids)
            label='1-4' if rank<=4 else '5-8' if rank<=8 else '9-16' if rank<=16 else '17-32' if rank<=32 else '33+'
            s['rank_histogram'][label]+=1
        if time.monotonic()-self.last>5:
            self.csvfile.flush()
            self.report['block_stats']=s
            self.report['cache_estimated_mib']=self.bytes/2**20
            self.save_report()
            self.last=time.monotonic()

def configure(path,q):
    nf,nc,npip=MESHES[q]
    text=(BASE/'config_smoke.nml').read_text()
    text=text.replace('n_face=2, n_z_cone=2, n_z_pipe=2',
                      f'n_face={nf}, n_z_cone={nc}, n_z_pipe={npip}')
    path.write_text(text)
    return [nf,nc,npip]

def run_child(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    report=dict(Q=args.q,status='running',stage='prepare',stage_started=time.time(),
                settings=asdict(hm.Settings(args.tol,args.leaf,args.rank,args.eta,args.phase,args.cache_mib)),
                threads=os.environ.get('OMP_NUM_THREADS'),samples={},block_stats_complete=False)
    def save():
        report['peak_rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        atomic(out/'result.json',report)
    def stage(name):
        report['stage']=name;report['stage_started']=time.time();save()
        print('STAGE '+name,flush=True)
    def measure(name,fn):
        stage(name);values=[];value=None
        for _ in range(args.repeats):
            t=time.perf_counter();value=fn();values.append(time.perf_counter()-t)
            report['samples'][name]=values.copy();save()
        report[name+'_seconds']=float(np.median(values));save()
        return value
    save()
    try:
        cfg=out/'config.nml';report['mesh']=configure(cfg,args.q)
        kernel=hm.Kernel(cfg)
        if kernel.q!=args.q:raise ValueError('Q does not match requested mesh')
        kernel.lib.bm_initial.argtypes=[C.c_void_p,C.c_void_p]
        kernel.lib.bm_initial.restype=None
        kernel.lib.bm_fields.argtypes=[C.c_void_p]*3
        kernel.lib.bm_fields.restype=None
        j=np.empty(18*kernel.q,complex);areas=np.empty(6*kernel.q)
        kernel.lib.bm_initial(j.ctypes.data,areas.ctypes.data)
        def fields(current):
            e=np.empty(567,complex);h=np.empty(567,complex)
            kernel.lib.bm_fields(current.ctypes.data,e.ctypes.data,h.ctypes.data)
            return e
        reference=measure('direct_operator',lambda:kernel.apply_direct(j))
        e_ref=measure('field_integral',lambda:fields(reference))
        stage('hm_build')
        identity=hashlib.sha256((cfg.read_text()+json.dumps(report['settings'],sort_keys=True)).encode()).hexdigest()
        t=time.perf_counter()
        h=MeasuredH(kernel,hm.Settings(**report['settings']),out/'cache',identity,report,save)
        report['hm_build_seconds']=time.perf_counter()-t
        report['hm_assembly_seconds']=h.stats.get('build_seconds')
        report['hm_cache_save_seconds']=h.save_seconds
        report['hm_stats']=h.stats
        report['block_stats_complete']=True
        s=report['block_stats']
        if s['covered_direct_pairs']+s['covered_low_rank_pairs']!=s['eligible_pairs']:
            raise AssertionError('Terminal blocks do not cover the eligible operator')
        s['rank_mean']=s['rank_sum']/s['low_rank_blocks'] if s['low_rank_blocks'] else None
        s['direct_block_fraction']=s['direct_blocks']/max(s['terminal_blocks'],1)
        s['direct_pair_fraction']=s['covered_direct_pairs']/s['eligible_pairs']
        s['factor_bytes']=16*s['factor_complex_entries']
        save()
        actual=measure('hm_operator',lambda:h.apply(j))
        stage('accuracy')
        e_actual=fields(actual)
        relative=lambda a,b:float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
        weighted=lambda x:float(np.sqrt(np.sum(np.sum(abs(x.reshape(-1,3))**2,axis=1)*areas)))
        report['errors']=dict(j=relative(actual,reference),e=relative(e_actual,e_ref),
            j_area_weighted=weighted(actual-reference)/max(weighted(reference),1e-300))
        report['break_even']=break_even(report['hm_build_seconds'],report['direct_operator_seconds'],report['hm_operator_seconds'])
        report['operator_speedup']=report['direct_operator_seconds']/report['hm_operator_seconds']
        report['status']='complete'
        save()
    except Exception as exc:
        report['status']='resource_limit' if isinstance(exc,MemoryError) else 'error'
        report['error']=str(exc)
        save();traceback.print_exc()
        return 1
    return 0

def rss_mib(pid):
    try:
        for line in Path(f'/proc/{pid}/status').read_text().splitlines():
            if line.startswith('VmRSS:'):return int(line.split()[1])/1024
    except FileNotFoundError:pass
    return 0

def orchestrate(args):
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock=(out/'.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    settings=dict(q_values=args.q_values,seconds_per_stage=args.seconds_per_stage,rss_mib=args.rss_mib,
        threads=args.threads,repeats=args.repeats,tol=args.tol,leaf=args.leaf,rank=args.rank,
        eta=args.eta,phase=args.phase,cache_mib=args.cache_mib)
    paths=[Path(__file__),ROOT/'timing_bridge.f90',ROOT/'build/libhm.so',BASE/'hmatrix.py',
           *sorted((BASE/'common').glob('*.f90')),*sorted((BASE/'c6_modal').glob('*.f90')),
           BASE/'acceleration/mod_bridge.f90',BASE/'config_smoke.nml']
    manifest=dict(settings=settings,sha256={str(p):hm.sha256(p) for p in paths},
                  cpu=Path('/proc/cpuinfo').read_text().split('model name')[1].splitlines()[0].strip(' :'))
    if (out/'manifest.json').exists():
        if json.loads((out/'manifest.json').read_text())!=manifest:
            raise RuntimeError('Benchmark conditions changed; choose a new output directory')
    else:atomic(out/'manifest.json',manifest)
    for q in args.q_values:
        qout=out/f'Q{q}';qout.mkdir(exist_ok=True)
        result_path=qout/'result.json'
        if result_path.exists() and json.loads(result_path.read_text()).get('status') in ['complete','timeout','resource_limit']:
            print(f'Q={q}: preserved existing result',flush=True);continue
        if (qout/'cache/manifest.json').exists():
            raise RuntimeError('Partial prior run has a completed cache: use a new output to time a fresh build')
        env=os.environ.copy();env.update(OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1',
                                        MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
        cmd=[sys.executable,str(Path(__file__).resolve()),'--child','--q',str(q),'--output',str(qout),
             '--repeats',str(args.repeats),'--tol',str(args.tol),'--leaf',str(args.leaf),
             '--rank',str(args.rank),'--eta',str(args.eta),'--phase',str(args.phase),
             '--cache-mib',str(args.cache_mib)]
        print(f'Q={q}: starting',flush=True)
        t=time.time();peak=0;last_stage=None;last_message=time.monotonic()
        with (qout/'run.log').open('w') as log:
            p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                while p.poll() is None:
                    time.sleep(0.5);peak=max(peak,rss_mib(p.pid))
                    r=json.loads(result_path.read_text()) if result_path.exists() else dict(stage='startup',stage_started=t)
                    if r['stage']!=last_stage:
                        print(f"Q={q}: {r['stage']}",flush=True);last_stage=r['stage']
                    elapsed=time.time()-r.get('stage_started',t)
                    reason=None
                    if elapsed>args.seconds_per_stage:reason='timeout'
                    if peak>args.rss_mib:reason='resource_limit'
                    if reason:
                        os.killpg(p.pid,signal.SIGKILL);p.wait()
                        r.update(status=reason,limit_stage=r['stage'],elapsed_stage_seconds=elapsed,
                                 supervisor_peak_rss_mib=peak)
                        atomic(result_path,r)
                        print(f'Q={q}: {reason} at {last_stage}',flush=True)
                        break
                    if time.monotonic()-last_message>30:
                        print(f'Q={q}: {last_stage}, {elapsed:.0f}s, peak RSS {peak:.0f} MiB',flush=True)
                        last_message=time.monotonic()
            except BaseException:
                if p.poll() is None:os.killpg(p.pid,signal.SIGKILL);p.wait()
                raise
        r=json.loads(result_path.read_text()) if result_path.exists() else dict(Q=q,status='error',error='No child report')
        if r.get('status')=='running':r.update(status='error',error=f'Child exited {p.returncode}')
        r['supervisor_peak_rss_mib']=peak
        atomic(result_path,r)
        print(f"Q={q}: {r['status']}",flush=True)
        subprocess.run([sys.executable,str(ROOT/'report.py'),str(out)],check=True)
    atomic(out/'DONE.json',dict(finished=time.time(),q_values=args.q_values))
    print('Scaling sweep finished',flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--child',action='store_true')
    p.add_argument('--q',type=int,choices=sorted(MESHES),default=512)
    p.add_argument('--q-values',nargs='+',type=int,choices=sorted(MESHES),default=list(MESHES))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seconds-per-stage',type=float,default=1800)
    p.add_argument('--rss-mib',type=float,default=2048)
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--tol',type=float,default=1e-5)
    p.add_argument('--leaf',type=int,default=8)
    p.add_argument('--rank',type=int,default=32)
    p.add_argument('--eta',type=float,default=3)
    p.add_argument('--phase',type=float,default=24)
    p.add_argument('--cache-mib',type=int,default=512)
    args=p.parse_args()
    if min(args.seconds_per_stage,args.rss_mib,args.threads,args.repeats)<=0:p.error('Limits must be positive')
    hm.Settings(args.tol,args.leaf,args.rank,args.eta,args.phase,args.cache_mib).validate()
    if args.child:raise SystemExit(run_child(args))
    orchestrate(args)

if __name__=='__main__':main()
