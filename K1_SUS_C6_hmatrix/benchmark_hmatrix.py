"""Small C6-only comparison. No production mesh is used."""
import json,os,tempfile,time
from pathlib import Path
import numpy as np
from hmatrix import Kernel,HMatrix,Settings
from run import ROOT,worker

records=[]
with tempfile.TemporaryDirectory(prefix='k1-hm-benchmark-') as temp:
    temp=Path(temp)
    cfg=temp/'small.nml'
    cfg.write_text((ROOT/'config_smoke.nml').read_text().replace(
        'n_face=2, n_z_cone=2, n_z_pipe=2','n_face=8, n_z_cone=16, n_z_pipe=16'))
    initial=worker(cfg,ROOT/'build/worker','init')
    direct=[initial]
    direct_times=[]
    for order in range(1,4):
        start=time.perf_counter()
        direct.append(worker(cfg,ROOT/'build/worker','step',direct[-1]['j']))
        direct_times.append(time.perf_counter()-start)
    for tol in [1e-3,1e-5]:
        kernel=Kernel(cfg)
        settings=Settings(tol=tol,leaf=8,max_rank=32,eta=3,phase_limit=24)
        start=time.perf_counter()
        h=HMatrix(kernel,settings,temp/f'cache-{tol}',f'benchmark-{tol}')
        setup=time.perf_counter()-start
        errors=[]; step_times=[]; apply_times=[]
        current=initial
        sum_fast=initial['j'].copy(); sum_ref=initial['j'].copy()
        sum_ef=initial['ei']+initial['e']; sum_er=sum_ef.copy()
        for order in range(1,4):
            start=time.perf_counter()
            next_j=h.apply(current['j'])
            apply_times.append(time.perf_counter()-start)
            current=worker(cfg,ROOT/'build/worker','fields',next_j)
            step_times.append(time.perf_counter()-start)
            sum_fast+=current['j']; sum_ref+=direct[order]['j']
            sum_ef+=current['e']; sum_er+=direct[order]['e']
            rel=lambda a,b:float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
            errors.append(dict(order=order,j=rel(current['j'],direct[order]['j']),
                e=rel(current['e'],direct[order]['e']),jsum=rel(sum_fast,sum_ref),esum=rel(sum_ef,sum_er)))
        stats=dict(h.stats)
        # Verify that cache reuse performs no kernel construction calls.
        before=kernel.entries_evaluated
        reloaded=HMatrix(kernel,settings,temp/f'cache-{tol}',f'benchmark-{tol}')
        assert kernel.entries_evaluated==before
        np.testing.assert_array_equal(h.apply(initial['j']),reloaded.apply(initial['j']))
        records.append(dict(Q=kernel.q,threads=os.environ.get('OMP_NUM_THREADS'),
            settings=settings.__dict__,setup_wall_seconds=setup,stats=stats,
            direct_step_seconds=direct_times,hm_apply_seconds=apply_times,
            hm_step_including_fields_seconds=step_times,errors=errors))
(ROOT/'benchmark_results.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2),flush=True)
