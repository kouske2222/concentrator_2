from pathlib import Path
import os
import subprocess
import sys
import tempfile
import time
import json
import numpy as np
from run import worker, ROOT
records=[]
with tempfile.TemporaryDirectory(prefix='k1-geometry-check-') as temp:
    temp=Path(temp)
    for nf,nc,npip in [(2,2,2),(8,12,12)]:
        cfg=temp/'config.nml'
        cfg.write_text((ROOT/'config_smoke.nml').read_text().replace(
            'n_face=2, n_z_cone=2, n_z_pipe=2',
            f'n_face={nf}, n_z_cone={nc}, n_z_pipe={npip}'))
        initial=worker(cfg,ROOT/'build/worker','init')
        states={}
        times={}
        for name in ['worker_old','worker']:
            samples=[]
            for rep in range(3):
                start=time.perf_counter()
                states[name]=worker(cfg,ROOT/'build'/name,'step',initial['j'])
                samples.append(time.perf_counter()-start)
            times[name]=float(np.median(samples))
        errors={key:float(np.linalg.norm(states['worker'][key]-states['worker_old'][key])/
                max(np.linalg.norm(states['worker_old'][key]),1e-300)) for key in ['j','e']}
        records.append(dict(mesh=[nf,nc,npip],Q=2*nf*(nc+npip),median_worker_seconds=times,
                            old_vs_exact_relative_error=errors))
    for out,steps in [('continuous',None),('restart',2),('restart',None)]:
        cmd=[sys.executable,str(ROOT/'run.py'),'--config',str(ROOT/'config_smoke.nml'),
             '--output',str(temp/out),'--max-order','3']
        if steps: cmd+=['--steps',str(steps)]
        subprocess.run(cmd,check=True)
    with np.load(temp/'continuous/order_000003/state.npz') as a, np.load(temp/'restart/order_000003/state.npz') as b:
        for key in ['j','jsum','e','esum']:
            np.testing.assert_allclose(a[key],b[key],rtol=1e-13,atol=1e-15)
    print('C6 uninterrupted/restart PASS')
(ROOT/'geometry_benchmark.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
