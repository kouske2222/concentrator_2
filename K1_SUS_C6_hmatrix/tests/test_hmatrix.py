import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from hmatrix import Kernel,HMatrix,Settings,aca,certify
from run import ROOT,worker

class ACATests(unittest.TestCase):
    def test_complex_low_rank_and_exhaustive_certificate(self):
        rng=np.random.default_rng(57)
        u=rng.normal(size=(90,4))+1j*rng.normal(size=(90,4))
        v=rng.normal(size=(4,75))+1j*rng.normal(size=(4,75))
        a=u@v
        factors=aca(*a.shape,lambda i:a[i],lambda j:a[:,j],1e-9,12)
        self.assertIsNotNone(factors)
        ok,norm,error=certify(*a.shape,lambda s,t:a[s:t],*factors,1e-9)
        self.assertTrue(ok)
        self.assertLess(np.linalg.norm(a-factors[0]@factors[1])/np.linalg.norm(a),1e-10)
        # A sparse defect outside ACA pivot rows cannot evade certification.
        a[33,47]+=10
        ok,_,_=certify(*a.shape,lambda s,t:a[s:t],*factors,1e-9)
        self.assertFalse(ok)

    def test_zero_samples_do_not_authorize_zero_block(self):
        a=np.zeros((90,80),complex); a[33,47]=1
        result=aca(*a.shape,lambda i:a[i],lambda j:a[:,j],1e-5,4)
        if result is not None:
            self.assertTrue(certify(*a.shape,lambda s,t:a[s:t],*result,1e-5)[0])

class HMTests(unittest.TestCase):
    def test_six_modes_direct_cache_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            k=Kernel(ROOT/'config_smoke.nml')
            h=HMatrix(k,Settings(leaf=2),Path(tmp)/'cache','smoke')
            rng=np.random.default_rng(81)
            j=rng.normal(size=18*k.q)+1j*rng.normal(size=18*k.q)
            direct=worker(ROOT/'config_smoke.nml',ROOT/'build/worker','step',j)
            actual=h.apply(j)
            np.testing.assert_allclose(actual,direct['j'],rtol=1e-10,atol=1e-10)
            again=HMatrix(k,Settings(leaf=2),Path(tmp)/'cache','smoke')
            self.assertTrue(again.stats['cache_loaded'])
            np.testing.assert_array_equal(actual,again.apply(j))
            with self.assertRaises(RuntimeError):
                HMatrix(k,Settings(leaf=2),Path(tmp)/'cache','changed')
            manifest=Path(tmp)/'cache/manifest.json'
            original=manifest.read_text()
            bad=json.loads(original); bad['count']-=1
            manifest.write_text(json.dumps(bad))
            with self.assertRaises(RuntimeError):
                HMatrix(k,Settings(leaf=2),Path(tmp)/'cache','smoke')
            manifest.write_text(original)
            path=Path(tmp)/'cache/blocks.npz'
            with path.open('r+b') as f:
                f.seek(20); f.write(b'BROKEN')
            with self.assertRaises(RuntimeError):
                HMatrix(k,Settings(leaf=2),Path(tmp)/'cache','smoke')

    def test_actual_compression_and_fallback(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            tmp=Path(tmp)
            cfg=tmp/'mesh.nml'
            cfg.write_text((ROOT/'config_smoke.nml').read_text().replace(
                'n_face=2, n_z_cone=2, n_z_pipe=2','n_face=6, n_z_cone=10, n_z_pipe=10'))
            k=Kernel(cfg)
            settings=Settings(tol=1e-3,leaf=8,max_rank=24,eta=3,phase_limit=24)
            h=HMatrix(k,settings,tmp/'cache','medium')
            self.assertGreater(h.stats['low_rank_blocks'],0)
            self.assertGreater(h.stats['direct_blocks'],0)
            self.assertLessEqual(h.stats['certified_error2'],settings.tol**2*h.stats['certified_norm2'])
            rng=np.random.default_rng(125)
            j=rng.normal(size=18*k.q)+1j*rng.normal(size=18*k.q)
            expected=worker(cfg,ROOT/'build/worker','step',j)['j']
            rel=np.linalg.norm(h.apply(j)-expected)/np.linalg.norm(expected)
            self.assertLess(rel,3e-3)
            # A rank cap of one must safely fall back to direct blocks.
            fallback=HMatrix(k,Settings(tol=1e-12,leaf=8,max_rank=1,eta=3,phase_limit=24),
                             tmp/'fallback','strict')
            self.assertGreater(fallback.stats['direct_blocks'],0)
            rel=np.linalg.norm(fallback.apply(j)-expected)/np.linalg.norm(expected)
            self.assertLess(rel,1e-10)

    def test_restart_and_setting_change_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            def launch(out,steps=None,tol='1e-6'):
                cmd=[sys.executable,str(ROOT/'run.py'),'--config',str(ROOT/'config_smoke.nml'),
                     '--output',str(tmp/out),'--operator','hmatrix','--hm-tol',tol,'--max-order','3']
                if steps:cmd+=['--steps',str(steps)]
                return subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            for out,steps in [('continuous',None),('resumed',2),('resumed',None)]:
                result=launch(out,steps)
                self.assertEqual(result.returncode,0,result.stdout)
            with np.load(tmp/'continuous/order_000003/state.npz') as a,np.load(tmp/'resumed/order_000003/state.npz') as b:
                for key in ['j','jsum','e','esum']:
                    np.testing.assert_allclose(a[key],b[key],rtol=1e-13,atol=1e-15)
            result=launch('resumed',tol='1e-5')
            self.assertNotEqual(result.returncode,0)
            self.assertIn('changed',result.stdout.lower())
