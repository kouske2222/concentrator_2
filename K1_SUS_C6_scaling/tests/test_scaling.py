import ctypes as C
import contextlib
import io
import tempfile
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scaling import ROOT,BASE,hm,break_even,MESHES,configure
sys.path.insert(0,str(BASE))
from run import worker

class ScalingTests(unittest.TestCase):
    def test_break_even_boundary(self):
        self.assertEqual(break_even(10,18,17)['real'],10)
        self.assertEqual(break_even(10,18,17)['strict_profit_order'],11)
        self.assertEqual(break_even(10,18,17)['non_loss_order'],10)
        self.assertEqual(break_even(10,18,3)['strict_profit_order'],1)
        self.assertIsNone(break_even(10,1,1)['real'])
        self.assertIsNone(break_even(10,1,2)['real'])
    def test_mesh_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            for q in MESHES:
                nf,nc,npip=configure(Path(tmp)/'cfg.nml',q)
                self.assertEqual(2*nf*(nc+npip),q)
    def test_bridge_matches_reference_worker(self):
        with contextlib.redirect_stdout(io.StringIO()):
            cfg=BASE/'config_smoke.nml'
            k=hm.Kernel(cfg)
            k.lib.bm_initial.argtypes=[C.c_void_p,C.c_void_p]
            k.lib.bm_fields.argtypes=[C.c_void_p]*3
            j=np.empty(18*k.q,complex);a=np.empty(6*k.q)
            k.lib.bm_initial(j.ctypes.data,a.ctypes.data)
            old=worker(cfg,BASE/'build/worker','init')
            np.testing.assert_array_equal(j,old['j'])
            np.testing.assert_array_equal(a,old['areas'])
            j1=k.apply_direct(j)
            ref=worker(cfg,BASE/'build/worker','step',j)
            np.testing.assert_array_equal(j1,ref['j'])
            e=np.empty(567,complex);h=np.empty(567,complex)
            k.lib.bm_fields(j1.ctypes.data,e.ctypes.data,h.ctypes.data)
            np.testing.assert_array_equal(e,ref['e'])
