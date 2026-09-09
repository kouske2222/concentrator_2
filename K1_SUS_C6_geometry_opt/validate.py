#!/usr/bin/env python3
"""Run AFTER building on WSL; tiny mesh verifies algorithms, not mesh accuracy."""
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from run import ROOT, worker


def main():
    exe = ROOT / 'build/worker'
    cfg = ROOT / 'config_smoke.nml'
    initial = worker(cfg, exe, 'init')
    modal = worker(cfg, exe, 'step', initial['j'])
    full = worker(cfg, exe, 'full', initial['j'])
    for key in ['j', 'e']:
        err = np.linalg.norm(modal[key]-full[key])/max(np.linalg.norm(full[key]), 1e-300)
        print(f'C6 vs full {key} relative error: {err:.3e}')
        assert err < 1e-10, (key, err)
    with tempfile.TemporaryDirectory(prefix='k1-validation-') as tmp:
        root = Path(tmp)
        def run(out, steps=None):
            cmd = [sys.executable, str(ROOT/'run.py'), '--config', str(cfg),
                   '--output', str(out), '--max-order', '3']
            if steps:
                cmd += ['--steps', str(steps)]
            subprocess.run(cmd, check=True)
        run(root/'continuous')
        run(root/'restart', 2)
        run(root/'restart')
        with np.load(root/'continuous/order_000003/state.npz') as a, np.load(root/'restart/order_000003/state.npz') as b:
            for key in ['j', 'jsum', 'e', 'esum']:
                np.testing.assert_allclose(a[key], b[key], rtol=1e-13, atol=1e-15)
        print('Uninterrupted vs restart: PASS')
        # PEC limit of the material path (large sigma), using a different temporary config.
        pec = root/'pec.nml'
        pec.write_text(cfg.read_text().replace('pec=.false.', 'pec=.true.'))
        limit = root/'limit.nml'
        limit.write_text(cfg.read_text().replace('conductivity=1.37e6', 'conductivity=1.0e30'))
        a, b = worker(pec, exe, 'init'), worker(limit, exe, 'init')
        np.testing.assert_allclose(a['j'], b['j'], rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(a['e'], b['e'], rtol=1e-10, atol=1e-14)
        print('PEC limit (J and radiated E): PASS')
    print('Software checks PASS. Physical mesh/power validation remains required.')


if __name__ == '__main__':
    main()
