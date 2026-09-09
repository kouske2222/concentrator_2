import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from run import save_order, load_latest, converged, metrics, main, ROOT


class RestartTests(unittest.TestCase):
    def state(self):
        return dict(j=np.ones(18, complex), jsum=np.ones(18, complex)*2,
                    e=np.ones(9, complex), esum=np.ones(9, complex)*2,
                    ei=np.ones(9, complex), areas=np.ones(6), xyz=np.zeros(9), jref=np.array(1.))

    def test_roundtrip_and_pending_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            save_order(out, 0, 'abc', self.state(), [{}])
            (out / '.pending-crash').mkdir()
            n, state, hist = load_latest(out, 'abc')
            self.assertEqual(n, 0)
            np.testing.assert_array_equal(state['jsum'], self.state()['jsum'])
            with self.assertRaises(RuntimeError):
                save_order(out, 0, 'abc', self.state(), [{}])

    def test_corrupt_latest_falls_back_and_preserves(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            save_order(out, 0, 'abc', self.state(), [{}])
            save_order(out, 1, 'abc', self.state(), [{}, {}])
            with (out / 'order_000001/state.npz').open('r+b') as f:
                f.truncate(31)
            self.assertEqual(load_latest(out, 'abc')[0], 0)
            self.assertTrue(list(out.glob('order_000001.corrupt-*')))
            save_order(out, 1, 'abc', self.state(), [{}, {}])

    def test_changed_conditions_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            save_order(out, 0, 'abc', self.state(), [{}])
            with self.assertRaises(RuntimeError):
                load_latest(out, 'different')

    def test_interrupted_save_keeps_previous(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            save_order(out, 0, 'abc', self.state(), [{}])
            with patch('run.np.savez', side_effect=OSError('simulated disk full')):
                with self.assertRaises(OSError):
                    save_order(out, 1, 'abc', self.state(), [{}, {}])
            self.assertEqual(load_latest(out, 'abc')[0], 0)
            self.assertFalse((out/'order_000001').exists())

    def test_controller_resume_with_stub_worker(self):
        # Software state-machine test: this stub is NOT an EM validation.
        def stub(config, executable, action, j=None):
            data = self.state()
            data['j'] = np.ones(18, complex) if j is None else j * (0.5 + 0.1j)
            data['e'] = np.ones(9, complex) * data['j'][0]
            return {k: data[k] for k in ['j', 'e', 'ei', 'areas', 'xyz']}
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            executable = root/'stub'
            executable.write_text('not an executable; calls mocked')
            def launch(output, steps=None):
                args = ['run.py', '--config', str(ROOT/'config_smoke.nml'), '--worker', str(executable),
                        '--output', str(root/output), '--max-order', '3']
                if steps:
                    args += ['--steps', str(steps)]
                with patch('sys.argv', args), patch('run.worker', side_effect=stub):
                    main()
            launch('continuous')
            launch('restart', 2)
            launch('restart')
            with np.load(root/'continuous/order_000003/state.npz') as a, np.load(root/'restart/order_000003/state.npz') as b:
                for key in ['j', 'jsum', 'e', 'esum']:
                    np.testing.assert_array_equal(a[key], b[key])

    def test_consecutive_and_minimum(self):
        good = dict(j_ratio=1e-5, j_initial=1e-5, e_ratio=1e-5)
        bad = dict(j_ratio=1., j_initial=1., e_ratio=1.)
        self.assertFalse(converged([good]*10, 1e-4, 1e-3, 5, 10))
        self.assertFalse(converged([good]*10+[bad], 1e-4, 1e-3, 5, 10))
        self.assertTrue(converged([good]*11, 1e-4, 1e-3, 5, 10))

    def test_area_weight_and_complex_cancellation(self):
        j = np.array([1j, 0, 0, 2, 0, 0])
        m = metrics(j, j*10, np.ones(3), np.ones(3)*100, np.array([4.,1.]), 10., np.ones(3))
        self.assertAlmostEqual(m['j_ratio'], .1)
        self.assertAlmostEqual(m['e_ratio'], .01)


if __name__ == '__main__':
    unittest.main()
