#!/usr/bin/env python3
"""Crash-safe, bounded-memory controller for the one-order Fortran worker."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
SCHEMA = 1


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def sync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_json(path, value):
    with Path(path).open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def save_order(out, order, fingerprint, arrays, history):
    """The directory rename is the commit point; partial writes are ignored."""
    target = out / f'order_{order:06d}'
    if target.exists():
        raise RuntimeError(f'Refusing to overwrite {target}')
    stage = Path(tempfile.mkdtemp(prefix='.pending-', dir=out))
    with (stage / 'state.npz').open('wb') as f:
        np.savez(f, **arrays)  # no compression CPU cost on the critical path
        f.flush()
        os.fsync(f.fileno())
    write_json(stage / 'meta.json', dict(schema=SCHEMA, order=order,
               fingerprint=fingerprint, sha256=digest(stage / 'state.npz'), history=history))
    sync_dir(stage)
    os.rename(stage, target)
    sync_dir(out)


def load_latest(out, fingerprint):
    for path in sorted(out.glob('order_[0-9][0-9][0-9][0-9][0-9][0-9]'), reverse=True):
        try:
            meta = json.loads((path / 'meta.json').read_text())
            if meta['schema'] != SCHEMA or meta['fingerprint'] != fingerprint:
                raise RuntimeError('Checkpoint conditions/code differ; use another output directory.')
            if meta['sha256'] != digest(path / 'state.npz'):
                raise ValueError('checksum mismatch')
            with np.load(path / 'state.npz', allow_pickle=False) as data:
                state = {k: data[k] for k in data.files}
            required = {'j', 'jsum', 'e', 'esum', 'ei', 'areas', 'xyz', 'jref'}
            if not required.issubset(state) or any(not np.all(np.isfinite(v)) for v in state.values()):
                raise ValueError('invalid state arrays')
            if state['j'].shape != state['jsum'].shape or state['e'].shape != state['esum'].shape:
                raise ValueError('invalid state shapes')
            if meta['order'] != int(path.name[6:]) or len(meta['history']) != meta['order'] + 1:
                raise ValueError('invalid order history')
            return meta['order'], state, meta['history']
        except RuntimeError:
            raise
        except (OSError, ValueError, KeyError, EOFError) as exc:
            # Preserve damaged files for inspection; recompute from a prior valid order.
            backup = path.with_name(path.name + f'.corrupt-{time.time_ns()}')
            path.rename(backup)
            sync_dir(out)
            print(f'WARNING: {path.name}: {exc}; preserved as {backup.name}', flush=True)
    return -1, None, []


def current_norm(j, areas):
    return float(np.sqrt(np.sum(np.sum(abs(j.reshape(-1, 3))**2, axis=1) * areas)))


def metrics(j, jsum, e, esum, areas, jref, ei):
    nj = current_norm(j, areas)
    ne = float(np.linalg.norm(e))
    return dict(j_ratio=nj / max(current_norm(jsum, areas), jref * 1e-12, 1e-300),
                j_initial=nj / max(jref, 1e-300),
                e_ratio=ne / max(float(np.linalg.norm(esum)), float(np.linalg.norm(ei)) * 1e-12, 1e-300))


def converged(history, tol_j, tol_e, consecutive, min_order):
    if len(history) - 1 < min_order or len(history) < consecutive + 1:
        return False
    return all(h['j_ratio'] < tol_j and h['j_initial'] < tol_j and h['e_ratio'] < tol_e
               for h in history[-consecutive:])


def worker(config, executable, action, j=None):
    # Temporary raw interchange only: never used as a restart checkpoint.
    with tempfile.TemporaryDirectory(prefix='k1-worker-') as tmp:
        tmp = Path(tmp)
        if j is not None:
            np.asarray(j, dtype=np.complex128).tofile(tmp / 'input.bin')
        subprocess.run([str(executable), str(config), action], cwd=tmp, check=True)
        result = {k: np.fromfile(tmp / name, dtype=dtype) for k, name, dtype in [
            ('j', 'current.bin', np.complex128), ('e', 'electric.bin', np.complex128),
            ('ei', 'incident.bin', np.complex128), ('areas', 'areas.bin', np.float64),
            ('xyz', 'probes.bin', np.float64)]}
        if (result['j'].size != 3 * result['areas'].size or result['j'].size == 0
                or result['e'].size != 567 or result['ei'].size != 567 or result['xyz'].size != 567
                or np.any(result['areas'] <= 0)
                or any(not np.all(np.isfinite(a)) for a in result.values())):
            raise RuntimeError('Invalid/nonfinite worker output; previous checkpoint preserved.')
        return result


def preflight(config, max_order):
    text = re.sub(r'!.*', '', config.read_text())
    def value(name):
        match = re.search(r'\b' + name + r'\s*=\s*([\d.eEdD+\-]+)', text, re.I)
        if not match:
            raise ValueError(f'Preflight requires explicit {name}')
        return float(match[1].lower().replace('d', 'e'))
    q = 2 * int(value('n_face')) * (int(value('n_z_cone')) + int(value('n_z_pipe')))
    n = 6 * q
    print(f'170 GHz / C6: Q={q:,} panels/sector, N={n:,} total')
    q_cone = 2 * int(value('n_face')) * int(value('n_z_cone'))
    q_pipe = q - q_cone
    print(f'Candidate pairs/order: {6*q*q:,} (matrix-free, O(Q^2) time)')
    print(f'Pair-map calls after coplanar exclusion: {6*q*q-q_cone*q_cone-q_pipe*q_pipe:,}')
    print(f'Current + cumulative checkpoint/order: about {2*3*n*16/2**30:.3f} GiB')
    print(f'Up to order {max_order}: about {(max_order+1)*2*3*n*16/2**30:.1f} GiB, plus small metadata')
    print('All 6 Fourier modes retained. RAM does not scale with reflection order.')
    print('Mesh/physical accuracy and runtime are NOT certified by this estimate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config_170ghz.nml')
    parser.add_argument('--output', type=Path, default=ROOT / 'results/170ghz')
    parser.add_argument('--worker', type=Path, default=ROOT / 'build/worker')
    parser.add_argument('--max-order', type=int, default=100)
    parser.add_argument('--steps', type=int, help='Stop after this many newly saved orders (including order 0)')
    parser.add_argument('--tol-j', type=float, default=1e-4)
    parser.add_argument('--tol-e', type=float, default=1e-3)
    parser.add_argument('--consecutive', type=int, default=5)
    parser.add_argument('--min-order', type=int, default=10)
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--operator', choices=['direct','hmatrix'], default='direct')
    parser.add_argument('--hm-tol', type=float, default=1e-6)
    parser.add_argument('--hm-leaf', type=int, default=16)
    parser.add_argument('--hm-max-rank', type=int, default=32)
    parser.add_argument('--hm-eta', type=float, default=2.0)
    parser.add_argument('--hm-phase-limit', type=float, default=8.0)
    parser.add_argument('--hm-memory-mib', type=int, default=512)
    args = parser.parse_args()
    if (args.max_order < 0 or args.max_order > 999999 or args.min_order < 1 or args.consecutive < 2
            or not 0 < args.tol_j < 1 or not 0 < args.tol_e < 1 or (args.steps is not None and args.steps < 1)):
        parser.error('Invalid order limits, tolerances, or steps')
    config, executable, out = args.config.resolve(), args.worker.resolve(), args.output.resolve()
    from hmatrix import Settings, Kernel, HMatrix
    from dataclasses import asdict
    hm_settings = Settings(args.hm_tol,args.hm_leaf,args.hm_max_rank,args.hm_eta,
                           args.hm_phase_limit,args.hm_memory_mib)
    try:
        hm_settings.validate()
    except ValueError as exc:
        parser.error(str(exc))
    preflight(config, args.max_order)
    if args.operator == 'hmatrix':
        print('H-matrix: near direct; far ACA with full streamed residual certification.')
        print('First cache build may be O(Q^2); compressed apply cost/rank are not known before building.')
        print('H-matrix settings: '+json.dumps(asdict(hm_settings)))
    if args.preflight:
        return
    if not executable.is_file():
        parser.error('Build the Fortran worker first: make')
    # Byte-exact configuration + executable + controller/source provenance.
    identity = dict(schema=SCHEMA, config_sha256=digest(config), executable_sha256=digest(executable),
                    controller_sha256=digest(Path(__file__)),
                    sources={str(p.relative_to(ROOT)): digest(p) for p in sorted(ROOT.rglob('*.f90'))},
                    tol_j=args.tol_j, tol_e=args.tol_e, consecutive=args.consecutive, min_order=args.min_order)
    identity['operator'] = args.operator
    identity['hm_python_sha256'] = digest(ROOT/'hmatrix.py')
    if args.operator == 'hmatrix':
        if not (ROOT/'build/libhm.so').is_file():
            parser.error('Build the H-matrix library first: make')
        identity['hm_library_sha256'] = digest(ROOT/'build/libhm.so')
        identity['hm_settings'] = asdict(hm_settings)
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True, exist_ok=True)
    with (out / '.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another process owns this output directory')
        manifest = out / 'run.json'
        if manifest.exists():
            if json.loads(manifest.read_text()) != identity:
                raise RuntimeError('Configuration, code, or convergence settings changed. Use a NEW output directory.')
        else:
            temp = out / '.run.json.tmp'
            write_json(temp, identity); os.replace(temp, manifest); sync_dir(out)
        order, state, history = load_latest(out, fingerprint)
        saved = 0
        accelerated = None
        print(f'Resume: last committed order = {order}', flush=True)
        while order < args.max_order:
            if (out / 'STOP').exists():
                print('STOP marker found; paused, not converged.', flush=True)
                return
            if converged(history, args.tol_j, args.tol_e, args.consecutive, args.min_order):
                break
            start = time.monotonic()
            if state is not None and args.operator == 'hmatrix':
                if accelerated is None:
                    accelerated = HMatrix(Kernel(config),hm_settings,out/'hm_cache',fingerprint)
                next_j = accelerated.apply(state['j'])
                new = worker(config,executable,'fields',next_j)
            else:
                new = worker(config, executable, 'init' if state is None else 'step', None if state is None else state['j'])
            if state is None:
                new['jsum'] = new['j'].copy()
                new['esum'] = new['ei'] + new['e']
                new['jref'] = np.array(current_norm(new['j'], new['areas']))
                if float(new['jref']) == 0:
                    raise RuntimeError('No illuminated wall current; check beam and geometry.')
            else:
                new['jsum'] = state['jsum'] + new['j']
                new['esum'] = state['esum'] + new['e']
                new['jref'] = state['jref']
            if any(not np.all(np.isfinite(v)) for v in new.values()):
                raise RuntimeError('Nonfinite accumulated state; last checkpoint remains valid.')
            order += 1
            h = metrics(new['j'], new['jsum'], new['e'], new['esum'], new['areas'], float(new['jref']), new['ei'])
            h.update(order=order, wall_seconds=time.monotonic() - start)
            history = history + [h]
            save_order(out, order, fingerprint, new, history)
            state = new; saved += 1
            print(json.dumps(h), flush=True)
            if args.steps is not None and saved >= args.steps:
                print('Requested step limit: paused; rerun the same command to resume.', flush=True)
                return
        status = 'CONVERGED (order truncation at monitored probes only)' if converged(
            history, args.tol_j, args.tol_e, args.consecutive, args.min_order) else 'MAX_ORDER: NOT CONVERGED'
        print(f'{status}; last committed order={order}', flush=True)


if __name__ == '__main__':
    def stop_signal(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop_signal)
    try:
        main()
    except KeyboardInterrupt:
        print('\nInterrupted. Last committed order is retained; current unfinished order will be repeated.')
        raise SystemExit(130)
