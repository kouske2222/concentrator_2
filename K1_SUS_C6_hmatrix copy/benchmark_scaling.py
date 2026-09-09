#!/usr/bin/env python3
"""
Scaling benchmark for K1_SUS_C6_hmatrix.

Purpose
-------
Phase 1:
    Check whether the H-matrix remains viable as Q increases.

Phase 2:
    Measure operator-only wall time:
        T_direct-op = Kernel.apply_direct(j)
        T_H-op      = HMatrix.apply(j)

    No observation-field calculation, subprocess launch, checkpoint I/O,
    or IPO accumulation is included in these operator timings.

Phase 3:
    Export H-matrix rank statistics and block geometry:
        - rank distribution
        - direct fallback fraction
        - block-count weighting
        - dense-entry weighting
        - algebraic apply-work proxy
        - distance in wavelengths
        - electrical cluster sizes k D_t and k D_s
        - chi = k D_t D_s / R

Phase 4:
    Compute the operator-level break-even application count:
        N_break = T_build / (T_direct-op - T_H-op)

Notes
-----
* The production physics kernel is NOT modified.
* The benchmark uses the same 170 GHz geometry/material configuration and
  scales only n_face, n_z_cone, and n_z_pipe.
* HMatrix performs exhaustive streamed residual certification, so initial
  construction can still scale approximately as O(Q^2).
* The H-matrix cache is deliberately built in a fresh temporary directory
  for each Q so that T_build is always a first-build time.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import os
import re
import tempfile
import time
import traceback
from dataclasses import asdict
from pathlib import Path
from typing import Callable

import numpy as np

from hmatrix import HMatrix, Kernel, Settings


ROOT = Path(__file__).resolve().parent

DEFAULT_TARGET_Q = [
    512,
    1024,
    2048,
    4096,
    8192,
    16384,
]


# ---------------------------------------------------------------------------
# Configuration utilities
# ---------------------------------------------------------------------------

def read_int_parameter(text: str, name: str) -> int:
    match = re.search(
        rf"\b{re.escape(name)}\s*=\s*([+-]?\d+)",
        re.sub(r"!.*", "", text),
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"Could not find integer parameter {name!r} in config")
    return int(match.group(1))


def replace_int_parameter(text: str, name: str, value: int) -> str:
    pattern = re.compile(
        rf"(\b{re.escape(name)}\s*=\s*)([+-]?\d+)",
        flags=re.IGNORECASE,
    )
    updated, count = pattern.subn(rf"\g<1>{int(value)}", text, count=1)
    if count != 1:
        raise ValueError(f"Could not replace parameter {name!r} in config")
    return updated


def production_mesh_parameters(config_path: Path) -> dict:
    text = config_path.read_text()
    n_face = read_int_parameter(text, "n_face")
    n_z_cone = read_int_parameter(text, "n_z_cone")
    n_z_pipe = read_int_parameter(text, "n_z_pipe")
    q = 2 * n_face * (n_z_cone + n_z_pipe)
    return {
        "n_face": n_face,
        "n_z_cone": n_z_cone,
        "n_z_pipe": n_z_pipe,
        "Q": q,
    }


def scaled_mesh_parameters(
    target_q: int,
    production: dict,
) -> dict:
    """
    Scale circumferential and axial resolutions by the same linear factor.

    Since
        Q = 2*n_face*(n_z_cone+n_z_pipe),
    scaling all mesh counts by s gives approximately Q ~ s^2 Q_prod.
    """
    if target_q < 1:
        raise ValueError("target_q must be positive")

    q_prod = production["Q"]
    scale = math.sqrt(target_q / q_prod)

    n_face = max(2, int(round(production["n_face"] * scale)))
    n_z_cone = max(2, int(round(production["n_z_cone"] * scale)))
    n_z_pipe = max(2, int(round(production["n_z_pipe"] * scale)))

    actual_q = 2 * n_face * (n_z_cone + n_z_pipe)

    return {
        "target_Q": int(target_q),
        "actual_Q_requested": int(actual_q),
        "scale": float(scale),
        "n_face": int(n_face),
        "n_z_cone": int(n_z_cone),
        "n_z_pipe": int(n_z_pipe),
    }


def write_scaled_config(
    base_config: Path,
    output_config: Path,
    mesh: dict,
) -> None:
    text = base_config.read_text()
    text = replace_int_parameter(text, "n_face", mesh["n_face"])
    text = replace_int_parameter(text, "n_z_cone", mesh["n_z_cone"])
    text = replace_int_parameter(text, "n_z_pipe", mesh["n_z_pipe"])
    output_config.write_text(text)


# ---------------------------------------------------------------------------
# Benchmark current and timing
# ---------------------------------------------------------------------------

def make_test_current(q: int, seed: int) -> np.ndarray:
    """
    Deterministic arbitrary complex current.

    Shape is the flattened physical C6 current:
        6 sectors * Q panels/sector * 3 local components = 18 Q.
    The arbitrary complex field excites all six Fourier modes.
    """
    rng = np.random.default_rng(seed)
    j = (
        rng.standard_normal(18 * q)
        + 1j * rng.standard_normal(18 * q)
    ).astype(np.complex128)

    norm = np.linalg.norm(j)
    if not np.isfinite(norm) or norm == 0.0:
        raise RuntimeError("Failed to create finite benchmark current")
    j /= norm
    return j


def timed_call(func: Callable[[], np.ndarray]) -> tuple[float, np.ndarray]:
    gc_was_enabled = gc.isenabled()
    try:
        gc.disable()
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
    finally:
        if gc_was_enabled:
            gc.enable()
    return float(elapsed), result


def benchmark_operator_pair(
    kernel: Kernel,
    hmatrix: HMatrix,
    j: np.ndarray,
    repeats: int,
    warmup: int,
) -> dict:
    """
    Compare exactly the two operator paths.

    direct:
        Kernel.apply_direct(j)

    H-matrix:
        HMatrix.apply(j)

    The order of the two measurements is alternated between repetitions
    to reduce systematic thermal/order bias.
    """
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    if warmup < 0:
        raise ValueError("warmup must be >= 0")

    for _ in range(warmup):
        _ = kernel.apply_direct(j)
        _ = hmatrix.apply(j)

    direct_times: list[float] = []
    h_times: list[float] = []
    direct_result = None
    h_result = None

    for rep in range(repeats):
        if rep % 2 == 0:
            td, direct_result = timed_call(lambda: kernel.apply_direct(j))
            th, h_result = timed_call(lambda: hmatrix.apply(j))
        else:
            th, h_result = timed_call(lambda: hmatrix.apply(j))
            td, direct_result = timed_call(lambda: kernel.apply_direct(j))

        direct_times.append(td)
        h_times.append(th)

    assert direct_result is not None
    assert h_result is not None

    denom = max(float(np.linalg.norm(direct_result)), 1.0e-300)
    relative_error = float(
        np.linalg.norm(h_result - direct_result) / denom
    )

    t_direct = float(np.median(direct_times))
    t_h = float(np.median(h_times))

    if t_h > 0.0:
        speedup = t_direct / t_h
    else:
        speedup = math.inf

    return {
        "direct_seconds": {
            "median": t_direct,
            "min": float(np.min(direct_times)),
            "max": float(np.max(direct_times)),
            "all": [float(x) for x in direct_times],
        },
        "hmatrix_seconds": {
            "median": t_h,
            "min": float(np.min(h_times)),
            "max": float(np.max(h_times)),
            "all": [float(x) for x in h_times],
        },
        "speedup_direct_over_hmatrix": float(speedup),
        "relative_operator_error_l2": relative_error,
    }


# ---------------------------------------------------------------------------
# H-matrix block statistics
# ---------------------------------------------------------------------------

def rank_category(rank: int | None) -> str:
    if rank is None:
        return "direct"
    if rank <= 4:
        return "1-4"
    if rank <= 8:
        return "5-8"
    if rank <= 16:
        return "9-16"
    return "17-32+"


def power_of_two_bin(value: float) -> str:
    """
    Human-readable logarithmic bin.

    Examples:
        0              -> "0"
        0.7            -> "[2^-1,2^0)"
        3.2            -> "[2^1,2^2)"
        inf            -> "inf"
    """
    if not np.isfinite(value):
        return "inf"
    if value <= 0.0:
        return "0"
    exponent = int(math.floor(math.log2(value)))
    return f"[2^{exponent},2^{exponent + 1})"


def block_geometry(
    kernel: Kernel,
    d: int,
    rows: np.ndarray,
    cols: np.ndarray,
) -> dict:
    """
    Reproduce the geometry quantities used by HMatrix._build.

    Important:
    Node.diameter is defined using sector 0 geometry for both target and source
    clusters.  The separation distance uses target sector 0 and source sector d.
    """
    target0 = kernel.xyz[0, rows, :]
    source0 = kernel.xyz[0, cols, :]
    source_d = kernel.xyz[d, cols, :]

    t_lo = target0.min(axis=0)
    t_hi = target0.max(axis=0)

    s0_lo = source0.min(axis=0)
    s0_hi = source0.max(axis=0)

    sd_lo = source_d.min(axis=0)
    sd_hi = source_d.max(axis=0)

    diameter_t = float(np.linalg.norm(t_hi - t_lo))
    diameter_s = float(np.linalg.norm(s0_hi - s0_lo))

    gap = np.maximum(
        0.0,
        np.maximum(t_lo - sd_hi, sd_lo - t_hi),
    )
    distance = float(np.linalg.norm(gap))

    kdt = float(kernel.k * diameter_t)
    kds = float(kernel.k * diameter_s)
    max_kd = max(kdt, kds)

    distance_lambda = float(kernel.k * distance / (2.0 * np.pi))

    if distance > 0.0:
        chi = float(
            kernel.k * diameter_t * diameter_s / distance
        )
    else:
        chi = math.inf

    return {
        "diameter_target_m": diameter_t,
        "diameter_source_m": diameter_s,
        "kD_target": kdt,
        "kD_source": kds,
        "max_kD": max_kd,
        "distance_m": distance,
        "distance_lambda": distance_lambda,
        "chi_kDtDs_over_R": chi,
        "distance_bin_lambda": power_of_two_bin(distance_lambda),
        "max_kD_bin": power_of_two_bin(max_kd),
        "chi_bin": power_of_two_bin(chi),
    }


def collect_block_records(
    kernel: Kernel,
    hmatrix: HMatrix,
) -> list[dict]:
    records: list[dict] = []

    for block_index, (d, rows, cols, u, v) in enumerate(hmatrix.blocks):
        nr = 3 * len(rows)
        nc = 3 * len(cols)

        if u is None:
            rank = None
            direct = True
            low_rank_factor_entries = 0
            algebraic_apply_work = nr * nc
        else:
            rank = int(u.shape[1])
            direct = False
            low_rank_factor_entries = rank * (nr + nc)
            algebraic_apply_work = rank * (nr + nc)

        dense_entries = nr * nc
        geometry = block_geometry(kernel, d, rows, cols)

        record = {
            "block_index": int(block_index),
            "d": int(d),
            "target_panels": int(len(rows)),
            "source_panels": int(len(cols)),
            "nrow_dof": int(nr),
            "ncol_dof": int(nc),
            "rank": "" if rank is None else int(rank),
            "rank_category": rank_category(rank),
            "direct_fallback": bool(direct),
            "dense_equivalent_entries": int(dense_entries),
            "low_rank_factor_entries": int(low_rank_factor_entries),
            "algebraic_apply_work_proxy": int(algebraic_apply_work),
        }
        record.update(geometry)
        records.append(record)

    return records


def weighted_rank_summary(
    block_records: list[dict],
) -> dict:
    categories = ["1-4", "5-8", "9-16", "17-32+", "direct"]

    summary = {
        cat: {
            "block_count": 0,
            "dense_equivalent_entries": 0,
            "apply_work_proxy": 0,
        }
        for cat in categories
    }

    for row in block_records:
        cat = row["rank_category"]
        summary[cat]["block_count"] += 1
        summary[cat]["dense_equivalent_entries"] += row[
            "dense_equivalent_entries"
        ]
        summary[cat]["apply_work_proxy"] += row[
            "algebraic_apply_work_proxy"
        ]

    total_blocks = sum(x["block_count"] for x in summary.values())
    total_dense = sum(
        x["dense_equivalent_entries"] for x in summary.values()
    )
    total_work = sum(x["apply_work_proxy"] for x in summary.values())

    for cat in categories:
        item = summary[cat]
        item["block_percent"] = (
            100.0 * item["block_count"] / total_blocks
            if total_blocks else 0.0
        )
        item["dense_entry_percent"] = (
            100.0 * item["dense_equivalent_entries"] / total_dense
            if total_dense else 0.0
        )
        item["apply_work_percent"] = (
            100.0 * item["apply_work_proxy"] / total_work
            if total_work else 0.0
        )

    return {
        "categories": summary,
        "totals": {
            "block_count": int(total_blocks),
            "dense_equivalent_entries": int(total_dense),
            "apply_work_proxy": int(total_work),
        },
    }


def geometry_rank_summary(
    block_records: list[dict],
) -> list[dict]:
    """
    Aggregate by:
        distance in wavelengths (power-of-two bin)
        max(k D_t, k D_s)      (power-of-two bin)

    This gives a compact table for studying how rank/fallback depends on
    separation and electrical cluster size.
    """
    groups: dict[tuple[str, str], list[dict]] = {}

    for row in block_records:
        key = (
            row["distance_bin_lambda"],
            row["max_kD_bin"],
        )
        groups.setdefault(key, []).append(row)

    output: list[dict] = []

    for (distance_bin, kd_bin), rows in groups.items():
        low_ranks = [
            int(row["rank"])
            for row in rows
            if not row["direct_fallback"]
        ]

        direct_count = sum(
            1 for row in rows if row["direct_fallback"]
        )

        dense_total = sum(
            row["dense_equivalent_entries"] for row in rows
        )

        direct_dense = sum(
            row["dense_equivalent_entries"]
            for row in rows
            if row["direct_fallback"]
        )

        output.append({
            "distance_bin_lambda": distance_bin,
            "max_kD_bin": kd_bin,
            "block_count": len(rows),
            "low_rank_block_count": len(rows) - direct_count,
            "direct_block_count": direct_count,
            "direct_block_percent": (
                100.0 * direct_count / len(rows)
                if rows else 0.0
            ),
            "direct_dense_entry_percent": (
                100.0 * direct_dense / dense_total
                if dense_total else 0.0
            ),
            "mean_low_rank": (
                float(np.mean(low_ranks))
                if low_ranks else math.nan
            ),
            "median_low_rank": (
                float(np.median(low_ranks))
                if low_ranks else math.nan
            ),
            "max_low_rank": (
                int(max(low_ranks))
                if low_ranks else ""
            ),
        })

    output.sort(
        key=lambda row: (
            row["distance_bin_lambda"],
            row["max_kD_bin"],
        )
    )
    return output


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        path.write_text("")
        return

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w") as f:
        json.dump(data, f, indent=2, allow_nan=True)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def flatten_summary_record(record: dict) -> dict:
    operator = record.get("operator", {})
    direct = operator.get("direct_seconds", {})
    hm = operator.get("hmatrix_seconds", {})
    hstats = record.get("hmatrix_stats", {})
    break_even = record.get("break_even", {})

    return {
        "target_Q": record.get("target_Q"),
        "actual_Q": record.get("actual_Q"),
        "n_face": record.get("mesh", {}).get("n_face"),
        "n_z_cone": record.get("mesh", {}).get("n_z_cone"),
        "n_z_pipe": record.get("mesh", {}).get("n_z_pipe"),
        "kernel_init_seconds": record.get("kernel_init_seconds"),
        "T_build_total_seconds": record.get(
            "T_build_total_seconds"
        ),
        "T_direct_operator_seconds": direct.get("median"),
        "T_hmatrix_operator_seconds": hm.get("median"),
        "operator_speedup": operator.get(
            "speedup_direct_over_hmatrix"
        ),
        "relative_operator_error_l2": operator.get(
            "relative_operator_error_l2"
        ),
        "N_break": break_even.get("N_break"),
        "first_profitable_application": break_even.get(
            "first_profitable_application"
        ),
        "low_rank_blocks": hstats.get("low_rank_blocks"),
        "direct_blocks": hstats.get("direct_blocks"),
        "zero_blocks": hstats.get("zero_blocks"),
        "stored_bytes": hstats.get("stored_bytes"),
        "scalar_kernel_evaluations": hstats.get(
            "scalar_kernel_evaluations"
        ),
        "status": record.get("status"),
        "error": record.get("error"),
    }


# ---------------------------------------------------------------------------
# Break-even
# ---------------------------------------------------------------------------

def calculate_break_even(
    t_build: float,
    t_direct: float,
    t_h: float,
) -> dict:
    """
    Operator-level comparison:

        direct total = N * T_direct
        H total      = T_build + N * T_H

    H is strictly faster when:
        N > T_build / (T_direct - T_H)
    """
    difference = t_direct - t_h

    if difference <= 0.0:
        return {
            "exists": False,
            "N_break": math.inf,
            "first_profitable_application": None,
            "reason": (
                "T_H >= T_direct, so no positive break-even "
                "application count exists."
            ),
        }

    n_break = t_build / difference
    first_profitable = int(math.floor(n_break)) + 1

    return {
        "exists": True,
        "N_break": float(n_break),
        "first_profitable_application": int(first_profitable),
        "reason": None,
    }


# ---------------------------------------------------------------------------
# One Q benchmark
# ---------------------------------------------------------------------------

def run_one_case(
    target_q: int,
    base_config: Path,
    production: dict,
    settings: Settings,
    output_dir: Path,
    repeats: int,
    warmup: int,
    seed: int,
) -> dict:
    mesh = scaled_mesh_parameters(target_q, production)

    case_dir = output_dir / f"Q_target_{target_q:06d}"
    case_dir.mkdir(parents=True, exist_ok=True)

    cfg = case_dir / "config_scaled.nml"
    write_scaled_config(base_config, cfg, mesh)

    print(
        "\n"
        + "=" * 78
        + f"\nTarget Q={target_q:,}"
        + f"\nmesh: n_face={mesh['n_face']}, "
          f"n_z_cone={mesh['n_z_cone']}, "
          f"n_z_pipe={mesh['n_z_pipe']}"
        + f"\nrequested actual Q={mesh['actual_Q_requested']:,}"
        + "\n"
        + "=" * 78,
        flush=True,
    )

    # Kernel initialization is recorded separately from H-matrix build.
    t0 = time.perf_counter()
    kernel = Kernel(cfg)
    kernel_init_seconds = time.perf_counter() - t0

    actual_q = int(kernel.q)
    if actual_q != mesh["actual_Q_requested"]:
        raise RuntimeError(
            "Q mismatch between Python scaling formula and Fortran mesh: "
            f"{mesh['actual_Q_requested']} vs {actual_q}"
        )

    print(
        f"Kernel initialized: Q={actual_q:,}, "
        f"{kernel_init_seconds:.6f} s",
        flush=True,
    )

    # Fresh cache for a genuine first-build measurement.
    with tempfile.TemporaryDirectory(
        prefix=f"k1-hm-scale-Q{actual_q}-"
    ) as temp_cache:
        cache_dir = Path(temp_cache) / "hm_cache"

        identity = {
            "benchmark": "benchmark_scaling",
            "target_Q": int(target_q),
            "actual_Q": actual_q,
            "mesh": mesh,
            "settings": asdict(settings),
        }

        entries_before = int(kernel.entries_evaluated)

        t0 = time.perf_counter()
        h = HMatrix(
            kernel,
            settings,
            cache_dir,
            identity,
        )
        t_build_total = time.perf_counter() - t0

        build_kernel_evaluations = (
            int(kernel.entries_evaluated) - entries_before
        )

        print(
            f"H build total: {t_build_total:.6f} s",
            flush=True,
        )

        # Phase 3 is intentionally outside T_build and T_operator timings.
        block_records = collect_block_records(kernel, h)
        rank_summary = weighted_rank_summary(block_records)
        geometry_summary = geometry_rank_summary(block_records)

        write_csv(
            case_dir / "block_stats.csv",
            block_records,
        )
        write_csv(
            case_dir / "rank_by_geometry.csv",
            geometry_summary,
        )

        # Phase 2
        j = make_test_current(actual_q, seed)

        operator = benchmark_operator_pair(
            kernel=kernel,
            hmatrix=h,
            j=j,
            repeats=repeats,
            warmup=warmup,
        )

        t_direct = operator["direct_seconds"]["median"]
        t_h = operator["hmatrix_seconds"]["median"]

        # Phase 4
        break_even = calculate_break_even(
            t_build=t_build_total,
            t_direct=t_direct,
            t_h=t_h,
        )

        hstats = dict(h.stats)
        hstats["build_kernel_evaluations_this_run"] = (
            build_kernel_evaluations
        )

        record = {
            "status": "ok",
            "target_Q": int(target_q),
            "actual_Q": actual_q,
            "mesh": mesh,
            "settings": asdict(settings),
            "kernel_init_seconds": float(kernel_init_seconds),
            "T_build_total_seconds": float(t_build_total),
            "hmatrix_stats": hstats,
            "operator": operator,
            "break_even": break_even,
            "rank_summary": rank_summary,
            "files": {
                "scaled_config": str(cfg),
                "block_stats_csv": str(
                    case_dir / "block_stats.csv"
                ),
                "rank_by_geometry_csv": str(
                    case_dir / "rank_by_geometry.csv"
                ),
            },
        }

        write_json(case_dir / "result.json", record)

        print(
            "operator median: "
            f"direct={t_direct:.6f} s, "
            f"H={t_h:.6f} s, "
            f"speedup={operator['speedup_direct_over_hmatrix']:.3f}x, "
            f"rel.err={operator['relative_operator_error_l2']:.3e}",
            flush=True,
        )

        if break_even["exists"]:
            print(
                f"N_break={break_even['N_break']:.3f}, "
                "first profitable application="
                f"{break_even['first_profitable_application']}",
                flush=True,
            )
        else:
            print(
                "No break-even: "
                + break_even["reason"],
                flush=True,
            )

        return record


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "config_170ghz.nml",
        help="Production/base namelist. Only mesh counts are scaled.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "benchmark_scaling_results",
        help="Benchmark output directory.",
    )

    parser.add_argument(
        "--q",
        type=int,
        nargs="+",
        default=DEFAULT_TARGET_Q,
        help="Target Q values. Actual Q will be nearby integer-mesh values.",
    )

    parser.add_argument(
        "--repeat",
        type=int,
        default=3,
        help="Number of operator timing repetitions.",
    )

    parser.add_argument(
        "--warmup",
        type=int,
        default=0,
        help=(
            "Full-size untimed warmup calls for each operator. "
            "Default 0 because direct warmup can itself be expensive."
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed for deterministic complex benchmark current.",
    )

    parser.add_argument(
        "--hm-tol",
        type=float,
        default=1.0e-6,
    )
    parser.add_argument(
        "--hm-leaf",
        type=int,
        default=16,
    )
    parser.add_argument(
        "--hm-max-rank",
        type=int,
        default=32,
    )
    parser.add_argument(
        "--hm-eta",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "--hm-phase-limit",
        type=float,
        default=8.0,
    )
    parser.add_argument(
        "--hm-memory-mib",
        type=int,
        default=4096,
        help=(
            "H-matrix cache budget. This is not the whole-process RAM limit."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    base_config = args.config.resolve()
    output_dir = args.output.resolve()

    if not base_config.is_file():
        raise FileNotFoundError(base_config)

    if args.repeat < 1:
        raise ValueError("--repeat must be >= 1")
    if args.warmup < 0:
        raise ValueError("--warmup must be >= 0")
    if any(q < 1 for q in args.q):
        raise ValueError("All --q values must be positive")

    settings = Settings(
        tol=args.hm_tol,
        leaf=args.hm_leaf,
        max_rank=args.hm_max_rank,
        eta=args.hm_eta,
        phase_limit=args.hm_phase_limit,
        memory_mib=args.hm_memory_mib,
    )
    settings.validate()

    production = production_mesh_parameters(base_config)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(
        "Production mesh: "
        f"n_face={production['n_face']}, "
        f"n_z_cone={production['n_z_cone']}, "
        f"n_z_pipe={production['n_z_pipe']}, "
        f"Q={production['Q']:,}",
        flush=True,
    )
    print(
        "H settings: "
        + json.dumps(asdict(settings), sort_keys=True),
        flush=True,
    )
    print(
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, "
        f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}",
        flush=True,
    )

    records: list[dict] = []

    for target_q in args.q:
        try:
            record = run_one_case(
                target_q=target_q,
                base_config=base_config,
                production=production,
                settings=settings,
                output_dir=output_dir,
                repeats=args.repeat,
                warmup=args.warmup,
                seed=args.seed,
            )
        except KeyboardInterrupt:
            print(
                "\nInterrupted by user. Completed cases are preserved.",
                flush=True,
            )
            raise
        except Exception as exc:
            traceback.print_exc()
            record = {
                "status": "failed",
                "target_Q": int(target_q),
                "actual_Q": None,
                "mesh": scaled_mesh_parameters(
                    target_q,
                    production,
                ),
                "settings": asdict(settings),
                "error": f"{type(exc).__name__}: {exc}",
            }

            case_dir = (
                output_dir / f"Q_target_{target_q:06d}"
            )
            case_dir.mkdir(parents=True, exist_ok=True)
            write_json(case_dir / "result.json", record)

        records.append(record)

        # Commit the global summary after every case.
        write_json(
            output_dir / "scaling_results.json",
            records,
        )
        write_csv(
            output_dir / "scaling_summary.csv",
            [flatten_summary_record(x) for x in records],
        )

    ok = sum(record.get("status") == "ok" for record in records)
    failed = len(records) - ok

    print(
        "\nBenchmark finished: "
        f"{ok} successful, {failed} failed.",
        flush=True,
    )
    print(
        f"Summary JSON: {output_dir / 'scaling_results.json'}",
        flush=True,
    )
    print(
        f"Summary CSV : {output_dir / 'scaling_summary.csv'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
