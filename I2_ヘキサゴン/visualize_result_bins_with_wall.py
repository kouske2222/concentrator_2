#!/usr/bin/env python3
"""
Visualize results/<timestamp>/field_map*.bin with theoretical wall-shape overlay.

Purpose:
  - Distinguish whether the electromagnetic field merely looks circular
    from whether the wall/mask shape is actually circular.
  - For XY maps, draw the expected soft-hexagon / circle wall boundary at each z.
  - For XZ maps, draw x = +/- r_wall(theta=0 or pi, z).

Expected Fortran binary format:
  obs_x(:), obs_y(:), obs_z(:), e_sq(:)
as float64 stream data, no header.

Examples:
  python3 visualize_result_bins_with_wall.py --latest --config mod_config.f90
  python3 visualize_result_bins_with_wall.py --dir results/260522_1530 --config mod_config.f90
  python3 visualize_result_bins_with_wall.py --latest --config mod_config.f90 --beta 200
  python3 visualize_result_bins_with_wall.py --latest --shape-only
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

POLYGON_N_SIDES = 6
POLYGON_NAME = "hexagon"


@dataclass
class Config:
    base_z: float = 1.0
    l_cone: float = 0.10
    l_pipe: float = 0.40
    r_cone_in: float = 0.05
    r_pipe: float = 0.028
    f_cap: float = 0.050
    field_outer_margin_factor_xy: float = 0.0
    beta: float = 80.0


@dataclass
class FieldMap:
    path: Path
    kind: str
    npts: int
    nx: int
    ny_or_nz: int
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    e2: np.ndarray
    X: np.ndarray
    Y_or_Z: np.ndarray
    E2: np.ndarray
    const_value: float
    e2_max: float


def natural_key(path: Path) -> Tuple:
    parts = re.split(r"(\d+)", path.name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _safe_eval_fortran_expr(expr: str, values: Dict[str, float]) -> float | None:
    """Evaluate simple numeric Fortran config expressions such as 0.05_dp or cfg%d_in."""
    expr = expr.split("!", 1)[0].strip()
    if not expr or "'" in expr or '"' in expr or ".true." in expr.lower() or ".false." in expr.lower():
        return None

    expr = re.sub(r"(?<=[0-9.])_[a-zA-Z][a-zA-Z0-9_]*", "", expr)
    expr = re.sub(r"(?<=[0-9])[dD](?=[+-]?[0-9])", "e", expr)
    expr = expr.replace("PI", "math.pi")

    def repl(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in values:
            raise KeyError(key)
        return str(values[key])

    try:
        expr = re.sub(r"cfg%([A-Za-z0-9_]+)", repl, expr)
        return float(eval(expr, {"__builtins__": {}}, {"math": math, "sqrt": math.sqrt}))
    except Exception:
        return None


def read_config(path: Path) -> Config:
    cfg = Config()
    if not path.exists():
        print(f"[warn] config file not found: {path}. Using defaults.")
        return cfg

    values: Dict[str, float] = {}
    pattern = re.compile(r"cfg%([A-Za-z0-9_]+)\s*=\s*(.+)$", flags=re.IGNORECASE)
    for raw in path.read_text(errors="ignore").splitlines():
        m = pattern.search(raw)
        if not m:
            continue
        key, expr = m.group(1), m.group(2)
        val = _safe_eval_fortran_expr(expr, values)
        if val is not None:
            values[key] = val

    for name in (
        "base_z",
        "l_cone",
        "l_pipe",
        "r_cone_in",
        "r_pipe",
        "f_cap",
        "field_outer_margin_factor_xy",
    ):
        if name in values:
            setattr(cfg, name, float(values[name]))
    return cfg

    text = path.read_text(errors="ignore")
    patterns = {
        "base_z": r"cfg%base_z\s*=\s*([^\n!]+)",
        "l_cone": r"cfg%l_cone\s*=\s*([^\n!]+)",
        "l_pipe": r"cfg%l_pipe\s*=\s*([^\n!]+)",
        "r_cone_in": r"cfg%r_cone_in\s*=\s*([^\n!]+)",
        "r_pipe": r"cfg%r_pipe\s*=\s*([^\n!]+)",
        "f_cap": r"cfg%f_cap\s*=\s*([^\n!]+)",
        "field_outer_margin_factor_xy": r"cfg%field_outer_margin_factor_xy\s*=\s*([^\n!]+)",
    }
    for name, pat in patterns.items():
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            try:
                setattr(cfg, name, parse_fortran_float(m.group(1)))
            except Exception as exc:
                print(f"[warn] failed to parse {name}: {m.group(1)!r} ({exc})")
    return cfg


def unique_count(arr: np.ndarray, decimals: int = 12) -> Tuple[np.ndarray, int]:
    vals = np.unique(np.round(arr[np.isfinite(arr)], decimals=decimals))
    return vals, vals.size


def finite_positive_max(arr: np.ndarray) -> float:
    vals = arr[np.isfinite(arr) & (arr > 0.0)]
    if vals.size == 0:
        return float("nan")
    return float(np.max(vals))


def read_field_map_binary(path: Path) -> FieldMap:
    raw = np.fromfile(path, dtype=np.float64)
    if raw.size == 0:
        raise ValueError(f"empty file: {path}")
    if raw.size % 4 != 0:
        raise ValueError(f"{path}: data length {raw.size} is not divisible by 4")

    npts = raw.size // 4
    x = raw[:npts]
    y = raw[npts:2*npts]
    z = raw[2*npts:3*npts]
    e2 = raw[3*npts:]

    ux, nx = unique_count(x)
    uy, ny = unique_count(y)
    uz, nz = unique_count(z)

    if ny == 1 and nx * nz == npts:
        X = x.reshape((nz, nx))
        Z = z.reshape((nz, nx))
        E2 = e2.reshape((nz, nx))
        return FieldMap(path, "XZ", npts, nx, nz, x, y, z, e2, X, Z, E2, float(uy[0]), finite_positive_max(e2))

    if nz == 1 and nx * ny == npts:
        X = x.reshape((ny, nx))
        Y = y.reshape((ny, nx))
        E2 = e2.reshape((ny, nx))
        return FieldMap(path, "XY", npts, nx, ny, x, y, z, e2, X, Y, E2, float(uz[0]), finite_positive_max(e2))

    raise ValueError(f"{path}: cannot infer grid: nx={nx}, ny={ny}, nz={nz}, npts={npts}")


def find_latest_results_dir(root: Path) -> Path:
    if not root.exists():
        raise FileNotFoundError(f"results root not found: {root}")
    candidates = [p for p in root.iterdir() if p.is_dir() and any(p.glob("field_map*.bin"))]
    if not candidates:
        raise FileNotFoundError(f"no directories with field_map*.bin under {root}")
    timestamp_like = [p for p in candidates if re.match(r"^\d{6}_\d{4}$", p.name)]
    if timestamp_like:
        return sorted(timestamp_like, key=lambda p: p.name)[-1]
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def collect_bin_files(result_dir: Path) -> List[Path]:
    files = sorted(result_dir.glob("field_map*.bin"), key=natural_key)
    if not files:
        raise FileNotFoundError(f"no field_map*.bin found in {result_dir}")
    return files


def smoothstep_g2(t: np.ndarray | float) -> np.ndarray | float:
    tt = np.clip(t, 0.0, 1.0)
    return 6.0*tt**5 - 15.0*tt**4 + 10.0*tt**3


def softmax_polygon_cos(theta: np.ndarray | float, beta: float) -> np.ndarray | float:
    theta_arr = np.asarray(theta, dtype=float)
    alphas = math.pi/POLYGON_N_SIDES + np.arange(POLYGON_N_SIDES, dtype=float)*2.0*math.pi/POLYGON_N_SIDES
    a = beta * np.cos(theta_arr[..., None] - alphas[None, :])
    amax = np.max(a, axis=-1)
    return (amax + np.log(np.sum(np.exp(a - amax[..., None]), axis=-1))) / beta


def soft_polygon_radius(theta: np.ndarray | float, r_vertex: float, beta: float) -> np.ndarray | float:
    m_vertex = softmax_polygon_cos(0.0, beta)
    m_theta = softmax_polygon_cos(theta, beta)
    return r_vertex * m_vertex / m_theta


def wall_radius_theta_z(theta: np.ndarray | float, z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    theta_arr = np.asarray(theta, dtype=float)
    z_arr = np.asarray(z, dtype=float)
    theta_b, z_b = np.broadcast_arrays(theta_arr, z_arr)
    rw = np.full(theta_b.shape, np.nan, dtype=float)

    z_pipe_start = cfg.base_z + cfg.l_cone
    z_pipe_end = cfg.base_z + cfg.l_cone + cfg.l_pipe
    z_cap_vertex = z_pipe_end + cfg.r_pipe**2 / (4.0*cfg.f_cap)

    cone = (z_b >= cfg.base_z) & (z_b <= z_pipe_start)
    pipe = (z_b > z_pipe_start) & (z_b <= z_pipe_end)
    cap = (z_b > z_pipe_end) & (z_b <= z_cap_vertex)

    if np.any(cone):
        t = (z_b[cone] - cfg.base_z) / cfg.l_cone
        s = smoothstep_g2(t)
        rpoly = soft_polygon_radius(theta_b[cone], cfg.r_cone_in, cfg.beta)
        rw[cone] = (1.0 - s)*rpoly + s*cfg.r_pipe
    if np.any(pipe):
        rw[pipe] = cfg.r_pipe
    if np.any(cap):
        arg = cfg.r_pipe**2 - 4.0*cfg.f_cap*(z_b[cap] - z_pipe_end)
        rw[cap] = np.sqrt(np.maximum(arg, 0.0))

    if np.isscalar(theta) and np.isscalar(z):
        return float(rw)
    return rw


def to_db(E2: np.ndarray, ref: float, floor_db: float) -> np.ndarray:
    C = np.full_like(E2, np.nan, dtype=float)
    mask = np.isfinite(E2) & (E2 > 0.0) & np.isfinite(ref) & (ref > 0.0)
    C[mask] = 10.0*np.log10(E2[mask]/ref)
    C = np.maximum(C, floor_db)
    return C


def normalized_linear(E2: np.ndarray, ref: float) -> np.ndarray:
    C = np.full_like(E2, np.nan, dtype=float)
    mask = np.isfinite(E2) & np.isfinite(ref) & (ref > 0.0)
    C[mask] = E2[mask]/ref
    return C


def classify_z(z: float, cfg: Config) -> str:
    t = (z - cfg.base_z) / cfg.l_cone
    if z < cfg.base_z:
        return "upstream/outside"
    if 0.0 <= t <= 1.0:
        return f"cone, t={t:.4f}"
    if z <= cfg.base_z + cfg.l_cone + cfg.l_pipe:
        return "pipe/circular"
    return "cap/downstream"


def plot_shape_only(cfg: Config, out_dir: Path, file_format: str = "svg", dpi: int = 200) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    theta = np.linspace(0.0, 2.0*math.pi, 721)
    z_list = [
        cfg.base_z,
        cfg.base_z + 0.25*cfg.l_cone,
        cfg.base_z + 0.50*cfg.l_cone,
        cfg.base_z + 0.75*cfg.l_cone,
        cfg.base_z + cfg.l_cone,
    ]

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    for z in z_list:
        r = wall_radius_theta_z(theta, z, cfg)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        ax.plot(x*1e3, y*1e3, label=f"z={z:.4f} m, {classify_z(z, cfg)}")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title("Theoretical wall cross sections")
    ax.grid(True)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / f"shape_only_xy_sections.{file_format}", dpi=dpi)
    plt.close(fig)

    # radius vs angle at cone entrance
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    r0 = wall_radius_theta_z(theta, cfg.base_z, cfg)
    ax.plot(theta, r0*1e3)
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(r"$r_{wall}$ [mm]")
    ax.set_title("Entrance soft-hexagon radius vs angle")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out_dir / f"shape_only_radius_vs_theta.{file_format}", dpi=dpi)
    plt.close(fig)


def plot_one_map_with_wall(
    fmap: FieldMap,
    out_path: Path,
    cfg: Config,
    mode: str,
    ref: float,
    db_floor: float,
    cmap: str,
    file_format: str,
    dpi: int,
) -> None:
    if mode == "db":
        C = to_db(fmap.E2, ref, db_floor)
        label = r"$10\log_{10}(|E|^2/|E|^2_{max})$ [dB]"
        vmin, vmax = db_floor, 0.0
    elif mode == "linear":
        C = normalized_linear(fmap.E2, ref)
        label = r"$|E|^2/|E|^2_{max}$ [-]"
        vmin, vmax = 0.0, 1.0
    elif mode == "raw":
        C = fmap.E2
        vals = C[np.isfinite(C)]
        vmin = float(np.nanmin(vals)) if vals.size else None
        vmax = float(np.nanmax(vals)) if vals.size else None
        label = r"$|E|^2$"
    else:
        raise ValueError(mode)

    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    im = ax.pcolormesh(fmap.X, fmap.Y_or_Z, C, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(label)

    if fmap.kind == "XY":
        z0 = fmap.const_value
        th = np.linspace(0.0, 2.0*math.pi, 721)
        rw = wall_radius_theta_z(th, z0, cfg)
        # Solid black wall boundary. Also draw margin boundary if margin is nonzero.
        ax.plot(rw*np.cos(th), rw*np.sin(th), color="black", linewidth=1.8, label="wall")
        if cfg.field_outer_margin_factor_xy != 0.0:
            rm = (1.0 + cfg.field_outer_margin_factor_xy)*rw
            ax.plot(rm*np.cos(th), rm*np.sin(th), color="black", linewidth=0.8, linestyle="--", label="wall + XY margin")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"{fmap.path.name}: XY, z={z0:.6g} m ({classify_z(z0, cfg)})")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="upper right", fontsize=8)
    elif fmap.kind == "XZ":
        zline = np.linspace(np.nanmin(fmap.z), np.nanmax(fmap.z), 900)
        r_plus = wall_radius_theta_z(0.0, zline, cfg)
        r_minus = wall_radius_theta_z(math.pi, zline, cfg)
        ax.plot(r_plus, zline, color="black", linewidth=1.6, label="wall +x")
        ax.plot(-r_minus, zline, color="black", linewidth=1.6, label="wall -x")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("z [m]")
        ax.set_title(f"{fmap.path.name}: XZ, y={fmap.const_value:.4g} m")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="upper right", fontsize=8)
    else:
        ax.set_xlabel("x [m]")
        ax.set_ylabel("coordinate [m]")
        ax.set_title(fmap.path.name)

    ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(f".{file_format}"), dpi=dpi)
    plt.close(fig)


def write_wall_diagnostics(maps: Sequence[FieldMap], cfg: Config, out_dir: Path) -> None:
    theta_vertex = 0.0
    theta_side = math.pi/POLYGON_N_SIDES
    with (out_dir / "wall_shape_diagnostics.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "kind", "z_or_y", "region", "r_vertex_m", "r_side_m", "contrast_percent"])
        for m in maps:
            if m.kind != "XY":
                continue
            z0 = m.const_value
            rv = wall_radius_theta_z(theta_vertex, z0, cfg)
            rs = wall_radius_theta_z(theta_side, z0, cfg)
            contrast = 100.0*(rv - rs)/rv if np.isfinite(rv) and rv != 0.0 else np.nan
            w.writerow([m.path.name, m.kind, f"{z0:.12e}", classify_z(z0, cfg), f"{rv:.12e}", f"{rs:.12e}", f"{contrast:.6f}"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize field_map*.bin with soft-hexagon wall overlay.")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--dir", type=Path, help="result directory, e.g. results/260522_1530")
    g.add_argument("--latest", action="store_true", help="use latest results/<timestamp> directory")
    parser.add_argument("--root", type=Path, default=Path("results"), help="root used by --latest")
    parser.add_argument("--config", type=Path, default=Path("mod_config.f90"), help="Fortran config file")
    parser.add_argument("--out", type=Path, default=None, help="output directory")
    parser.add_argument("--beta", type=float, default=None, help="override softmax beta")
    parser.add_argument("--shape-only", action="store_true", help="plot theoretical wall shapes only")
    parser.add_argument("--plot", choices=["db", "linear", "raw", "both", "all"], default="both")
    parser.add_argument("--local-scale", action="store_true", help="normalize each map by local maximum")
    parser.add_argument("--db-floor", type=float, default=-60.0)
    parser.add_argument("--cmap", default="rainbow")
    parser.add_argument("--format", default="svg", choices=["svg", "png", "pdf"])
    parser.add_argument("--dpi", type=int, default=200)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = read_config(args.config)
    if args.beta is not None:
        cfg.beta = args.beta

    if args.shape_only:
        out_dir = args.out if args.out is not None else Path("geometry_wall_check")
        print("Config:", cfg)
        print(f"Polygon sides: {POLYGON_N_SIDES} ({POLYGON_NAME})")
        plot_shape_only(cfg, out_dir, args.format, args.dpi)
        print(f"Shape-only figures saved in: {out_dir}")
        return

    result_dir = args.dir if args.dir is not None else find_latest_results_dir(args.root)
    if not result_dir.exists():
        raise FileNotFoundError(f"result directory not found: {result_dir}")
    out_dir = args.out if args.out is not None else result_dir / "bin_visualization_with_wall"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = collect_bin_files(result_dir)
    maps = [read_field_map_binary(p) for p in files]

    print("Config:", cfg)
    print(f"Polygon sides: {POLYGON_N_SIDES} ({POLYGON_NAME})")
    print(f"Result directory : {result_dir}")
    print(f"Output directory : {out_dir}")
    print(f"Found {len(maps)} maps")
    for m in maps:
        if m.kind == "XY":
            rv = wall_radius_theta_z(0.0, m.const_value, cfg)
            rs = wall_radius_theta_z(math.pi/POLYGON_N_SIDES, m.const_value, cfg)
            contrast = 100.0*(rv-rs)/rv if np.isfinite(rv) and rv != 0.0 else np.nan
            print(f"  {m.path.name}: XY z={m.const_value:.6f} m, {classify_z(m.const_value, cfg)}, wall contrast={contrast:.3f}%")
        else:
            print(f"  {m.path.name}: {m.kind}, grid={m.nx} x {m.ny_or_nz}, max={m.e2_max:.6e}")

    all_e2 = np.concatenate([m.e2 for m in maps])
    global_ref = finite_positive_max(all_e2)
    if not np.isfinite(global_ref) or global_ref <= 0.0:
        raise ValueError("no positive finite |E|^2 values found")

    modes = ["db", "linear"] if args.plot == "both" else ["db", "linear", "raw"] if args.plot == "all" else [args.plot]

    for m in maps:
        for mode in modes:
            ref = m.e2_max if args.local_scale or mode == "raw" else global_ref
            scale_tag = "local" if args.local_scale or mode == "raw" else "global"
            out_stem = out_dir / f"{m.path.stem}_{mode}_{scale_tag}_wall"
            plot_one_map_with_wall(m, out_stem, cfg, mode, ref, args.db_floor, args.cmap, args.format, args.dpi)

    plot_shape_only(cfg, out_dir, args.format, args.dpi)
    write_wall_diagnostics(maps, cfg, out_dir)
    print("Done.")
    print(f"Figures saved in: {out_dir}")


if __name__ == "__main__":
    main()
