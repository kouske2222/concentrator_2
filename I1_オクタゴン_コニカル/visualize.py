#!/usr/bin/env python3
"""
Visualize the G2 soft-octagon-to-circle cone geometry used for the closed-cap PO model.

This script does not run the electromagnetic solver. It only reproduces the wall-shape
formula so that the rounded octagon and its G2 connection to the circular pipe can be
checked quickly.

Typical use:
    python3 visualize_softoct_geometry.py --config mod_config.f90
    python3 visualize_softoct_geometry.py --config mod_config.f90 --beta 80 --out-dir results/260522_1530/geometry_check
"""
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  Needed by matplotlib for 3D projection.


@dataclass
class GeometryConfig:
    base_z: float = 1.0
    r_cone_in: float = 0.05
    r_pipe: float = 0.028
    l_cone: float = 0.10
    l_pipe: float = 0.30
    f_cap: float = 0.014
    N_z_cone: int = 350
    N_theta_max: int = 720

    @property
    def z_cone_start(self) -> float:
        return self.base_z

    @property
    def z_cone_end(self) -> float:
        return self.base_z + self.l_cone

    @property
    def z_pipe_end(self) -> float:
        return self.base_z + self.l_cone + self.l_pipe

    @property
    def z_cap_vertex(self) -> float:
        return self.z_pipe_end + self.r_pipe**2 / (4.0 * self.f_cap)


def _safe_eval_fortran_expr(expr: str, values: Dict[str, float]) -> float | None:
    """Evaluate simple numeric Fortran config expressions such as 0.05_dp or cfg%a**2."""
    expr = expr.split("!", 1)[0].strip()
    if not expr or "'" in expr or '"' in expr or ".true." in expr.lower() or ".false." in expr.lower():
        return None

    # Convert Fortran double/single precision suffixes to plain Python numbers.
    expr = re.sub(r"(?<=[0-9.])_[a-zA-Z][a-zA-Z0-9_]*", "", expr)
    expr = expr.replace("PI", "pi")
    expr = expr.replace("pi", "math.pi")

    # Replace cfg%name with already parsed numeric values.
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


def load_config(path: Path) -> GeometryConfig:
    cfg = GeometryConfig()
    if not path.exists():
        print(f"[WARN] Config file not found: {path}. Using built-in defaults.")
        return cfg

    values: Dict[str, float] = {}
    pattern = re.compile(r"cfg%([A-Za-z0-9_]+)\s*=\s*(.+)$")
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = pattern.search(raw)
        if not m:
            continue
        key, expr = m.group(1), m.group(2)
        val = _safe_eval_fortran_expr(expr, values)
        if val is None:
            continue
        values[key] = val

    for key in ("base_z", "r_cone_in", "r_pipe", "l_cone", "l_pipe", "f_cap"):
        if key in values:
            setattr(cfg, key, float(values[key]))
    for key in ("N_z_cone", "N_theta_max"):
        if key in values:
            setattr(cfg, key, int(round(values[key])))
    return cfg


def smootherstep(t: np.ndarray | float) -> np.ndarray | float:
    u = np.clip(t, 0.0, 1.0)
    return u**3 * (10.0 - 15.0*u + 6.0*u**2)


def smootherstep_d1(t: np.ndarray | float) -> np.ndarray | float:
    u = np.clip(t, 0.0, 1.0)
    ds = 30.0 * u**2 * (u - 1.0)**2
    return np.where((u <= 0.0) | (u >= 1.0), 0.0, ds)


def smootherstep_d2(t: np.ndarray | float) -> np.ndarray | float:
    u = np.clip(t, 0.0, 1.0)
    d2s = 60.0*u - 180.0*u**2 + 120.0*u**3
    return np.where((u <= 0.0) | (u >= 1.0), 0.0, d2s)


def smoothmax_octagon_cos(theta: np.ndarray, beta: float) -> np.ndarray:
    """Smooth approximation to max_k cos(theta - alpha_k)."""
    theta = np.asarray(theta)
    alphas = math.pi/8.0 + np.arange(8) * math.pi/4.0
    vals = beta * np.cos(theta[..., None] - alphas[None, ...])
    vmax = np.max(vals, axis=-1, keepdims=True)
    return ((vmax[..., 0] + np.log(np.sum(np.exp(vals - vmax), axis=-1))) / beta)


def smooth_octagon_radius(theta: np.ndarray, r_circ: float, beta: float) -> np.ndarray:
    # Normalize by theta=0 so that r(0)=r_circ. This preserves the intended vertex radius.
    m = smoothmax_octagon_cos(theta, beta)
    m_vertex = float(smoothmax_octagon_cos(np.array([0.0]), beta)[0])
    return r_circ * m_vertex / m


def blended_radius(theta: np.ndarray, t: np.ndarray, cfg: GeometryConfig, beta: float) -> np.ndarray:
    s = smootherstep(t)
    r8 = smooth_octagon_radius(theta, cfg.r_cone_in, beta)
    return (1.0 - s) * r8 + s * cfg.r_pipe


def blended_radius_dz(theta: np.ndarray, t: np.ndarray, cfg: GeometryConfig, beta: float) -> np.ndarray:
    r8 = smooth_octagon_radius(theta, cfg.r_cone_in, beta)
    return smootherstep_d1(t) * (cfg.r_pipe - r8) / cfg.l_cone


def blended_radius_dz2(theta: np.ndarray, t: np.ndarray, cfg: GeometryConfig, beta: float) -> np.ndarray:
    r8 = smooth_octagon_radius(theta, cfg.r_cone_in, beta)
    return smootherstep_d2(t) * (cfg.r_pipe - r8) / (cfg.l_cone**2)


def wall_radius_theta_z(theta: np.ndarray, z: np.ndarray, cfg: GeometryConfig, beta: float) -> np.ndarray:
    theta_arr, z_arr = np.broadcast_arrays(theta, z)
    r = np.full(theta_arr.shape, np.nan, dtype=float)

    cone = (z_arr >= cfg.z_cone_start) & (z_arr <= cfg.z_cone_end)
    pipe = (z_arr > cfg.z_cone_end) & (z_arr <= cfg.z_pipe_end)
    cap = (z_arr > cfg.z_pipe_end) & (z_arr <= cfg.z_cap_vertex)

    if np.any(cone):
        t = (z_arr[cone] - cfg.base_z) / cfg.l_cone
        r[cone] = blended_radius(theta_arr[cone], t, cfg, beta)
    if np.any(pipe):
        r[pipe] = cfg.r_pipe
    if np.any(cap):
        arg = cfg.r_pipe**2 - 4.0*cfg.f_cap*(z_arr[cap] - cfg.z_pipe_end)
        r[cap] = np.sqrt(np.maximum(arg, 0.0))
    return r


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"saved: {path}")


def plot_cross_sections(out_dir: Path, cfg: GeometryConfig, beta: float, ntheta: int) -> None:
    theta = np.linspace(0.0, 2.0*math.pi, ntheta, endpoint=True)
    t_values = [0.0, 0.05, 0.25, 0.50, 0.75, 0.95, 1.0]

    plt.figure(figsize=(7.2, 7.2))
    for t in t_values:
        r = blended_radius(theta, t, cfg, beta)
        plt.plot(1e3*r*np.cos(theta), 1e3*r*np.sin(theta), label=f"t={t:.2f}")
    circle = cfg.r_pipe * np.ones_like(theta)
    plt.plot(1e3*circle*np.cos(theta), 1e3*circle*np.sin(theta), linestyle="--", label="pipe circle")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.title("Cross-sections of soft octagon to circle transition")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    savefig(out_dir / "01_cross_sections_xy.svg")


def plot_xz_cuts(out_dir: Path, cfg: GeometryConfig, beta: float, nz: int) -> None:
    z = np.linspace(cfg.z_cone_start, cfg.z_pipe_end, nz)
    theta_list = [0.0, math.pi/16.0, math.pi/8.0]
    labels = ["vertex direction theta=0", "between", "side-center theta=pi/8"]

    plt.figure(figsize=(8.0, 5.0))
    for th, lab in zip(theta_list, labels):
        r = wall_radius_theta_z(np.full_like(z, th), z, cfg, beta)
        plt.plot(1e3*(z - cfg.base_z), 1e3*r, label=lab)
        plt.plot(1e3*(z - cfg.base_z), -1e3*r, linestyle="--")
    plt.axvline(1e3*cfg.l_cone, linestyle=":", label="cone-pipe joint")
    plt.xlabel("z - base_z [mm]")
    plt.ylabel("wall coordinate x or radial cut [mm]")
    plt.title("Longitudinal cuts: G2 transition into circular pipe")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    savefig(out_dir / "02_xz_longitudinal_cuts.svg")


def plot_radius_heatmap(out_dir: Path, cfg: GeometryConfig, beta: float, ntheta: int, nz: int) -> None:
    theta = np.linspace(0.0, 2.0*math.pi, ntheta)
    z = np.linspace(cfg.z_cone_start, cfg.z_cone_end, nz)
    TH, ZZ = np.meshgrid(theta, z)
    TT = (ZZ - cfg.base_z) / cfg.l_cone
    RR = blended_radius(TH, TT, cfg, beta)

    plt.figure(figsize=(8.5, 4.8))
    pc = plt.pcolormesh(1e3*(z - cfg.base_z), theta/math.pi, 1e3*RR.T, shading="auto", rasterized=True)
    plt.xlabel("z - base_z [mm]")
    plt.ylabel("theta / pi [-]")
    plt.title("Radius map in cone section")
    cb = plt.colorbar(pc)
    cb.set_label("r(theta,z) [mm]")
    savefig(out_dir / "03_radius_theta_z_heatmap.svg")


def plot_g2_derivatives(out_dir: Path, cfg: GeometryConfig, beta: float, nz: int) -> Tuple[float, float]:
    eps = 1e-9
    t = np.linspace(eps, 1.0 - eps, nz)
    zmm = 1e3 * cfg.l_cone * t
    theta_list = [0.0, math.pi/8.0]
    labels = ["vertex direction theta=0", "side-center theta=pi/8"]

    max_abs_dz_at_ends = 0.0
    max_abs_dz2_at_ends = 0.0

    plt.figure(figsize=(8.0, 5.0))
    for th, lab in zip(theta_list, labels):
        theta = np.full_like(t, th)
        drdz = blended_radius_dz(theta, t, cfg, beta)
        plt.plot(zmm, drdz, label=lab)
        for te in [0.0, 1.0]:
            max_abs_dz_at_ends = max(max_abs_dz_at_ends, abs(float(blended_radius_dz(np.array([th]), np.array([te]), cfg, beta)[0])))
    plt.axhline(0.0, linestyle=":")
    plt.xlabel("z - base_z [mm]")
    plt.ylabel("dr/dz [-]")
    plt.title("First derivative check: slope goes to zero at both ends")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    savefig(out_dir / "04_g2_check_first_derivative.svg")

    plt.figure(figsize=(8.0, 5.0))
    for th, lab in zip(theta_list, labels):
        theta = np.full_like(t, th)
        d2rdz2 = blended_radius_dz2(theta, t, cfg, beta)
        plt.plot(zmm, d2rdz2, label=lab)
        for te in [0.0, 1.0]:
            max_abs_dz2_at_ends = max(max_abs_dz2_at_ends, abs(float(blended_radius_dz2(np.array([th]), np.array([te]), cfg, beta)[0])))
    plt.axhline(0.0, linestyle=":")
    plt.xlabel("z - base_z [mm]")
    plt.ylabel("d2r/dz2 [1/m]")
    plt.title("Second derivative check: curvature goes to zero at both ends")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    savefig(out_dir / "05_g2_check_second_derivative.svg")

    return max_abs_dz_at_ends, max_abs_dz2_at_ends


def plot_3d_surface(out_dir: Path, cfg: GeometryConfig, beta: float, ntheta: int, nz: int) -> None:
    theta = np.linspace(0.0, 2.0*math.pi, ntheta)
    z = np.linspace(cfg.z_cone_start, cfg.z_pipe_end, nz)
    TH, ZZ = np.meshgrid(theta, z)
    RR = wall_radius_theta_z(TH, ZZ, cfg, beta)
    X = RR * np.cos(TH)
    Y = RR * np.sin(TH)

    fig = plt.figure(figsize=(8.2, 6.6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(1e3*X, 1e3*Y, 1e3*(ZZ - cfg.base_z), rstride=max(1, nz//30), cstride=max(1, ntheta//48), linewidth=0.5)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z - base_z [mm]")
    ax.set_title("3D wireframe: soft octagon cone connected to circular pipe")
    ax.set_box_aspect((2*cfg.r_cone_in, 2*cfg.r_cone_in, cfg.l_cone + cfg.l_pipe))
    savefig(out_dir / "06_3d_wireframe_cone_pipe.svg")


def plot_panel_like_sampling(out_dir: Path, cfg: GeometryConfig, beta: float) -> None:
    # This imitates the adaptive theta sampling used by the Fortran panel generator.
    ds_target = 2.0 * math.pi * cfg.r_cone_in / float(cfg.N_theta_max)
    dz = cfg.l_cone / float(cfg.N_z_cone)
    xs, ys, zs = [], [], []

    # Keep the plot light even when N_z_cone is large.
    stride_z = max(1, cfg.N_z_cone // 70)
    for i in range(0, cfg.N_z_cone, stride_z):
        zc = (i + 0.5) * dz
        t = zc / cfg.l_cone
        r_envelope = (1.0 - smootherstep(t)) * cfg.r_cone_in + smootherstep(t) * cfg.r_pipe
        n_th = max(int(round(2.0 * math.pi * r_envelope / ds_target)), 12)
        stride_th = max(1, n_th // 96)
        dtheta = 2.0 * math.pi / n_th
        for j in range(0, n_th, stride_th):
            th = (j + 0.5) * dtheta
            r = blended_radius(np.array([th]), np.array([t]), cfg, beta)[0]
            xs.append(r * math.cos(th))
            ys.append(r * math.sin(th))
            zs.append(cfg.base_z + zc)

    fig = plt.figure(figsize=(8.0, 6.6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(1e3*np.array(xs), 1e3*np.array(ys), 1e3*(np.array(zs) - cfg.base_z), s=1.0)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z - base_z [mm]")
    ax.set_title("Panel-center-like sampling of cone section")
    ax.set_box_aspect((2*cfg.r_cone_in, 2*cfg.r_cone_in, cfg.l_cone))
    savefig(out_dir / "07_panel_like_sampling_cone.svg")


def write_summary(out_dir: Path, cfg: GeometryConfig, beta: float, max_d1: float, max_d2: float) -> None:
    theta = np.linspace(0.0, 2.0*math.pi, 4097)
    r0 = smooth_octagon_radius(theta, cfg.r_cone_in, beta)
    sharp_inradius = cfg.r_cone_in * math.cos(math.pi/8.0)

    text = f"""Geometry visualization summary
================================

Formula:
  r(theta,t) = (1 - s(t)) r8_soft(theta) + s(t) r_pipe
  s(t)      = 6 t^5 - 15 t^4 + 10 t^3

Soft octagon:
  m_beta(theta) = (1/beta) log sum_k exp(beta cos(theta - alpha_k))
  alpha_k       = pi/8 + k pi/4,  k=0..7
  r8_soft(theta)= r_cone_in * m_beta(0) / m_beta(theta)

Parameters:
  beta       = {beta:.6g}
  base_z     = {cfg.base_z:.9g} m
  r_cone_in  = {cfg.r_cone_in:.9g} m
  r_pipe     = {cfg.r_pipe:.9g} m
  l_cone     = {cfg.l_cone:.9g} m
  l_pipe     = {cfg.l_pipe:.9g} m
  f_cap      = {cfg.f_cap:.9g} m
  N_z_cone   = {cfg.N_z_cone:d}
  N_theta_max= {cfg.N_theta_max:d}

Octagon check at t=0:
  max radius                         = {np.nanmax(r0):.9e} m
  min radius                         = {np.nanmin(r0):.9e} m
  min/max ratio of soft octagon       = {np.nanmin(r0)/np.nanmax(r0):.9f}
  sharp regular-octagon inradius/R    = {math.cos(math.pi/8.0):.9f}
  sharp regular-octagon inradius      = {sharp_inradius:.9e} m

G2 endpoint check:
  max |dr/dz| at t=0 or t=1 among checked directions     = {max_d1:.9e}
  max |d2r/dz2| at t=0 or t=1 among checked directions   = {max_d2:.9e} 1/m

Interpretation:
  01_cross_sections_xy.svg checks whether the inlet is a rounded octagon and gradually becomes circular.
  02_xz_longitudinal_cuts.svg checks whether vertex and side-center directions join the pipe smoothly.
  04 and 05 check the G2 condition: slope and second derivative vanish at both ends.
"""
    path = out_dir / "geometry_check_summary.txt"
    path.write_text(text, encoding="utf-8")
    print(f"saved: {path}")


def default_output_dir() -> Path:
    stamp = datetime.now().strftime("%y%m%d_%H%M")
    return Path("results") / f"{stamp}_geometry_check"


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize the G2 soft-octagon-to-circle concentrator geometry.")
    parser.add_argument("--config", default="mod_config.f90", help="Fortran config file to read. Default: mod_config.f90")
    parser.add_argument("--out-dir", default=None, help="Output directory. Default: results/YYMMDD_HHMM_geometry_check")
    parser.add_argument("--beta", type=float, default=80.0, help="Softmax beta. Larger means closer to a sharp octagon. Default: 80")
    parser.add_argument("--ntheta", type=int, default=721, help="Number of theta samples for smooth plots. Default: 721")
    parser.add_argument("--nz", type=int, default=401, help="Number of z samples for smooth plots. Default: 401")
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    out_dir = Path(args.out_dir) if args.out_dir else default_output_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Geometry visualization")
    print(f"  config : {Path(args.config).resolve()}")
    print(f"  out_dir: {out_dir.resolve()}")
    print(f"  beta   : {args.beta}")

    plot_cross_sections(out_dir, cfg, args.beta, args.ntheta)
    plot_xz_cuts(out_dir, cfg, args.beta, args.nz)
    plot_radius_heatmap(out_dir, cfg, args.beta, args.ntheta, args.nz)
    max_d1, max_d2 = plot_g2_derivatives(out_dir, cfg, args.beta, args.nz)
    plot_3d_surface(out_dir, cfg, args.beta, min(args.ntheta, 181), min(args.nz, 161))
    plot_panel_like_sampling(out_dir, cfg, args.beta)
    write_summary(out_dir, cfg, args.beta, max_d1, max_d2)

    print("done")


if __name__ == "__main__":
    main()
