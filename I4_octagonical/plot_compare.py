#!/usr/bin/env python3
"""Create common-scale full/C8 comparison figures from Fortran CSV output."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_csv(path: Path) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", names=True)


def field_components(data: np.ndarray) -> list[np.ndarray]:
    ex = data["Ex_re"] + 1j * data["Ex_im"]
    ey = data["Ey_re"] + 1j * data["Ey_im"]
    ez = data["Ez_re"] + 1j * data["Ez_im"]
    return [np.abs(ex), np.abs(ey), np.abs(ez), data["E_mag"]]


def db(values: np.ndarray, reference: float) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(values / max(reference, 1e-300), 1e-6))


def reshape_plane(data: np.ndarray, values: np.ndarray, horizontal: str, vertical: str):
    a = np.unique(data[horizontal])
    b = np.unique(data[vertical])
    return a, b, values.reshape(len(b), len(a))


def metrics(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def create_figures(result: Path) -> None:
    out = result / "figures"
    out.mkdir(exist_ok=True)
    full = load_csv(result / "xz_full.csv")
    modal = load_csv(result / "xz_modal.csv")
    vf = field_components(full)
    vm = field_components(modal)
    reference = np.nanmax(np.where(full["mask"] > 0, vf[-1], np.nan))
    labels = ["|Ex|", "|Ey|", "|Ez|", "|E|"]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
    for column, label in enumerate(labels):
        for row, (data, values, method) in enumerate(((full, vf[column], "full"), (modal, vm[column], "C8 modal"))):
            x, z, image = reshape_plane(data, np.where(data["mask"] > 0, values, np.nan), "x_m", "z_m")
            im = axes[row, column].imshow(db(image, reference), origin="lower", aspect="auto",
                                          extent=[x[0]*1e3, x[-1]*1e3, z[0]*1e3, z[-1]*1e3],
                                          vmin=-60, vmax=0, cmap="magma")
            axes[row, column].set_title(f"{method} {label}")
            axes[row, column].set_xlabel("x [mm]")
            if column == 0:
                axes[row, column].set_ylabel("z [mm]")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="normalized magnitude [dB]")
    fig.savefig(out / "xz_full_vs_c8.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    order_path = result / "order_comparison.csv"
    if order_path.exists():
        order = load_csv(order_path)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].semilogy(order["order"], order["full_ratio"], "o-", label="full")
        axes[0].semilogy(order["order"], order["modal_ratio"], "x--", label="all modes")
        axes[0].semilogy(order["order"], order["fast_ratio"], "+:", label="m=1,7")
        axes[0].set(xlabel="reflection order", ylabel=r"$||J^{(r)}||/||J_{cum}^{(r)}||$")
        axes[0].grid(True, which="both", alpha=.3); axes[0].legend()
        axes[1].semilogy(order["order"], order["current_rel_error"], "o-")
        axes[1].axhline(1e-10, color="k", ls="--", lw=.8)
        axes[1].set(xlabel="reflection order", ylabel="full-modal relative error")
        axes[1].grid(True, which="both", alpha=.3)
        fig.tight_layout(); fig.savefig(out / "order_comparison.png", dpi=180); plt.close(fig)

    modes = load_csv(result / "mode_norms.csv")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(modes["mode"] - .15, modes["initial_norm"], width=.3, label="initial")
    ax.bar(modes["mode"] + .15, modes["final_norm"], width=.3, label="last order")
    ax.set(xlabel="C8 mode m", ylabel="DFT current norm")
    ax.legend(); fig.tight_layout(); fig.savefig(out / "mode_norms.png", dpi=180); plt.close(fig)

    af = load_csv(result / "axis_full.csv")
    am = load_csv(result / "axis_modal.csv")
    exf = af["Ex_re"] + 1j * af["Ex_im"]
    exm = am["Ex_re"] + 1j * am["Ex_im"]
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(af["z_m"]*1e3, af["E_mag"], label="full")
    axes[0].plot(am["z_m"]*1e3, am["E_mag"], "--", label="C8")
    axes[0].set_ylabel("|E| [V/m]"); axes[0].legend()
    axes[1].plot(af["z_m"]*1e3, np.unwrap(np.angle(exf)), label="full")
    axes[1].plot(am["z_m"]*1e3, np.unwrap(np.angle(exm)), "--", label="C8")
    axes[1].set(xlabel="z [mm]", ylabel="phase(Ex) [rad]")
    fig.tight_layout(); fig.savefig(out / "axis_magnitude_phase.png", dpi=180); plt.close(fig)

    current = load_csv(result / "surface_current.csv")
    vmax = max(current["Jfull_mag"].max(), current["Jmodal_mag"].max())
    fig = plt.figure(figsize=(11, 4))
    for index, (key, title) in enumerate((("Jfull_mag", "full"), ("Jmodal_mag", "C8 modal")), 1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        sc = ax.scatter(current["x_m"]*1e3, current["y_m"]*1e3, current["z_m"]*1e3,
                        c=current[key], s=8, vmin=0, vmax=vmax, cmap="viridis")
        ax.set(title=title, xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]")
        fig.colorbar(sc, ax=ax, shrink=.65, label="|J| [A/m]")
    fig.tight_layout(); fig.savefig(out / "surface_current.png", dpi=180); plt.close(fig)

    for plane in range(1, 4):
        pf = load_csv(result / f"xy{plane}_full.csv")
        pm = load_csv(result / f"xy{plane}_modal.csv")
        ref = np.nanmax(np.where(pf["mask"] > 0, pf["E_mag"], np.nan))
        error = np.abs(pm["E_mag"] - pf["E_mag"]) / np.maximum(pf["E_mag"], 1e-300)
        fig = plt.figure(figsize=(13.2, 4.0), constrained_layout=True)
        gs = fig.add_gridspec(1, 5, width_ratios=[1.0, 1.0, 0.045, 1.0, 0.045])
        axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 3])]
        cax_mag = fig.add_subplot(gs[0, 2])
        cax_err = fig.add_subplot(gs[0, 4])
        for ax, data, title in ((axes[0], pf, "full"), (axes[1], pm, "C8 modal")):
            x, y, image = reshape_plane(data, np.where(data["mask"] > 0, data["E_mag"], np.nan), "x_m", "y_m")
            im = ax.imshow(db(image, ref), origin="lower", extent=[x[0]*1e3,x[-1]*1e3,y[0]*1e3,y[-1]*1e3],
                           vmin=-60, vmax=0, cmap="magma")
            ax.set(title=title, xlabel="x [mm]", ylabel="y [mm]")
            ax.set_aspect("equal", adjustable="box")
        fig.colorbar(im, cax=cax_mag, label="normalized |E| [dB]")
        x, y, image = reshape_plane(pf, np.where(pf["mask"] > 0, error, np.nan), "x_m", "y_m")
        err_im = axes[2].imshow(np.log10(np.maximum(image, 1e-16)), origin="lower",
                                extent=[x[0]*1e3,x[-1]*1e3,y[0]*1e3,y[-1]*1e3],
                                vmin=-16, vmax=-8, cmap="viridis")
        axes[2].set(title="relative |E| error", xlabel="x [mm]", ylabel="y [mm]")
        axes[2].set_aspect("equal", adjustable="box")
        fig.colorbar(err_im, cax=cax_err, label=r"log$_{10}$ relative error")
        fig.savefig(out / f"xy{plane}_full_vs_c8.png", dpi=180)
        plt.close(fig)

        # The y=0 transverse cut uses exactly the same grid points for both methods.
        center = np.isclose(pf["y_m"], 0.0)
        xcut = pf["x_m"][center] * 1e3
        order = np.argsort(xcut)
        fig, axes = plt.subplots(2, 1, figsize=(7, 5.2), sharex=True)
        axes[0].plot(xcut[order], pf["E_mag"][center][order], label="full")
        axes[0].plot(xcut[order], pm["E_mag"][center][order], "--", label="C8 modal")
        axes[0].set(ylabel="|E| [V/m]", title=f"XY plane {plane}: y=0 transverse cut")
        axes[0].legend()
        axes[1].semilogy(xcut[order], np.maximum(error[center][order], 1e-16))
        axes[1].axhline(1e-10, color="k", ls="--", lw=.8)
        axes[1].set(xlabel="x [mm]", ylabel="relative |E| error")
        axes[1].grid(True, which="both", alpha=.3)
        fig.tight_layout()
        fig.savefig(out / f"xy{plane}_transverse_cut.png", dpi=180)
        plt.close(fig)

    info = metrics(result / "metrics.txt")
    labels = ["full operator", "C8 blocks"]
    memory = [float(info["full_operator_mib"]), float(info["c8_blocks_mib"])]
    build = [float(info["full_operator_build_s"]), float(info["c8_blocks_build_s"])]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    axes[0].bar(labels, memory); axes[0].set_ylabel("stored operator [MiB]")
    axes[1].bar(labels, build); axes[1].set_ylabel("build time [s]")
    fig.tight_layout(); fig.savefig(out / "cost_comparison.png", dpi=180); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=Path("results/validation"))
    args = parser.parse_args()
    create_figures(args.result)


if __name__ == "__main__":
    main()
