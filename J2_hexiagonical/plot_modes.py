#!/usr/bin/env python3
"""Visualize the physical meaning and induced exit field of every C6 mode."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


M = 6
MODE_INFO = {
    0: ("common / C6 symmetric", "0 deg per sector"),
    1: ("first-order +", "+60 deg per sector"),
    2: ("second-order +", "+120 deg per sector"),
    3: ("alternating", "180 deg per sector"),
    4: ("second-order -", "-120 deg per sector"),
    5: ("first-order -", "-60 deg per sector"),
}


def load_csv(path: Path) -> np.ndarray:
    return np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))


def pure_mode_weights(mode: int) -> np.ndarray:
    sectors = np.arange(M)
    return np.exp(1j * 2.0 * np.pi * mode * sectors / M)


def draw_hex_mode(ax: plt.Axes, mode: int, cmap: mpl.colors.Colormap) -> None:
    weights = pure_mode_weights(mode)
    vertices = np.column_stack(
        (np.cos(np.arange(M) * 2.0 * np.pi / M - np.pi / M),
         np.sin(np.arange(M) * 2.0 * np.pi / M - np.pi / M))
    )
    vertices = np.vstack((vertices, vertices[0]))
    for sector in range(M):
        phase = np.angle(weights[sector])
        color = cmap((phase + np.pi) / (2.0 * np.pi))
        p0, p1 = vertices[sector], vertices[sector + 1]
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=12,
                solid_capstyle="butt")
        center = 0.5 * (p0 + p1)
        ax.text(1.17 * center[0], 1.17 * center[1],
                f"p={sector}\n{np.degrees(phase):+.0f}°",
                ha="center", va="center", fontsize=7)
    ax.plot(vertices[:, 0], vertices[:, 1], color="0.2", lw=0.7)
    ax.set_aspect("equal")
    ax.set_xlim(-1.42, 1.42)
    ax.set_ylim(-1.32, 1.32)
    ax.axis("off")
    label, step = MODE_INFO[mode]
    ax.set_title(f"m={mode}: {label}\n{step}", fontsize=10)


def plot_mode_meaning(output: Path) -> Path:
    cmap = mpl.colormaps["twilight"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2))
    for mode, ax in enumerate(axes.flat):
        draw_hex_mode(ax, mode, cmap)
    sm = mpl.cm.ScalarMappable(
        norm=mpl.colors.Normalize(vmin=-180.0, vmax=180.0), cmap=cmap
    )
    sm.set_array([])
    fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.028, pad=0.035,
                 label="local-current phase [deg]")
    fig.suptitle(
        r"Pure $C_6$ current modes: "
        r"$J_p^{(m)}=J_0^{(m)}e^{i2\pi mp/6}$",
        fontsize=14,
    )
    output.mkdir(parents=True, exist_ok=True)
    path = output / "c6_mode_physical_meaning.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def reshape_plane(data: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.unique(data["x_m"])
    y = np.unique(data["y_m"])
    return x, y, values.reshape(len(y), len(x))


def field_phasors(data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ex = data["Ex_re"] + 1j * data["Ex_im"]
    ey = data["Ey_re"] + 1j * data["Ey_im"]
    ez = data["Ez_re"] + 1j * data["Ez_im"]
    emag = np.sqrt(np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2)
    return ex, ey, ez, emag


def phase_aligned_transverse_vectors(
    ex: np.ndarray, ey: np.ndarray, ez: np.ndarray, emag: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    peak = int(np.nanargmax(emag))
    reference_vector = np.array([ex[peak], ey[peak], ez[peak]])
    dominant = reference_vector[int(np.argmax(np.abs(reference_vector)))]
    rotation = np.exp(-1j * np.angle(dominant)) if abs(dominant) > 0.0 else 1.0
    ux = np.real(ex * rotation)
    uy = np.real(ey * rotation)
    length = np.hypot(ux, uy)
    good = length > max(float(np.nanmax(length)) * 1.0e-6, 1.0e-300)
    ux = np.divide(ux, length, out=np.zeros_like(ux), where=good)
    uy = np.divide(uy, length, out=np.zeros_like(uy), where=good)
    return ux, uy


def load_summary(result: Path) -> dict[int, np.void]:
    path = result / "mode_field_summary.csv"
    if not path.exists():
        return {}
    data = load_csv(path)
    return {int(row["mode"]): row for row in data}


def plot_induced_fields(result: Path, output: Path) -> Path | None:
    mode_paths = [result / f"mode_m{mode}_exit.csv" for mode in range(M)]
    available = [path.exists() for path in mode_paths]
    if not any(available):
        return None
    if not all(available):
        missing = ", ".join(path.name for path, ok in zip(mode_paths, available) if not ok)
        raise FileNotFoundError(f"Incomplete C6 mode-field output; missing: {missing}")

    datasets = [load_csv(path) for path in mode_paths]
    phasors = [field_phasors(data) for data in datasets]
    summary = load_summary(result)
    global_reference = max(float(np.nanmax(parts[3])) for parts in phasors)
    global_reference = max(global_reference, 1.0e-300)

    fig, axes = plt.subplots(2, 3, figsize=(12, 8.2), constrained_layout=True)
    image = None
    for mode, (ax, data, parts) in enumerate(zip(axes.flat, datasets, phasors)):
        ex, ey, ez, emag = parts
        mask = data["mask"] > 0
        display = np.where(
            mask,
            20.0 * np.log10(np.maximum(emag / global_reference, 1.0e-6)),
            np.nan,
        )
        x, y, image_values = reshape_plane(data, display)
        image = ax.imshow(
            image_values,
            origin="lower",
            extent=[x[0] * 1e3, x[-1] * 1e3, y[0] * 1e3, y[-1] * 1e3],
            vmin=-60.0,
            vmax=0.0,
            cmap="magma",
        )

        mode_peak = float(np.nanmax(np.where(mask, emag, np.nan)))
        if mode_peak > global_reference * 1.0e-12:
            ux, uy = phase_aligned_transverse_vectors(ex, ey, ez, emag)
            _, _, ux_grid = reshape_plane(data, np.where(mask, ux, np.nan))
            _, _, uy_grid = reshape_plane(data, np.where(mask, uy, np.nan))
            stride = max(1, len(x) // 11)
            xx, yy = np.meshgrid(x * 1e3, y * 1e3)
            ax.quiver(
                xx[::stride, ::stride],
                yy[::stride, ::stride],
                ux_grid[::stride, ::stride],
                uy_grid[::stride, ::stride],
                color="cyan",
                pivot="mid",
                scale=18,
                width=0.004,
                headwidth=3,
            )
        else:
            ax.text(
                0.5, 0.5, "not excited\n(zero within tolerance)",
                transform=ax.transAxes, ha="center", va="center",
                color="white", fontsize=10,
                bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none"},
            )

        norm_text = ""
        if mode in summary:
            norm_value = float(summary[mode]["cumulative_dft_current_norm"])
            norm_text = f"\nDFT current norm={norm_value:.3e}"
        ax.set_title(f"m={mode}: {MODE_INFO[mode][0]}{norm_text}", fontsize=9)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_aspect("equal", adjustable="box")

    assert image is not None
    fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.028, pad=0.025,
                 label=r"induced $|E_m|$ relative to all-mode maximum [dB]")
    fig.suptitle(
        "Exit-plane electric field induced by each cumulative C6 surface-current mode\n"
        "(incident field excluded; cyan arrows: phase-aligned transverse direction)",
        fontsize=13,
    )
    output.mkdir(parents=True, exist_ok=True)
    path = output / "c6_mode_induced_exit_fields.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot C6 mode phase meaning and mode-resolved induced electric fields."
    )
    parser.add_argument("--result", type=Path, default=Path("results/validation_split"))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output if args.output is not None else args.result / "figures"

    meaning_path = plot_mode_meaning(output)
    field_path = plot_induced_fields(args.result, output)
    print(f"wrote {meaning_path}")
    if field_path is None:
        print(
            "mode field CSV files were not found; run the modal solver with "
            "write_mode_fields=.true. to create the induced-field figure"
        )
    else:
        print(f"wrote {field_path}")


if __name__ == "__main__":
    main()
