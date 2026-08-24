# Validation report

## Environment status

- Source inspected at fixed commit: `d37311da872406c1a439f1c184f12d80fb8269d7`
- Modified Fortran sources: 9 files parsed successfully with `fparser 0.2.4`
- Native Fortran build: not run because this container has no `gfortran`, `ifx`, `ifort`, `flang`, `lfortran`, or `nvfortran`
- Numerically equivalent Python reference: executed successfully

No native Fortran timing or compiler result is fabricated. Run `make debug` on a machine with GNU Fortran before using the target model.

## Lightweight numerical reference

The repeated reference run used the same panel geometry, constant-vector panel kernel, visibility gate, C8 local-coordinate rotation, DFT signs, current reconstruction, and field postprocessor as the Fortran implementation.

| Item | Result |
|---|---:|
| Frequency | 3 GHz |
| Triangular panels | 480 |
| Vector degrees of freedom | 1440 |
| Full stored complex elements | 2,073,600 |
| C8 stored complex elements | 259,200 |
| Full storage | 31.640625 MiB |
| C8 storage | 3.955078 MiB |
| Full operator build | 0.5737 s |
| C8 block build | 0.1297 s |
| Block-circulant error | 9.63e-16 |
| Cumulative full/all-mode current error | 6.29e-16 |
| Cumulative full/m=1,7 current error | 6.32e-16 |
| XZ electric-field error | 7.13e-16 |
| Maximum XY electric-field error | 6.17e-16 |
| Axis electric-field error | 4.73e-16 |
| Exit electric-field error | 4.86e-16 |
| Full/modal exit power | 2.1436793e-6 / 2.1436793e-6 W |

The excitation mode norms were `m=1: 7.38545e-2`, `m=7: 7.38532e-2`; every other mode was below `1.0e-17`. The maximum reflection order was reached before the `1e-4` stopping condition, so this run validates Full/C8 algebra and common field reconstruction, not converged IPO physics.

Reference arrays, metrics, and figures are in `results/python_reference_updated`.

## 94 GHz preflight

With approximately lambda/6 panels, the input namelist uses 323,180 panels per sector and 2,585,440 panels overall. Estimated dense storage is about 896,461 GiB for Full and 112,058 GiB for the eight C8 blocks. Dense execution is therefore intentionally rejected.

`run_94ghz_sweep` writes 135 frequency/polarization/offset/angle conditions. At a fixed frequency the propagation operator can be reused; changing frequency requires rebuilding it. A practical target solver must replace dense `K_d` storage/application with GPU matrix-free evaluation, FMM/MLFMA, or admissible-block H-matrix/ACA.
