# C8 modal PO/IPO implementation

The complete build, formulation, validation results, limitations, DFT signs, and 94 GHz scaling discussion are in [README_full_vs_c8.md](README_full_vs_c8.md).

Quick start:

```bash
make debug
./build/compare_full_vs_c8 config_validation.nml
python3 plot_compare.py --result results/validation
./build/run_94ghz_sweep config_94ghz.nml
```

The Full dense and Full matrix-free paths use the same 3x3 vector panel-pair map. The C8 path constructs only the eight relative-sector blocks, solves all eight DFT modes, and also evaluates the `m=1,7` fast path. Automatic mode selection is available as `mode_policy=2` with `mode_tolerance`.

The 94 GHz wavelength-resolved case is a preflight by design: C8 symmetry alone does not make the dense operator feasible. Use the reported application points for GPU, FMM/MLFMA, or H-matrix/ACA acceleration.
