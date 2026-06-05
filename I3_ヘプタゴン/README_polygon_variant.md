# heptagonal_cone

This is a 7-sided polygonal-cone variant generated from the uploaded octagonal closed-cap concentrator code.

Main changes:

- `mod_geometry.f90`: soft polygon radius uses `N_sides = 7` instead of 8.
- `mod_config.f90`: output root is `results_heptagonal`.
- `mod_config.f90`: IPO checkpoint directory is `ipo_checkpoint_closed_cap_heptagonal` to avoid accidentally reusing octagonal checkpoints.
- `Makefile`: executable name is `bin/closed_cap_heptagonal` and visualization root is `results_heptagonal`.
- `visualize.py` and `visualize_result_bins_with_wall.py`: theoretical wall shape and labels are updated for heptagon geometry.
- `main_closed_cap.f90`: one duplicate `end do` in `build_split_xy_plane_list` was removed so the code can compile.

Use:

```bash
make clean
make
make run
make visualize_bins
```

If you change geometry-related parameters, run:

```bash
make clean_checkpoints
```

Otherwise humanity will compare a heptagon wall to old surface-current data and pretend physics is at fault.
