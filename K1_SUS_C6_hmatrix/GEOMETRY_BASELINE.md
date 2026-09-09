# K1_SUS_C6_geometry_opt

Independent implementation of geometry optimizations 2 and 3, based on
K1_SUS_C6. The original directory, executable and checkpoints were not edited.
170 GHz, all six C6 modes, SUS surface impedance and equivalent magnetic
current are retained. No empirical damping, distance cutoff or beam cutoff.

## Changes

1. Visibility: replace ten interior samples by a piecewise-convex test.
   The taper and pipe are individually convex. A wall-centroid segment
   within either piece is inside it; a crossing segment is visible exactly
   when its junction intersection lies in the hexagon. Each side-plane
   inequality is affine along each segment piece, so endpoints and junction
   suffice. The original 2e-12 m geometric tolerance is retained.
   visible_pair requires valid mesh wall centroids; it is not an arbitrary
   endpoint containment API. Straight and expanding tapers are supported.
2. Candidate exclusion: omit the same-sector, same-piece coplanar block
   before the C6 source loop. Such pairs are grazing and have zero response
   under the existing IPO front/back rule. Cross-piece pairs remain.
   A conservative unnormalized sign rejection precedes sqrt/division;
   survivors still pass the original normalized 1e-12 angular test.

Candidate exclusion introduces no additional physical approximation.
Analytic visibility corrects missed occlusions: it is NOT numerically
equivalent to the old ten-sample approximation. Use a new output directory.
Asymptotic complexity remains O(Q^2).

## Validation performed (2026-09-09)

From this directory run: bash test_geometry.sh

The script builds and runs small C6 tests (no full 3D solver calls):
- Contracting, straight and expanding meshes, Q=42 per sector.
  All 63,504 geometry candidate pairs per mesh are checked individually.
  These are pair-map checks, not a full 3D propagation solve.
- Exact visibility against 101 samples including endpoints, plus a
  constructed narrow junction obstruction missed by old sampling.
- Candidate exclusion versus an unculled analytic reference: maximum map
  difference and relative six-mode current difference are both zero.
  Input excites all six modes and all three vector components.
- Seven Python state-management tests, TE/TM surface-response tests,
  and continuous versus interrupted/resumed C6 runs through order 3.

On the contracting Q=42 mesh, old sampling misses 324 directed occlusions;
180 pair maps change (others already fail front/back tests).
See geometry_test.log and geometry_benchmark.json.

Q=384, two OpenMP threads, three repetitions, median one-step worker time
including mesh setup, field probes, subprocess and interchange:
old 0.352 s, optimized 0.209 s, approximately 1.68x in this small test.
This is NOT a production speed prediction or routine-level profile.
Changed occlusion also changes work performed.
One-step relative differences versus the old sampled model:
current 1.145%, scattered probe electric field 1.195%.
These differences do not establish physical accuracy of either result.
Q=16 results are recorded but are dominated by startup overhead.

## Python and restart

Use ../.venv/bin/python; the Makefile defaults to that interpreter.
Original crash-safe per-order saving and restart logic is retained.
Original checkpoints must not be reused: source/executable identity and
visibility operator have changed. The copied production configuration
is present for inspection; no production mesh computation was run.
Heavy runs require explicit user instruction.

Preflight distinguishes candidates from calls after coplanar exclusion.
Storage estimates are unchanged. With production mesh settings, exclusion
removes Q_cone^2 + Q_pipe^2 candidates, approximately 8.39% of 6Q^2.
Further front/back and visibility rejection occurs inside pair_map_local.

## Limits and provenance

Physical/order/mesh convergence and power balance remain unverified.
Original beam waist position, SUS constants, centroid quadrature,
self-term treatment and PO assumptions are unchanged.
See reference/README_original.md for the original assumptions.
That archived README contains historical test-status and Python setup
instructions; use this README for the independent copy.

reference/operator.txt and reference/modal.txt preserve original operators.
Renamed modules isolate old sampling and analytic visibility without
culling for regression. source_manifest.json records original file hashes.
