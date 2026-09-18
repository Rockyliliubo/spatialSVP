# spatialSVP 0.3.0

* Matched reference sets preserve the original signed scoring weights within
  expression bins. Unit-weight gene sets retain identical sampled members
  under the same seed. Weighted footprint analyses should be rerun.
* Documentation distinguishes absolute association, competitive enrichment,
  and spatial organization after expression-component adjustment. Principal
  component removal can suppress coherent biological signal and is a
  sensitivity analysis rather than a generally preferred preprocessing step.
* Analysis scripts include direct method comparisons, residualization
  sensitivity, independent-donor controls and memory-bounded Moran scoring
  verified against the full-matrix implementation.

# spatialSVP 0.2.1

* `matched_null()` now fails with an informative error when an expression
  bin cannot support exact same-size matching; it no longer returns smaller
  matched sets as a fallback.

* `matched_enrichment()` gains `stat_matched_mad` and `z_matched` columns
  (robust standardized effect). The ratio `stat / stat_matched_med` can flip
  sign or explode because the reference median of a centered statistic can
  sit at or below zero; use `z_matched` or rank-based effects for
  quantitative comparisons, ratios for display only.

* `.coords_xy()` keeps spot names when it rebuilds the coordinate data
  frame, so activity and coordinate inputs stay aligned after subsetting.

* `Depends:` now declares `R (>= 4.3)` (was `R (>= 4.2.0)`).

# spatialSVP 0.2.0

## Behaviour changes

* `test_spatial_variability(method = "gp")` now defaults to permutation
  p-values (`n_perm = 199`), like the other back ends. The analytic
  mixture-chi-square p-value is uncalibrated for the maximum over the
  lengthscale x rho grid and is available only via explicit `n_perm = 0`,
  with a warning. gp results report `lengthscale = NA` when the
  intercept-only null model wins.
* `residualize` in `compute_pathway_activity()` is now a three-tier choice
  `"none" | "libsize" | "pcs"` (was logical). Logical values keep working
  with a deprecation warning (`TRUE -> "libsize"`). New `n_pcs` argument
  for the `"pcs"` tier.
* Moran's I statistics carry the conventional `n/(n-1)` factor (textbook
  scale; p-values unchanged).
* `mean_z` normalizes by `sum(abs(mor))` so weighted resources (e.g.
  PROGENy) keep their scale; identical output when `mor = 1`.
* covariance back end: per-scale statistics are studentized with
  permutation-pooled null moments, so `dominant_scale` recovers the true
  signal scale instead of collapsing onto the largest kernel.
  `dominant_scale` is a coarse descriptor; use gp `lengthscale` for sharp
  scale estimates.
* `svp_input()` uses the Seurat v5 `layer` argument (silences the `slot`
  deprecation warning).

## New features

* `matched_null()` + `matched_enrichment()`: the RELATIVE null — is a
  pathway more spatially organized than size- and expression-matched random
  gene sets? On structured tissue the absolute null can be rejected for many
  ordinary aggregates, so the matched test answers a distinct competitive
  question.
  Statistics only; cost scales with the number of gene sets, not with
  permutation counts.
* `get_progeny_sets()`: PROGENy perturbation-response footprints fetched
  at runtime via decoupleR (resources are never redistributed).
* New scoring back end `method = "fgsea"`: per-spot fgsea enrichment
  scores (rank-based; ignores `mor` weights; no permutations).
* `simulate_svp_data()` gains `coord_scale` for pixel-scale coordinate
  simulations.
