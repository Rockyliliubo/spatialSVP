# Calibration tests for the gp back end. The length-scale search does not
# follow a single-boundary mixture, so gp p-values use shared permutations
# by default. The analytic p-value remains available through n_perm = 0 as
# a diagnostic with a warning. These tests use Visium-like pixel-scale
# coordinates (276 px spacing).

pixel_coords <- function(n_side = 20, spacing = 276, seed = 1) {
  g <- expand.grid(x = seq_len(n_side), y = seq_len(n_side))
  set.seed(seed)
  data.frame(x = g$x * spacing + stats::runif(nrow(g), -20, 20),
             y = g$y * spacing + stats::runif(nrow(g), -20, 20))
}

noise_activity <- function(n, p, seed) {
  set.seed(seed)
  A <- matrix(stats::rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("N", seq_len(p))))
  A
}

test_that("gp permutation p is valid (conservative, never liberal) under the null", {
  # A variance-component LRT carries a point mass at 0 under the null, so
  # the tie-conservative permutation p piles mass at exactly 1.0: the test
  # is conservative by construction and must never be liberal. (A uniform
  # KS test is the wrong target here: the null p distribution is discrete
  # with mass at 1, and multi-pathway calls share permutations by design.)
  coords <- pixel_coords(n_side = 15)  # 225 spots keep the loop fast
  ps <- vapply(seq_len(40), function(s) {
    A <- noise_activity(nrow(coords), 1, seed = 100 + s)
    test_spatial_variability(A, coords, method = "gp", n_perm = 99,
                             seed = s)$p_value
  }, numeric(1))
  expect_gt(mean(ps == 1), 0.2)     # boundary mass present (conservative)
  expect_lte(mean(ps < 0.05), 0.15) # type-I under control
  expect_lte(mean(ps < 0.01), 0.06)
})

test_that("gp analytic p warns; rho = 0 restores boundary-null behaviour", {
  coords <- pixel_coords()
  A <- noise_activity(nrow(coords), 60, seed = 43)
  expect_warning(
    res <- test_spatial_variability(A, coords, method = "gp", n_perm = 0),
    "uncalibrated")
  # Boundary theory: under the null a variance component pinned at the
  # boundary has P(LRT = 0) ~ 0.5, so the exact p = 0.5 mass should sit
  # near one half; the accepted range below allows finite-sample variation.
  mass <- mean(abs(res$p_value - 0.5) < 1e-12)
  expect_gt(mass, 0.25)
  expect_lt(mass, 0.72)
  # the alternative nests the null: the LRT is non-negative by construction
  expect_true(all(res$stat >= 0))
})

test_that("gp enforces the n_perm >= 19 guard like the other back ends", {
  coords <- pixel_coords()
  A <- noise_activity(nrow(coords), 5, seed = 44)
  expect_error(test_spatial_variability(A, coords, method = "gp", n_perm = 5),
               "19")
})

test_that("gp default (n_perm = NULL) uses permutations, not the analytic p", {
  coords <- pixel_coords()
  A <- noise_activity(nrow(coords), 5, seed = 45)
  res <- test_spatial_variability(A, coords, method = "gp", n_perm = NULL,
                                  seed = 1)
  # permutation p-values with n_perm = 199 live on the 1/200 grid and can
  # never equal the analytic mixture's exact 0.5 mass point at lrt = 0
  expect_true(all(res$p_value %in% ((1:200) / 200)))
})
