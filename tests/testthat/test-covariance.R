# Tests for the covariance back end's cross-scale standardization: the raw
# max-over-scales statistic was decided by kernel geometry (dominant_scale
# pinned at one grid point whatever the implanted scale); after per-scale
# standardization the dominant scale must track the true signal scale.

test_that(".stat_covmax handles a single pathway (orientation regression)", {
  set.seed(1)
  coords <- data.frame(x = rep(1:10, each = 10) * 276,
                       y = rep(1:10, 10) * 276)
  Z <- matrix(scale(rnorm(100)), ncol = 1)
  oc <- spatialSVP:::.stat_covmax(Z, spatialSVP:::.cov_setup(coords))
  expect_length(oc$stat, 1L)
  expect_length(oc$scale, 1L)
  expect_true(is.finite(oc$stat))
})

test_that("matched_enrichment covariance needs >= 5 scored sets", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 2,
                           n_spatial_sets = 1, set_size = 8, seed = 1)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 1,
                     n_bins = 2L, min_size = 4L, seed = 1)
  act <- compute_pathway_activity(sim$counts, mn$network,
                                  coords = sim$coords, method = "mean_z",
                                  preprocess = "none", min_size = 4L)
  # 2 originals + 2 matched replicates = 4 columns < 5
  expect_error(matched_enrichment(act, mn, method = "covariance"),
               "5")
})

test_that("dominant_scale responds to the implanted scale (not a constant)", {
  # Visium-like grid: 25 x 25 spots at 276 px spacing
  gx <- rep(1:25, each = 25)
  gy <- rep(1:25, 25)
  coords <- data.frame(x = gx * 276, y = gy * 276)
  ctr <- c(median(coords$x), median(coords$y))
  bump <- function(sigma) {
    d2 <- (coords$x - ctr[1])^2 + (coords$y - ctr[2])^2
    exp(-d2 / (2 * sigma^2))
  }
  set.seed(2)
  A <- cbind(bump(400), bump(1300), bump(3000),
             matrix(rnorm(625 * 3), nrow = 625))
  colnames(A) <- c("sig400", "sig1300", "sig3000", "n1", "n2", "n3")
  res <- test_spatial_variability(A, coords, method = "covariance",
                                  n_perm = 29, seed = 1)
  ds <- res$dominant_scale[match(c("sig400", "sig1300", "sig3000"),
                                 res$pathway)]
  # regression case: without studentization every implanted scale reported the SAME
  # grid point. Now: the extreme scales separate and the order is monotone
  # (the statistic is a coarse descriptor with ~2 grid points resolution;
  # the sharp scale estimator is the gp lengthscale)
  expect_gt(length(unique(ds)), 1L)
  expect_lt(ds[1], ds[3])
  expect_true(all(diff(ds) >= 0))
  # and sits inside the physical grid range
  nn <- 276
  span <- max(diff(range(coords$x)), diff(range(coords$y)))
  expect_true(all(ds >= nn & ds <= span))
})
