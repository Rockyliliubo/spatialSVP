# Equivariance and scale-calibration tests for the spatial back ends:
# Multiplying coordinates by s > 0 must multiply lengthscale estimates by s
# and leave statistics/p-values unchanged; morans_i (rank-based adjacency)
# must be fully invariant. Pixel-scale fixtures pin the grids to Visium-like
# coordinate magnitudes.

test_that("gp lengthscale scales with coordinates; stats and p-values do not", {
  sim <- simulate_svp_data(n_spots = 150, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 5)
  s <- 250
  a1 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords, method = "mean_z")
  a2 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords * s, method = "mean_z")
  r1 <- test_spatial_variability(a1, method = "gp", seed = 1)
  r2 <- test_spatial_variability(a2, method = "gp", seed = 1)
  expect_equal(r2$lengthscale, r1$lengthscale * s, tolerance = 1e-6)
  expect_equal(r2$stat, r1$stat, tolerance = 1e-6)
  expect_equal(r2$p_value, r1$p_value, tolerance = 1e-6)
})

test_that("covariance dominant_scale scales with coordinates", {
  sim <- simulate_svp_data(n_spots = 150, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 5)
  s <- 400
  a1 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords, method = "mean_z")
  a2 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords * s, method = "mean_z")
  r1 <- test_spatial_variability(a1, method = "covariance", n_perm = 19,
                                 seed = 1)
  r2 <- test_spatial_variability(a2, method = "covariance", n_perm = 19,
                                 seed = 1)
  expect_equal(r2$dominant_scale, r1$dominant_scale * s, tolerance = 1e-6)
  expect_equal(r2$stat, r1$stat, tolerance = 1e-6)
})

test_that("morans_i is fully invariant under coordinate rescaling", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 9)
  s <- 276
  a1 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords, method = "mean_z")
  a2 <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                 coords = sim$coords * s, method = "mean_z")
  r1 <- test_spatial_variability(a1, method = "morans_i", n_perm = 19,
                                 seed = 1)
  r2 <- test_spatial_variability(a2, method = "morans_i", n_perm = 19,
                                 seed = 1)
  expect_equal(r2$stat, r1$stat, tolerance = 1e-8)
  expect_equal(r2$p_value, r1$p_value, tolerance = 1e-8)
})

test_that("gp separates planted sets on Visium-like pixel-scale coordinates", {
  sim <- simulate_svp_data(n_spots = 225, n_genes = 120, n_sets = 8,
                           n_spatial_sets = 3, set_size = 8,
                           pattern = "bump", effect = 2, seed = 21)
  s <- 250
  coords <- sim$coords * s
  act <- compute_pathway_activity(sim$counts, sim$gene_sets, coords = coords,
                                  method = "mean_z")
  res <- test_spatial_variability(act, method = "gp", seed = 1)
  sp <- sim$truth$set[sim$truth$spatial]
  nul <- sim$truth$set[!sim$truth$spatial]
  expect_lt(median(res$p_value[res$pathway %in% sp]),
            median(res$p_value[res$pathway %in% nul]))
  # lengthscale estimates must sit inside the physical coordinate range:
  # at least the spot spacing, at most the tissue span
  d2 <- as.matrix(dist(coords))^2
  diag(d2) <- Inf
  nn <- sqrt(median(apply(d2, 1, min)))
  span <- max(diff(range(coords$x)), diff(range(coords$y)))
  expect_true(all(res$lengthscale >= 0.9 * nn, na.rm = TRUE))
  expect_true(all(res$lengthscale <= span, na.rm = TRUE))
})
