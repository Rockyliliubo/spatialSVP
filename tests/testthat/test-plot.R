sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
                         n_spatial_sets = 2, set_size = 8, seed = 21)
act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                coords = sim$coords, method = "mean_z")

test_that("plot_pathway_spatial returns a ggplot", {
  p <- plot_pathway_spatial(act, "SPATIAL_1")
  expect_s3_class(p, "ggplot")
  expect_error(plot_pathway_spatial(act, "NOPE"), "not found")
  expect_silent(print(p))
})

test_that("plot_volcano works for both result types", {
  res <- test_spatial_variability(act, method = "morans_i",
                                  n_perm = 39, seed = 21)
  p1 <- plot_volcano(res)
  expect_s3_class(p1, "ggplot")
  expect_silent(print(p1))

  fdr <- permutation_fdr(act, method = "morans_i", n_perm = 49, seed = 21)
  p2 <- plot_volcano(fdr, p_col = "p_perm", fdr_col = "p_adj_wy")
  expect_s3_class(p2, "ggplot")
  expect_silent(print(p2))

  expect_error(plot_volcano(res, p_col = "nope"), "missing")
})

test_that("plot_lengthscale_map uses gp lengthscale or covariance scale", {
  res_gp <- test_spatial_variability(act, method = "gp", seed = 21)
  p1 <- plot_lengthscale_map(res_gp)
  expect_s3_class(p1, "ggplot")
  expect_silent(print(p1))

  res_cov <- test_spatial_variability(act, method = "covariance",
                                      n_perm = 39, seed = 21)
  p2 <- plot_lengthscale_map(res_cov)
  expect_s3_class(p2, "ggplot")

  res_mi <- test_spatial_variability(act, method = "morans_i",
                                     n_perm = 39, seed = 21)
  expect_error(plot_lengthscale_map(res_mi), "length-scale")
})

test_that("printing results does not recurse", {
  res <- test_spatial_variability(act, method = "morans_i",
                                  n_perm = 39, seed = 22)
  out <- capture.output(print(res))
  expect_gt(length(out), 3L)
  sub <- res[order(res$p_value), ]
  out2 <- capture.output(print(sub))
  expect_gt(length(out2), 3L)
})
