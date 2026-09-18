sim_small <- function(seed = 1, spatial = 2, sets = 6) {
  simulate_svp_data(n_spots = 400, n_genes = 80, n_sets = sets,
                    n_spatial_sets = spatial, set_size = 8, seed = seed)
}
# name-based indexing: pathway column order must not be assumed
sp_sets <- function(sim) sim$truth$set[sim$truth$spatial]
nl_sets <- function(sim) sim$truth$set[!sim$truth$spatial]

test_that("morans_i detects planted spatial sets and calibrates nulls", {
  sim <- sim_small(seed = 1)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  res <- test_spatial_variability(act, method = "morans_i",
                                  n_perm = 99, seed = 1)
  expect_s3_class(res, "svp_results")
  sp <- res$p_value[res$pathway %in% sp_sets(sim)]
  nl <- res$p_value[res$pathway %in% nl_sets(sim)]
  expect_lt(mean(sp), mean(nl))
  expect_gte(sum(res$p_adj_BH[res$pathway %in% sp_sets(sim)] < 0.1), 2L)
})

test_that("permutation p-values are calibrated under the global null", {
  sim <- sim_small(seed = 2, spatial = 0)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  res <- test_spatial_variability(act, method = "morans_i",
                                  n_perm = 99, seed = 2)
  p <- res$p_value
  expect_lte(mean(p < 0.05), 0.25)
  ks <- suppressWarnings(stats::ks.test(p, "punif"))
  expect_gt(ks$p.value, 0.01)
})

test_that("covariance back end detects planted sets and reports scales", {
  sim <- sim_small(seed = 3)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  res <- test_spatial_variability(act, method = "covariance",
                                  n_perm = 99, seed = 3)
  sp <- res$p_value[res$pathway %in% sp_sets(sim)]
  nl <- res$p_value[res$pathway %in% nl_sets(sim)]
  expect_lt(mean(sp), mean(nl))
  expect_true(all(is.finite(res$dominant_scale)))
  expect_true(all(res$dominant_scale > 0))
})

test_that("gp back end flags spatial sets and gives finite lengthscales", {
  sim <- sim_small(seed = 4)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  res <- test_spatial_variability(act, method = "gp", n_perm = 99, seed = 4)
  sp_rank <- rank(res$p_value)[res$pathway %in% sp_sets(sim)]
  expect_lte(median(sp_rank), 3)
  # planted (strong-signal) sets must beat the null model and thus carry a
  # finite lengthscale; null sets may select rho = 0 and report NA
  sp_ls <- res$lengthscale[res$pathway %in% sp_sets(sim)]
  expect_true(all(is.finite(sp_ls)))
  expect_true(all(sp_ls > 0))
  expect_true(all(res$spatial_var_fraction >= 0 &
                    res$spatial_var_fraction < 1))
  # permutation p-values are in [0, 1]
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
})

test_that("gp analytic ordering agrees with an independent permutation back end", {
  sim <- sim_small(seed = 5, spatial = 3)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  # n_perm = 0 is the diagnostic analytic path and warns
  expect_warning(
    ra <- test_spatial_variability(act, method = "gp", n_perm = 0, seed = 5),
    "uncalibrated")
  rm <- test_spatial_variability(act, method = "morans_i",
                                 n_perm = 49, seed = 5)
  m <- match(ra$pathway, rm$pathway)
  ct <- suppressWarnings(stats::cor.test(ra$p_value, rm$p_value[m],
                                         method = "spearman",
                                         exact = FALSE))
  expect_gt(unname(ct$estimate), 0.5)
})

test_that("blocked analyses run and respect block structure", {
  sim <- simulate_svp_data(n_spots = 400, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8,
                           n_blocks = 2, seed = 6)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords,
                                  block = sim$block, method = "mean_z")
  res <- test_spatial_variability(act, method = "morans_i",
                                  n_perm = 49, seed = 6)
  expect_false(anyNA(res$p_value))
  # significance survives within-block testing for at least one planted set
  expect_gte(sum(res$p_adj_BH[res$pathway %in% sp_sets(sim)] < 0.1), 1L)
})

test_that("zero-variance pathways are reported as NA, not dropped silently", {
  A <- cbind(rep(1, 101), matrix(rnorm(202), 101, 2))
  coords <- data.frame(x = rnorm(101), y = rnorm(101))
  res <- test_spatial_variability(A, coords, method = "morans_i",
                                  n_perm = 29, seed = 1)
  expect_true(is.na(res$stat[1]))
  expect_false(anyNA(res$stat[2:3]))
  expect_error(test_spatial_variability(cbind(rep(1, 100), rep(2, 100)),
                                        data.frame(x = 1:100, y = 1:100),
                                        method = "morans_i", n_perm = 29),
               "variance")
})
