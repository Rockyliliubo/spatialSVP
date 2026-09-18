test_that("Westfall-Young adjusted p-values dominate raw permutation p-values", {
  sim <- simulate_svp_data(n_spots = 300, n_genes = 80, n_sets = 10,
                           n_spatial_sets = 3, set_size = 8, seed = 11)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  fdr <- permutation_fdr(act, method = "morans_i", n_perm = 99, seed = 11)
  expect_s3_class(fdr, "svp_results")
  expect_true(all(fdr$p_adj_wy >= fdr$p_perm - 1e-12, na.rm = TRUE))
  expect_true(all(fdr$p_perm >= 1 / 100))
  expect_true(all(fdr$p_perm <= 1))
})

test_that("Westfall-Young controls family-wise error under the global null", {
  sim <- simulate_svp_data(n_spots = 300, n_genes = 160, n_sets = 20,
                           n_spatial_sets = 0, set_size = 8, seed = 12)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  fdr <- permutation_fdr(act, method = "morans_i", n_perm = 99, seed = 12)
  expect_lte(sum(fdr$p_adj_wy < 0.05, na.rm = TRUE), 1L)
  expect_lte(sum(fdr$p_adj_BH < 0.05, na.rm = TRUE), 4L)
})

test_that("planted sets are recovered by Westfall-Young correction", {
  sim <- simulate_svp_data(n_spots = 400, n_genes = 80, n_sets = 10,
                           n_spatial_sets = 4, set_size = 8, seed = 13)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z", preprocess = "none")
  fdr <- permutation_fdr(act, method = "morans_i", n_perm = 99, seed = 13)
  spn <- sim$truth$set[sim$truth$spatial]
  nln <- sim$truth$set[!sim$truth$spatial]
  sp <- fdr$p_perm[fdr$pathway %in% spn]
  nl <- fdr$p_perm[fdr$pathway %in% nln]
  expect_lt(mean(sp), mean(nl))
  expect_gte(sum(fdr$p_adj_wy[fdr$pathway %in% spn] < 0.1), 2L)
})

test_that("blocked Westfall-Young runs and is reproducible", {
  sim <- simulate_svp_data(n_spots = 400, n_genes = 80, n_sets = 8,
                           n_spatial_sets = 2, set_size = 8,
                           n_blocks = 2, seed = 14)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, block = sim$block,
                                  method = "mean_z", preprocess = "none")
  f1 <- permutation_fdr(act, method = "morans_i", n_perm = 49, seed = 14)
  f2 <- permutation_fdr(act, method = "morans_i", n_perm = 49, seed = 14)
  expect_equal(f1$p_perm, f2$p_perm)
  expect_false(anyNA(f1$p_perm))
  expect_error(permutation_fdr(act, method = "morans_i", n_perm = 5),
               "n_perm")
})
