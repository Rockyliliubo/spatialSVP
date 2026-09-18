test_that("simulate_svp_data produces consistent structures", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  expect_s3_class(sim$coords, "data.frame")
  expect_equal(nrow(sim$coords), 100)
  expect_equal(nrow(sim$counts), 60)
  expect_equal(ncol(sim$counts), 100)
  expect_true(all(c("source", "target", "mor") %in% colnames(sim$gene_sets)))
  expect_equal(nrow(sim$truth), 6)
  expect_equal(sum(sim$truth$spatial), 2)
  expect_null(sim$block)

  # null genes dedicated to null sets: sets must be gene-disjoint
  genes <- split(sim$gene_sets$target, sim$gene_sets$source)
  expect_equal(length(unique(unlist(genes))),
               sum(lengths(genes)))

  # reproducibility
  sim2 <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                            n_spatial_sets = 2, set_size = 8, seed = 1)
  expect_equal(sim$counts, sim2$counts)
  expect_equal(sim$coords, sim2$coords)
})

test_that("simulate_svp_data handles blocks and edge cases", {
  sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 8,
                           n_spatial_sets = 2, set_size = 8,
                           n_blocks = 2, pattern = "mix", seed = 2)
  expect_s3_class(sim$block, "factor")
  expect_equal(length(sim$block), 200)
  expect_length(levels(sim$block), 2)
  expect_true(all(sim$truth$pattern[sim$truth$spatial] %in%
                    c("bump", "gradient")))

  expect_error(simulate_svp_data(n_genes = 10, n_sets = 6, set_size = 8),
               "n_genes")
  expect_error(simulate_svp_data(n_sets = 4, n_spatial_sets = 5),
               "n_spatial_sets")
})

test_that("coord_scale rescales coordinates and truth centers", {
  s1 <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                          n_spatial_sets = 2, set_size = 8, seed = 4)
  s2 <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                          n_spatial_sets = 2, set_size = 8, coord_scale = 276,
                          seed = 4)
  expect_equal(s2$coords, s1$coords * 276)
  expect_equal(s2$truth$center_x, s1$truth$center_x * 276)
  expect_equal(s2$counts, s1$counts)  # signal is defined on the unit grid
})

test_that("the tests are one-sided: checkerboard anticorrelation is not flagged", {
  # the package tests positive spatial
  # autocorrelation. A checkerboard (negative autocorrelation) must not be
  # flagged; note the symmetrized k = 6 kNN graph mixes in same-colour
  # diagonal neighbours, which dilutes the negative statistic towards zero
  # (a pure 4-neighbour grid would give Moran's I ~ -0.29)
  g <- expand.grid(x = 1:15, y = 1:15)
  coords <- data.frame(x = g$x, y = g$y)
  checker <- ((g$x + g$y) %% 2) * 2 - 1
  smooth <- as.numeric(scale(g$x))
  set.seed(1)
  A <- cbind(checker + rnorm(225, 0, 0.01), smooth + rnorm(225, 0, 0.01),
             matrix(rnorm(225 * 4), ncol = 4))
  colnames(A) <- c("checker", "smooth", paste0("n", 1:4))
  res <- test_spatial_variability(A, coords, method = "morans_i",
                                  n_perm = 49, seed = 1)
  expect_gt(res$p_value[res$pathway == "checker"], 0.05)
  expect_lt(res$p_value[res$pathway == "smooth"], 0.05)
  expect_gt(res$stat[res$pathway == "smooth"],
            res$stat[res$pathway == "checker"])
})
