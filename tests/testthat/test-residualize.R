test_that("residualization matches lm() residuals on log10 library size", {
  set.seed(7)
  base <- matrix(rpois(10 * 30, 8), nrow = 10, ncol = 30,
                 dimnames = list(paste0("g", 1:10), paste0("s", 1:30)))
  scale <- rep(c(0.5, 1, 4), each = 10)
  counts <- sweep(base, 2L, scale, "*")
  storage.mode(counts) <- "integer"
  lc <- spatialSVP:::.logcpm(counts)
  got <- spatialSVP:::.resid_libsize(lc, counts)
  x <- log10(colSums(counts))
  for (g in paste0("g", 1:10)) {
    fit <- stats::lm(lc[g, ] ~ x)
    expect_equal(unname(got[g, ]), unname(stats::residuals(fit)),
                 tolerance = 1e-8)
  }
})

test_that("constant library size makes residualization a no-op after z-scoring", {
  vals <- c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89)
  set.seed(3)
  counts <- vapply(seq_len(20), function(...) vals[sample.int(10)],
                   numeric(10))
  rownames(counts) <- paste0("g", 1:10)
  colnames(counts) <- paste0("s", 1:20)
  coords <- data.frame(x = stats::rnorm(20), y = stats::rnorm(20))
  sets <- data.frame(source = "S", target = paste0("g", 1:6), mor = 1)
  a0 <- compute_pathway_activity(counts, sets, coords = coords,
                                 min_size = 5L, residualize = "none")
  a1 <- compute_pathway_activity(counts, sets, coords = coords,
                                 min_size = 5L, residualize = "libsize")
  expect_equal(a0$activity, a1$activity)
})

test_that("residualize is rejected for user-normalized input", {
  counts <- matrix(rpois(100, 5), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:10)))
  coords <- data.frame(x = stats::rnorm(10), y = stats::rnorm(10))
  sets <- data.frame(source = "S", target = paste0("g", 1:6), mor = 1)
  expect_error(compute_pathway_activity(counts, sets, coords = coords,
                                        preprocess = "none",
                                        residualize = "libsize"),
               "logcpm")
  expect_error(compute_pathway_activity(counts, sets, coords = coords,
                                        preprocess = "none",
                                        residualize = "pcs"),
               "logcpm")
})

test_that("legacy logical residualize maps onto the three-tier API", {
  counts <- matrix(rpois(200, 5), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:20)))
  coords <- data.frame(x = stats::rnorm(20), y = stats::rnorm(20))
  sets <- data.frame(source = "S", target = paste0("g", 1:6), mor = 1)
  expect_warning(
    a_true <- compute_pathway_activity(counts, sets, coords = coords,
                                       min_size = 5L, residualize = TRUE),
    "deprecated")
  a_lib <- compute_pathway_activity(counts, sets, coords = coords,
                                    min_size = 5L, residualize = "libsize")
  expect_equal(a_true$activity, a_lib$activity)
  expect_equal(a_lib$residualize, "libsize")
  a_false <- compute_pathway_activity(counts, sets, coords = coords,
                                      min_size = 5L, residualize = FALSE)
  expect_equal(a_false$residualize, "none")
})

test_that("library residualization breaks compositional coupling", {
  sim <- simulate_svp_data(n_spots = 225, n_genes = 140, n_sets = 10,
                           n_spatial_sets = 2, set_size = 10,
                           effect = 2.5, seed = 11)
  null_sets <- sim$truth$set[!sim$truth$spatial]
  spatial_sets <- sim$truth$set[sim$truth$spatial]

  run <- function(resid) {
    act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                    coords = sim$coords,
                                    method = "mean_z",
                                    residualize = resid)
    permutation_fdr(act, method = "morans_i", n_perm = 99, seed = 1)
  }
  f0 <- run("none")
  f1 <- run("libsize")

  fp0 <- mean(f0$p_adj_wy[f0$pathway %in% null_sets] < 0.05)
  fp1 <- mean(f1$p_adj_wy[f1$pathway %in% null_sets] < 0.05)
  # without residualization the composition effect fires
  expect_gt(fp0, 0.5)
  # residualization removes most of the inherited structure
  expect_lt(fp1, 0.2)
  # planted spatial sets remain detected after residualization
  expect_true(all(f1$p_adj_wy[f1$pathway %in% spatial_sets] < 0.05))
})

test_that("residualization integrates with the decoupler back end", {
  skip_if_not_installed("decoupleR")
  sim <- simulate_svp_data(n_spots = 100, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords,
                                  method = "decoupler_mlm",
                                  residualize = "libsize")
  expect_false(anyNA(act$activity))
  expect_equal(act$residualize, "libsize")
  expect_output(print(act), "resid = libsize")
})

test_that("pcs residualization equals an explicit PC projection", {
  set.seed(17)
  counts <- matrix(rpois(20 * 40, 8), nrow = 20,
                   dimnames = list(paste0("g", 1:20), paste0("s", 1:40)))
  lc <- spatialSVP:::.logcpm(counts)
  got <- spatialSVP:::.resid_pcs(lc, n_pcs = 4)
  # independent path: SVD principal scores + per-gene linear regression
  X <- t(lc - rowMeans(lc))                # spots x genes, gene-centered
  sv <- svd(X, nu = 4, nv = 0)
  sc <- sv$u
  for (g in paste0("g", c(1, 7, 20))) {
    fit <- stats::lm(lc[g, ] ~ sc - 1)     # scores are mean-orthogonal
    expect_equal(unname(got[g, ]), unname(lc[g, ] - stats::fitted(fit)),
                 tolerance = 1e-8)
  }
  # gene means are preserved
  expect_equal(rowMeans(got), rowMeans(lc), tolerance = 1e-10)
})

test_that("pcs residualization removes a shared spatial factor", {
  set.seed(23)
  n <- 100
  coords <- data.frame(x = stats::runif(n), y = stats::runif(n))
  grad <- coords$x - mean(coords$x)        # a strong shared spatial axis
  logexpr <- matrix(stats::rnorm(60 * n, 5, 0.8), nrow = 60)
  logexpr <- logexpr + outer(stats::runif(60, 0.5, 1.5), grad * 3)
  counts <- matrix(stats::rpois(60 * n, exp(logexpr)), nrow = 60,
                   dimnames = list(paste0("g", 1:60), paste0("s", 1:n)))
  lc <- spatialSVP:::.logcpm(counts)
  m1 <- spatialSVP:::.resid_libsize(lc, counts)
  m2 <- spatialSVP:::.resid_pcs(m1, n_pcs = 3)
  cor_after <- cor(t(m2), grad)[, 1]
  expect_lt(mean(abs(cor_after)), 0.05)
})

test_that("pcs tier runs end to end, caps n_pcs, and prints", {
  sim <- simulate_svp_data(n_spots = 120, n_genes = 100, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 31)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z",
                                  residualize = "pcs", n_pcs = 5)
  expect_false(anyNA(act$activity))
  expect_equal(act$residualize, "pcs")
  expect_output(print(act), "resid = pcs")
  # n_pcs larger than the spot count is capped, not an error
  act_big <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                      coords = sim$coords, method = "mean_z",
                                      residualize = "pcs", n_pcs = 500L)
  expect_false(anyNA(act_big$activity))
})

test_that(".pc_scores RSpectra path gives an orthonormal top-k subspace", {
  testthat::skip_if_not_installed("RSpectra")
  set.seed(41)
  n <- 8100  # past the full-eigen cutoff
  g <- 40
  mc <- matrix(stats::rnorm(g * n), nrow = g)
  # plant two strong latent axes so the top subspace is well separated
  # (signal must dominate the noise bulk ~100x for >99% subspace overlap)
  a1 <- stats::rnorm(n)
  a2 <- stats::rnorm(n)
  mc <- mc + outer(stats::rnorm(g, 30, 5), a1) + outer(stats::rnorm(g, 20, 3), a2)
  mc <- mc - rowMeans(mc)
  sc <- spatialSVP:::.pc_scores(mc, 5L)
  expect_equal(dim(sc), c(n, 5L))
  # orthonormal
  expect_equal(crossprod(sc), diag(5), tolerance = 1e-6)
  # the planted axes lie inside the returned subspace (projection recovers
  # essentially all of their variance)
  rec <- sc %*% crossprod(sc, cbind(a1 / sd(a1), a2 / sd(a2)))
  expect_gt(cor(rec[, 1], a1 / sd(a1))^2, 0.99)
  expect_gt(cor(rec[, 2], a2 / sd(a2))^2, 0.99)
})

test_that("coords keep spot names through compute_pathway_activity", {
  set.seed(31)
  counts <- matrix(rpois(200 * 60, 2), nrow = 200, ncol = 60,
                   dimnames = list(paste0("G", 1:200), paste0("BC", 1:60)))
  coords <- data.frame(x = runif(60), y = runif(60))
  rownames(coords) <- paste0("BC", 1:60)
  gs <- data.frame(source = "S1", target = paste0("G", 1:20), mor = 1)
  act <- compute_pathway_activity(counts, gs, coords = coords,
                                  method = "mean_z", min_size = 5L)
  expect_identical(rownames(act$coords), paste0("BC", 1:60))
  expect_identical(rownames(act$activity), paste0("BC", 1:60))
})
