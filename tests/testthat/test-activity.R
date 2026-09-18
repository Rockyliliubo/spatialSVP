test_that("mean_z scoring matches a hand-computed expectation", {
  set.seed(1)
  counts <- matrix(c(
    1, 1, 1, 10,
    1, 1, 1, 10,
    10, 10, 10, 1,
    5, 5, 5, 5,
    2, 8, 2, 8,
    8, 2, 8, 2
  ), nrow = 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  sets <- data.frame(source = c("A", "A", "B", "B"),
                     target = c("g1", "g2", "g3", "g4"),
                     mor = c(1, 1, 1, 1))
  coords <- data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2))
  act <- compute_pathway_activity(counts, sets, coords = coords,
                                  method = "mean_z", preprocess = "none",
                                  min_size = 2L)
  expect_s3_class(act, "svp_activity")
  expect_equal(colnames(act$activity), c("A", "B"))
  expect_equal(rownames(act$activity), colnames(counts))
  # z-scores per gene; constant genes get z = 0 by package convention
  z1 <- as.numeric(scale(c(1, 1, 1, 10)))
  z3 <- as.numeric(scale(c(10, 10, 10, 1)))
  z4 <- rep(0, 4)
  expect_equal(unname(act$activity[, "A"]), (z1 + z1) / 2, tolerance = 1e-8)
  expect_equal(unname(act$activity[, "B"]), (z3 + z4) / 2, tolerance = 1e-8)
  # g1 and g2 identical -> mean equals z1
  expect_equal(unname(act$activity[, "A"]), z1, tolerance = 1e-8)
})

test_that("network filtering drops small/unknown sets and aggregates dups", {
  counts <- matrix(rpois(200, 5), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:20)))
  coords <- data.frame(x = rnorm(20), y = rnorm(20))
  sets <- data.frame(
    source = c(rep("OK", 6), rep("small", 2), "ghost", "OK"),
    target = c(paste0("g", 1:6), paste0("g", 1:2), "g_absent", "g1"),
    mor = c(rep(1, 8), 1, 1))
  act <- compute_pathway_activity(counts, sets, coords = coords,
                                  method = "mean_z", min_size = 4L)
  expect_equal(colnames(act$activity), "OK")
  # duplicated (OK, g1) row must have been aggregated, not double-counted
  idx <- which(act$network$source == "OK" & act$network$target == "g1")
  expect_length(idx, 1)
  expect_false("ghost" %in% act$network$target)
})

test_that("list and matrix gene-set formats are accepted", {
  counts <- matrix(rpois(200, 5), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("s", 1:20)))
  coords <- data.frame(x = rnorm(20), y = rnorm(20))
  a1 <- compute_pathway_activity(counts, list(S = paste0("g", 1:6)),
                                 coords = coords, min_size = 5L)
  wm <- matrix(0, nrow = 10, ncol = 1,
               dimnames = list(paste0("g", 1:10), "S"))
  wm[paste0("g", 1:6), "S"] <- 1
  a2 <- compute_pathway_activity(counts, wm, coords = coords, min_size = 5L)
  expect_equal(a1$activity, a2$activity)
})

test_that("decoupler back end runs and separates planted sets", {
  skip_if_not_installed("decoupleR")
  sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords,
                                  method = "decoupler_mlm")
  expect_s3_class(act, "svp_activity")
  expect_false(anyNA(act$activity))
  # spatial sets should vary more across spots than null sets
  vv <- apply(act$activity, 2, var)
  spn <- sim$truth$set[sim$truth$spatial]
  nln <- sim$truth$set[!sim$truth$spatial]
  expect_true(mean(vv[colnames(act$activity) %in% spn]) >
                mean(vv[colnames(act$activity) %in% nln]))
})

test_that("gsva back end runs when GSVA is available", {
  skip_if_not_installed("GSVA")
  sim <- simulate_svp_data(n_spots = 100, n_genes = 80, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "gsva")
  expect_false(anyNA(act$activity))
})

test_that("fgsea back end matches fgseaSimple and separates planted sets", {
  skip_if_not_installed("fgsea")
  sim <- simulate_svp_data(n_spots = 150, n_genes = 100, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "fgsea",
                                  min_size = 5L)
  expect_s3_class(act, "svp_activity")
  expect_false(anyNA(act$activity))
  # independent reference: fgseaSimple on one spot's z-scores
  z <- spatialSVP:::.zscore_rows(spatialSVP:::.logcpm(sim$counts))
  stats <- z[, 3]
  names(stats) <- rownames(z)
  sets <- split(spatialSVP:::.filter_network(
    spatialSVP:::.as_network(sim$gene_sets), rownames(sim$counts),
    min_size = 5L)$target,
    spatialSVP:::.filter_network(
      spatialSVP:::.as_network(sim$gene_sets), rownames(sim$counts),
      min_size = 5L)$source)
  ref <- fgsea::fgseaSimple(sets, stats, nperm = 10, minSize = 5)
  got <- act$activity[3, names(sets)]
  expect_equal(unname(got), ref$ES, tolerance = 1e-8)
  # planted spatial sets vary more across spots than null sets
  vv <- apply(act$activity, 2, stats::var)
  spn <- sim$truth$set[sim$truth$spatial]
  nln <- sim$truth$set[!sim$truth$spatial]
  expect_true(mean(vv[colnames(act$activity) %in% spn]) >
                mean(vv[colnames(act$activity) %in% nln]))
})

test_that("get_progeny_sets validates arguments (offline-safe)", {
  expect_error(get_progeny_sets("zebrafish"), "one of")
})

test_that("block is validated and stored", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 60, n_sets = 6,
                           n_spatial_sets = 2, set_size = 8, seed = 1)
  blk <- rep(c("p1", "p2"), each = 50)
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, block = blk)
  expect_s3_class(act$block, "factor")
  expect_error(compute_pathway_activity(sim$counts, sim$gene_sets,
                                        coords = sim$coords,
                                        block = blk[1:10]),
               "match")
})
