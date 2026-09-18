# Tests for matched_null() / matched_enrichment(): the
# gene-set-level relative null. The reference class is defined by set size,
# expression-bin composition and signed score weights.

test_that("matched_null reproduces size and expression-bin composition", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 200, n_sets = 6,
                           n_spatial_sets = 2, set_size = 10, seed = 1)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 5, n_bins = 5,
                     min_size = 5, seed = 2)
  expect_s3_class(mn, "svp_matched")
  expect_output(print(mn), "svp_matched")
  # matched replicate sizes equal the parent size
  tab <- table(mn$network$source)
  parents <- split(as.character(mn$network$target),
                   as.character(mn$network$source))
  real <- split(as.character(sim$gene_sets$target), sim$gene_sets$source)
  for (src in names(mn$parents)) {
    expect_equal(length(parents[[src]]), length(real[[mn$parents[[src]]]]))
  }
  # bin composition matches the parent one-to-one
  bins <- as.integer(cut(rank(mn$expr_stat), breaks = 5,
                         include.lowest = TRUE))
  names(bins) <- names(mn$expr_stat)
  for (src in head(names(mn$parents), 3)) {
    expect_equal(unname(sort(table(bins[parents[[src]]]))),
                 unname(sort(table(bins[real[[mn$parents[[src]]]]]))))
  }
  # matched members never overlap their own parent set
  for (src in names(mn$parents)) {
    expect_length(intersect(parents[[src]], real[[mn$parents[[src]]]]), 0)
  }
  # with_original: network starts with the originals
  expect_true(all(sim$gene_sets$target %in% mn$network$target))
})

test_that("matched_null is reproducible under a seed", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 160, n_sets = 4,
                           n_spatial_sets = 1, set_size = 10, seed = 3)
  m1 <- matched_null(sim$gene_sets, sim$counts, n_match = 4, seed = 7)
  m2 <- matched_null(sim$gene_sets, sim$counts, n_match = 4, seed = 7)
  expect_identical(m1$network, m2$network)
})

test_that("matched_null preserves signed weights within expression bins", {
  set.seed(301)
  expr <- matrix(rpois(400 * 80, lambda = 8), 400, 80,
    dimnames = list(paste0("g", seq_len(400)), paste0("s", seq_len(80))))
  net <- data.frame(source = "weighted", target = rownames(expr)[seq_len(30)],
                    mor = rep(c(-3, -0.5, 0.2, 1, 4), 6))
  weighted <- matched_null(net, expr, n_match = 7, n_bins = 5, seed = 44)
  plain <- matched_null(transform(net, mor = 1), expr,
                         n_match = 7, n_bins = 5, seed = 44)
  expect_identical(weighted$network$target, plain$network$target)
  expect_identical(weighted$network$source, plain$network$source)
  bins <- setNames(as.integer(cut(rank(weighted$expr_stat), breaks = 5,
    include.lowest = TRUE)), names(weighted$expr_stat))
  source_pairs <- sort(paste(bins[net$target], net$mor, sep = ":"))
  for (id in names(weighted$parents)) {
    d <- weighted$network[weighted$network$source == id, ]
    expect_identical(sort(paste(bins[d$target], d$mor, sep = ":")), source_pairs)
    expect_false(anyDuplicated(d$target) > 0)
  }
  expect_true(all(plain$network$mor == 1))
})

test_that("matched_null fails when exact size matching is impossible", {
  expr <- matrix(seq_len(80), nrow = 8,
                 dimnames = list(paste0("g", 1:8), paste0("s", 1:10)))
  sets <- data.frame(source = "A", target = paste0("g", 1:5), mor = 1)
  expect_error(
    matched_null(sets, expr, preprocess = "none", n_match = 2,
                 n_bins = 2, min_size = 1, seed = 1),
    "Cannot construct exact matched sets"
  )
})

test_that("matched_enrichment separates planted spatial sets from nulls", {
  sim <- simulate_svp_data(n_spots = 400, n_genes = 200, n_sets = 8,
                           n_spatial_sets = 3, set_size = 10, effect = 2,
                           seed = 11)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 20, seed = 1)
  act <- compute_pathway_activity(sim$counts, mn$network,
                                  coords = sim$coords, method = "mean_z",
                                  preprocess = "none")
  enr <- matched_enrichment(act, mn)
  expect_s3_class(enr, "svp_results")
  expect_true(all(enr$n_matched == 20L))
  expect_true(all(enr$p_matched >= 0 & enr$p_matched <= 1, na.rm = TRUE))
  sp <- sim$truth$set[sim$truth$spatial]
  nl <- sim$truth$set[!sim$truth$spatial]
  expect_lt(median(enr$p_matched[enr$pathway %in% sp]),
            median(enr$p_matched[enr$pathway %in% nl]))
  # matched controls should sit far below the observed spatial stat
  expect_gt(median(enr$stat[enr$pathway %in% sp]),
            median(enr$stat_matched_med[enr$pathway %in% sp]))
  expect_output(print(enr), "svp_results")
})

test_that("matched_enrichment: null sets get non-significant matched p", {
  # global-null simulation: no set should look enriched vs its own matched
  # replicates beyond chance (loose bound; the relative null is calibrated
  # by construction because replicates are exchangeable with the parent
  # under the null)
  sim <- simulate_svp_data(n_spots = 300, n_genes = 200, n_sets = 6,
                           n_spatial_sets = 0, set_size = 10, seed = 12)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 19, seed = 1)
  act <- compute_pathway_activity(sim$counts, mn$network,
                                  coords = sim$coords, method = "mean_z",
                                  preprocess = "none")
  enr <- matched_enrichment(act, mn)
  expect_lte(sum(enr$p_matched < 0.05, na.rm = TRUE), 2L)
})

test_that("matched_enrichment validates its inputs", {
  sim <- simulate_svp_data(n_spots = 100, n_genes = 100, n_sets = 4,
                           n_spatial_sets = 1, set_size = 8, seed = 5)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 3, seed = 1)
  # activity scored WITHOUT the matched columns must error
  act <- compute_pathway_activity(sim$counts, sim$gene_sets,
                                  coords = sim$coords, method = "mean_z")
  expect_error(matched_enrichment(act, mn), "matched")
  expect_error(matched_enrichment(act, list()), "svp_matched")
})

test_that("matched_enrichment works with covariance and gp statistics", {
  sim <- simulate_svp_data(n_spots = 225, n_genes = 160, n_sets = 6,
                           n_spatial_sets = 2, set_size = 10, effect = 2,
                           seed = 21)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 9, seed = 1)
  act <- compute_pathway_activity(sim$counts, mn$network,
                                  coords = sim$coords, method = "mean_z",
                                  preprocess = "none")
  e1 <- matched_enrichment(act, mn, method = "covariance")
  e2 <- matched_enrichment(act, mn, method = "gp")
  expect_false(anyNA(e1$p_matched))
  expect_false(anyNA(e2$p_matched))
  sp <- sim$truth$set[sim$truth$spatial]
  expect_lte(median(rank(e1$p_matched)[e1$pathway %in% sp]), 3)
  expect_lte(median(rank(e2$p_matched)[e2$pathway %in% sp]), 3)
})

test_that("matched_enrichment reports MAD-based robust effect (z_matched)", {
  sim <- simulate_svp_data(n_spots = 300, n_genes = 200, n_sets = 6,
                           n_spatial_sets = 2, set_size = 10, seed = 13)
  mn <- matched_null(sim$gene_sets, sim$counts, n_match = 19, seed = 1)
  act <- compute_pathway_activity(sim$counts, mn$network,
                                  coords = sim$coords, method = "mean_z",
                                  preprocess = "none")
  enr <- matched_enrichment(act, mn)
  expect_true(all(c("stat_matched_mad", "z_matched") %in% colnames(enr)))
  # z_matched is by construction (stat - med) / mad wherever mad > 0
  has_mad <- is.finite(enr$stat_matched_mad) & enr$stat_matched_mad > 0
  expect_true(any(has_mad))
  expect_equal(enr$z_matched[has_mad],
               (enr$stat[has_mad] - enr$stat_matched_med[has_mad]) /
                 enr$stat_matched_mad[has_mad],
               tolerance = 1e-10)
  # degenerate reference class (mad == 0) must yield NA, never Inf/0
  expect_true(all(is.na(enr$z_matched[!has_mad])))
  # planted spatial sets should sit at positive z (above their reference)
  sp <- sim$truth$set[sim$truth$spatial]
  expect_gt(stats::median(enr$z_matched[enr$pathway %in% sp], na.rm = TRUE),
            0)
})
