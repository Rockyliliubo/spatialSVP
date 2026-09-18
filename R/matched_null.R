#' Matched random gene sets: the gene-set-level relative null
#'
#' On structured tissue the absolute null "this score vector has no spatial
#' structure" can be false for many ordinary gene sets. Cell density, library
#' composition, regional programs and cell-type geography can make random
#' gene aggregates spatially associated even when they do not encode a named
#' biological process. The competitive question is therefore distinct: is a
#' pathway MORE spatially variable than a random gene set of the same size
#' and the same marginal expression distribution?
#'
#' \code{matched_null()} draws, for every input gene set, \code{n_match}
#' random sets of the same size whose members match the real members
#' one-to-one within expression bins (rank-based bins of mean log-CPM,
#' robust to the zero-inflated tie mass of sparse matrices). Real members
#' are excluded from their own sampling pools. Matched sets preserve the
#' parent's signed weights within expression bins; unweighted sets retain
#' unit weights. They are returned in the same tidy format as the input,
#' ready to be scored by
#' \code{\link{compute_pathway_activity}} alongside the real sets.
#' Matching is fail-closed: if any expression bin lacks enough eligible genes
#' to reproduce a parent set exactly, the function stops with an informative
#' error instead of returning smaller matched sets.
#'
#' \code{matched_enrichment()} then computes, per real set, an empirical
#' one-sided p-value against the spatial statistics of its own matched
#' replicates -- significance means "more spatially variable than matched
#' controls". BH-adjusted reference p-values are also returned. Their
#' inferential interpretation depends on the tested set being exchangeable
#' with its references under the competitive null; size and expression
#' matching alone do not ensure this for correlated biological gene sets.
#' Only the spatial statistic is computed (no coordinate permutations),
#' so the cost scales with the number of sets, not with permutation counts.
#'
#' @param gene_sets Gene sets in any format accepted by
#'   \code{\link{compute_pathway_activity}} (tidy data.frame, named list, or
#'   genes x sets weight matrix). Matched replicates preserve the parent's
#'   \code{mor} values within expression bins.
#' @param expr Genes x spots expression matrix used for the matching
#'   statistic (and typically the same matrix that will be scored).
#' @param preprocess \code{"logcpm"} (default): match on row means of
#'   log-CPM; \code{"none"}: match on row means of \code{expr} as supplied
#'   (already normalized).
#' @param n_match Number of matched random replicates per gene set
#'   (default 50; the smallest achievable empirical p is 1/(n_match+1).
#'   Roughly 20x P matches is sufficient for a single floor p-value to cross
#'   a 5 percent BH threshold across P sets; fewer can suffice when several
#'   sets share small p-values. Small n_match is mainly useful for screening
#'   and ranking).
#' @param n_bins Number of rank-based expression bins (default 10).
#' @param min_size Minimum number of set genes present in \code{expr} for a
#'   set to enter matching.
#' @param with_original If \code{TRUE} (default), the returned
#'   \code{network} already contains the original sets (unchanged sources)
#'   plus the matched replicates, so it can be scored as-is.
#' @param seed Optional random seed.
#'
#' @return \code{matched_null()} returns a \code{svp_matched} object: a list
#'   with \code{network} (tidy data.frame for scoring), \code{parents}
#'   (named character vector mapping each matched source to its parent set),
#'   \code{originals} (the retained original set ids), \code{expr_stat}
#'   (the per-gene matching statistic), \code{n_match} and \code{n_bins}.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 120, n_sets = 8,
#'                          n_spatial_sets = 3, set_size = 8, seed = 1)
#' mn <- matched_null(sim$gene_sets, sim$counts, n_match = 20, seed = 1)
#' mn
#' act <- compute_pathway_activity(sim$counts, mn$network,
#'                                 coords = sim$coords, method = "mean_z")
#' enr <- matched_enrichment(act, mn)
#' head(enr[order(enr$p_matched), ])
#'
#' @export
matched_null <- function(gene_sets, expr, preprocess = c("logcpm", "none"),
                         n_match = 50L, n_bins = 10L, min_size = 5L,
                         with_original = TRUE, seed = NULL) {
  preprocess <- match.arg(preprocess)
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) {
    stop("'expr' must have gene ids as rownames.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  net <- .as_network(gene_sets)
  net <- .filter_network(net, rownames(expr), min_size = min_size)
  if (nrow(net) == 0L) {
    stop("No gene set retains at least ", min_size,
         " genes present in 'expr'.", call. = FALSE)
  }

  mexpr <- if (preprocess == "logcpm") rowMeans(.logcpm(expr)) else {
    rowMeans(expr)
  }
  # rank-based bins: cut() on raw quantile breaks collapses on the large
  # zero-tie mass of sparse matrices, and cut(labels = FALSE) drops names;
  # rank() gives balanced, tie-robust bins.
  bins <- as.integer(cut(rank(mexpr), breaks = n_bins, include.lowest = TRUE))
  bins <- stats::setNames(bins, rownames(expr))

  members <- split(as.character(net$target), as.character(net$source))
  originals <- names(members)
  # pools per expression bin computed once (the inner setdiff against the
  # parent members stays linear per replicate)
  pools <- lapply(seq_len(n_bins), function(b) rownames(expr)[bins == b])
  # list-accumulate one chunk per replicate, then a single unlist: growing
  # vectors with c() is O(n^2) and costs ~30 min at n_match = 1990 scale
  chunk_src <- vector("list", length(originals) * n_match)
  chunk_tgt <- vector("list", length(originals) * n_match)
  chunk_par <- vector("list", length(originals) * n_match)
  chunk_mor <- vector("list", length(originals) * n_match)
  n_chunks <- 0L
  for (s in originals) {
    g <- members[[s]]
    tab <- table(bins[g])
    parent_net <- net[as.character(net$source) == s, , drop = FALSE]
    source_weights <- parent_net$mor[match(g, parent_net$target)]
    # Sampling visits expression bins in names(tab) order. Carry each source
    # weight to the corresponding sampled slot; this preserves sign and
    # weight magnitude without changing gene draws or consuming extra RNG.
    matched_weights <- unlist(lapply(names(tab), function(lv)
      source_weights[bins[g] == as.integer(lv)]), use.names = FALSE)
    for (lv in names(tab)) {
      pool <- setdiff(pools[[as.integer(lv)]], g)
      k <- unname(tab[[lv]])
      if (length(pool) < k) {
        stop("Cannot construct exact matched sets for '", s,
             "': expression bin ", lv, " requires ", k,
             " eligible genes but only ", length(pool), " are available. ",
             "Reduce 'n_bins', supply a larger expression matrix, or remove ",
             "the affected gene set.", call. = FALSE)
      }
    }
    for (r in seq_len(n_match)) {
      pick <- unlist(lapply(names(tab), function(lv) {
        pool <- setdiff(pools[[as.integer(lv)]], g)
        k <- tab[[lv]]
        sample(pool, k, replace = FALSE)
      }))
      if (length(pick) == 0L) next
      n_chunks <- n_chunks + 1L
      chunk_src[[n_chunks]] <- rep(paste0("SVP_MATCH__", s, "__", r),
                                   length(pick))
      chunk_tgt[[n_chunks]] <- pick
      chunk_par[[n_chunks]] <- s
      chunk_mor[[n_chunks]] <- as.numeric(matched_weights)
    }
  }
  if (n_chunks == 0L) {
    stop("No matched replicates could be drawn (pools exhausted?).",
         call. = FALSE)
  }
  take <- seq_len(n_chunks)
  matched_sources <- unlist(chunk_src[take], use.names = FALSE)
  matched_targets <- unlist(chunk_tgt[take], use.names = FALSE)
  parents <- stats::setNames(unlist(chunk_par[take], use.names = FALSE),
                      vapply(take, function(i) chunk_src[[i]][1L], ""))

  matched_df <- data.frame(source = matched_sources,
                           target = matched_targets,
                           mor = unlist(chunk_mor[take], use.names = FALSE),
                           stringsAsFactors = FALSE)
  network <- if (with_original) {
    orig_df <- data.frame(source = as.character(net$source),
                          target = net$target, mor = net$mor,
                          stringsAsFactors = FALSE)
    rbind(orig_df, matched_df)
  } else {
    matched_df
  }
  out <- list(network = network, parents = parents, originals = originals,
              expr_stat = mexpr, n_match = as.integer(n_match),
              n_bins = as.integer(n_bins))
  class(out) <- "svp_matched"
  out
}

#' @export
print.svp_matched <- function(x, ...) {
  cat("svp_matched: ", length(x$originals), " gene sets x ", x$n_match,
      " matched replicates (", length(unique(x$parents)), " parents with ",
      "replicates; ", x$n_bins, " expression bins)\n", sep = "")
  invisible(x)
}

#' Matched-null enrichment test
#'
#' @param activity A \code{svp_activity} object (or spots x pathways matrix)
#'   whose columns include BOTH the original sets and the matched replicates
#'   generated by \code{\link{matched_null}} (score \code{mn$network}).
#' @param matched The \code{svp_matched} object from
#'   \code{\link{matched_null}}.
#' @param method Spatial statistic back end (see
#'   \code{\link{test_spatial_variability}}); only the statistic is
#'   computed, no coordinate permutations.
#' @param coords Spot coordinates; defaults to the object's stored coords.
#' @param block Optional block factor per spot.
#'
#' @return A data.frame of class \code{svp_results} with one row per
#'   original gene set and columns \code{pathway}, \code{stat} (observed
#'   spatial statistic), \code{stat_matched_med} (median statistic of the
#'   matched replicates), \code{stat_matched_mad} (their MAD around that
#'   median), \code{z_matched} (robust standardized effect:
#'   \code{(stat - stat_matched_med) / stat_matched_mad}; preferred over the
#'   plain ratio, whose denominator can sit at or below zero and flip sign),
#'   \code{n_matched}, \code{p_matched} (one-sided
#'   empirical p-value: the fraction of matched replicates at least as
#'   spatially variable as the observed set), and \code{p_adj_BH} (BH across
#'   the original sets). See \code{\link{matched_null}} for the semantics.
#'
#' @export
matched_enrichment <- function(activity, matched,
                               method = c("morans_i", "covariance", "gp"),
                               coords = NULL, block = NULL) {
  if (!inherits(matched, "svp_matched")) {
    stop("'matched' must be a svp_matched object from matched_null().",
         call. = FALSE)
  }
  method <- match.arg(method)
  ex <- .extract_activity(activity, coords, block)
  A <- ex$A
  coords <- ex$coords
  block <- ex$block
  cn <- colnames(A) %||% paste0("pathway", seq_len(ncol(A)))
  colnames(A) <- cn

  have <- matched$originals[matched$originals %in% cn]
  missing <- setdiff(matched$originals, cn)
  if (length(missing) > 0L) {
    warning("Original sets not found among the activity columns (dropped ",
            "during scoring?): ", paste(missing, collapse = ", "),
            call. = FALSE)
  }
  kids <- matched$parents[matched$parents %in% have]
  kids <- kids[names(kids) %in% cn]
  if (length(kids) == 0L) {
    stop("No matched replicate columns found in 'activity'. Score ",
         "matched$network with compute_pathway_activity() first.",
         call. = FALSE)
  }

  sd_col <- apply(A, 2L, stats::sd)
  ok <- is.finite(sd_col) & sd_col > 0
  Z <- scale(A[, ok, drop = FALSE])
  if (identical(method, "covariance") && ncol(Z) < 5L) {
    stop("method 'covariance' standardizes each kernel scale across ",
         "pathways and needs >= 5 scored sets.",
         call. = FALSE)
  }
  mask <- .same_block_mask(block)
  stat <- switch(method,
    morans_i = {
      W <- .adjacency_w(coords, k = 6L, mask = mask)
      S0 <- sum(W)
      if (S0 == 0) {
        stop("No spot pairs within blocks to build adjacency.",
             call. = FALSE)
      }
      .stat_morans(Z, W, S0)
    },
    covariance = .stat_covmax(Z, .cov_setup(coords, block))$stat,
    gp = .gp_stats(Z, .gp_setup(coords, block))$stat)
  stat_full <- rep(NA_real_, ncol(A))
  names(stat_full) <- cn
  stat_full[ok] <- stat

  res <- data.frame(pathway = have, stat = unname(stat_full[have]),
                    stringsAsFactors = FALSE)
  by_parent <- split(names(kids), kids)
  res$stat_matched_med <- NA_real_
  res$stat_matched_mad <- NA_real_
  res$z_matched <- NA_real_
  res$n_matched <- 0L
  res$p_matched <- NA_real_
  for (s in have) {
    ch <- by_parent[[s]]
    if (is.null(ch)) next
    t_r <- stat_full[ch]
    t_r <- t_r[is.finite(t_r)]
    t_obs <- stat_full[[s]]
    res$n_matched[res$pathway == s] <- length(t_r)
    if (length(t_r) == 0L || !is.finite(t_obs)) next
    med <- stats::median(t_r)
    res$stat_matched_med[res$pathway == s] <- med
    mad <- stats::mad(t_r, center = med)
    res$stat_matched_mad[res$pathway == s] <- mad
    # robust standardized effect: the median can sit at or below zero for
    # centered statistics, which makes plain ratios sign-unstable; MAD is
    # non-negative and only degenerate when the reference class collapses
    res$z_matched[res$pathway == s] <- if (is.finite(mad) && mad > 0)
      (t_obs - med) / mad else NA_real_
    res$p_matched[res$pathway == s] <- (1 + sum(t_r >= t_obs)) /
      (length(t_r) + 1)
  }
  res$p_adj_BH <- stats::p.adjust(res$p_matched, method = "BH")
  res$n_spots <- nrow(A)
  class(res) <- c("svp_results", "data.frame")
  attr(res, "method") <- paste0(method, "_matched")
  attr(res, "block") <- !is.null(block)
  res
}
