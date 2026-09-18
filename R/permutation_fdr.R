#' Gene-set-level permutation testing (Westfall-Young maxT)
#'
#' Pathway activities derived from one shared expression matrix are highly
#' correlated (shared genes), and hundreds of pathways are tested at once.
#' Benjamini-Hochberg control depends on conditions on the true-null p-values
#' and their dependence. This function applies ONE shared, block-aware
#' coordinate permutation to all
#' pathways per replicate and computes single-step Westfall-Young maxT
#' adjusted p-values, whose joint null automatically reproduces the empirical
#' cross-pathway correlation structure.
#'
#' @param activity A \code{svp_activity} object or a spots x pathways matrix
#'   (with \code{coords}).
#' @param coords Spot coordinates; defaults to the object's stored coords.
#' @param method Statistic to permute: \code{"morans_i"} (default, fastest),
#'   \code{"covariance"}, or \code{"gp"} (slow: re-fits the GP grid per
#'   permutation).
#' @param n_perm Number of shared permutations (default 499).
#' @param block Optional block factor (section/patient) per spot; permutations
#'   shuffle spots only within blocks.
#' @param seed Optional random seed.
#'
#' @return A data.frame of class \code{svp_results} with columns
#'   \code{pathway}, \code{stat}, \code{p_perm} (per-pathway permutation
#'   p-value), \code{p_adj_BH} (BH on \code{p_perm}), and
#'   \code{p_adj_wy} (Westfall-Young single-step adjusted p-value).
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' fdr <- permutation_fdr(act, method = "morans_i", n_perm = 49, seed = 1)
#' fdr[order(fdr$p_perm), ]
#'
#' @export
permutation_fdr <- function(activity, coords = NULL,
                            method = c("morans_i", "covariance", "gp"),
                            n_perm = 499L, block = NULL, seed = NULL) {
  method <- match.arg(method)
  ex <- .extract_activity(activity, coords, block)
  A <- ex$A
  coords <- ex$coords
  block <- ex$block
  if (!is.null(seed)) set.seed(seed)
  if (n_perm < 19L) stop("Use n_perm >= 19.", call. = FALSE)

  pathways <- colnames(A) %||% paste0("pathway", seq_len(ncol(A)))
  colnames(A) <- pathways
  sd_col <- apply(A, 2L, stats::sd)
  keep <- is.finite(sd_col) & sd_col > 0
  Z <- scale(A[, keep, drop = FALSE])

  mask <- .same_block_mask(block)
  dom_scale <- NULL
  if (identical(method, "covariance")) {
    # studentized multi-scale core: per-scale null moments pooled
    # over the shared permutation draws; max over scales on studentized
    # statistics; symmetric transform for observed and null draws.
    cov_setup <- .cov_setup(coords, block)
    cs <- .cov_studentized(Z, block, n_perm, cov_setup)
    stat_obs <- cs$stat
    null <- cs$null
    dom_scale <- cov_setup$scales[cs$scale_idx]
  } else {
    setup <- switch(method,
      morans_i = {
        W <- .adjacency_w(coords, k = 6L, mask = mask)
        list(W = W, S0 = sum(W))
      },
      gp = {
        if (nrow(Z) > 3000L) {
          warning("gp with ", nrow(Z), " spots x ", n_perm,
                  " permutations is slow; consider method='morans_i'.",
                  call. = FALSE)
        }
        .gp_setup(coords, block)
      })

    stat_one <- function(Zm) {
      switch(method,
        morans_i = .stat_morans(Zm, setup$W, setup$S0),
        gp = .gp_stats(Zm, setup)$stat)
    }
    stat_obs <- stat_one(Z)
    null <- .perm_shared_stats(Z, coords, block, mask, n_perm,
                               method = method, setup = setup)
  }

  p_perm <- .perm_pvals(stat_obs, null)
  maxB <- apply(null, 1L, max, na.rm = TRUE)
  p_wy <- (1 + colSums(maxB >= matrix(stat_obs, n_perm, length(stat_obs),
                                      byrow = TRUE))) / (n_perm + 1)

  res <- data.frame(pathway = pathways, stat = NA_real_,
                    p_perm = NA_real_, stringsAsFactors = FALSE)
  res$stat[keep] <- stat_obs
  res$p_perm[keep] <- p_perm
  res$p_adj_BH <- stats::p.adjust(res$p_perm, method = "BH")
  res$p_adj_wy <- NA_real_
  res$p_adj_wy[keep] <- p_wy
  if (!is.null(dom_scale)) {
    res$dominant_scale <- NA_real_
    res$dominant_scale[keep] <- dom_scale
  }
  res$n_spots <- nrow(A)
  class(res) <- c("svp_results", "data.frame")
  attr(res, "method") <- paste0(method, "_wy")
  attr(res, "block") <- !is.null(block)
  res
}
