#' Test pathway activity for spatial variability
#'
#' Tests, for every pathway, whether its per-spot activity vector shows
#' statistically significant spatial structure. Three back ends are provided:
#' \itemize{
#' \item \code{"morans_i"}: Moran's I on a symmetrized k-nearest-neighbour
#'   graph; p-values from block-aware permutations (default \code{n_perm}).
#' \item \code{"covariance"}: SPARK-X-inspired multi-scale squared-exponential
#'   covariance statistics; the max is taken over kernel scales AFTER
#'   per-scale studentization by the pooled permutation null moments,
#'   because the raw statistic's variance differs by orders of magnitude
#'   across scales and the unstandardized argmax tracks kernel geometry, not
#'   signal. Reports the dominant scale -- a coarse
#'   descriptor (about 2 grid points of resolution); the sharp per-pathway
#'   length scale estimator is the \code{"gp"} back end.
#' \item \code{"gp"}: SpatialDE-inspired Gaussian-process likelihood-ratio
#'   test with a lengthscale x mixing-proportion grid (the grid includes
#'   rho = 0 so the alternative nests the null); p-values from block-aware
#'   permutations like the other back ends; reports per-pathway lengthscale
#'   and spatial variance fraction.
#' }
#'
#' When \code{block} is set (sections/patients), spatial kernels are
#' block-diagonal and permutations are restricted within blocks.
#'
#' @details
#' All three back ends test ONE-SIDED alternatives: positive spatial
#' autocorrelation at scales at or above the spot spacing, which is the
#' relevant signature for tissue-region biology. High-frequency or
#' checkerboard-like patterns (negative autocorrelation) are deliberately
#' out of scope; a perfectly anticorrelated pattern scores low, not high.
#'
#' @param activity A \code{svp_activity} object from
#'   \code{\link{compute_pathway_activity}}, or a spots x pathways matrix
#'   (then \code{coords} is required).
#' @param coords Spot coordinates; defaults to the object's stored coords.
#' @param method Testing back end (see above).
#' @param n_perm Number of shared permutations. Default 199 for all back
#'   ends (p-values are permutation-calibrated). For
#'   \code{"gp"} only, \code{n_perm = 0} explicitly selects the analytic
#'   50:50 chi-square mixture p-value, which is uncalibrated for the max
#'   over the lengthscale x rho grid and is retained for diagnostics only
#'   (a warning is emitted).
#' @param block Optional block factor (section/patient) per spot; defaults to
#'   the object's stored block.
#' @param seed Optional random seed for reproducible permutations.
#'
#' @return A data.frame of class \code{svp_results} with columns
#'   \code{pathway}, \code{stat}, \code{p_value}, \code{p_adj_BH}, plus
#'   back-end-specific columns (\code{lengthscale},
#'   \code{spatial_var_fraction} for \code{"gp"}; \code{dominant_scale} for
#'   \code{"covariance"}). For \code{"gp"}, when the intercept-only null
#'   model wins (no grid point beats it), \code{spatial_var_fraction} is 0
#'   and \code{lengthscale} is \code{NA} -- there is no spatial scale to
#'   report.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' res <- test_spatial_variability(act, method = "morans_i",
#'                                 n_perm = 39, seed = 1)
#' res[order(res$p_value), ]
#'
#' @export
test_spatial_variability <- function(activity, coords = NULL,
                                     method = c("morans_i", "covariance",
                                                "gp"),
                                     n_perm = NULL, block = NULL,
                                     seed = NULL) {
  method <- match.arg(method)
  A <- .extract_activity(activity, coords, block)
  coords <- A$coords
  block <- A$block
  A <- A$A
  if (!is.null(seed)) set.seed(seed)
  if (is.null(n_perm)) n_perm <- 199L

  pathways <- colnames(A) %||% paste0("pathway", seq_len(ncol(A)))
  colnames(A) <- pathways
  sd_col <- apply(A, 2L, stats::sd)
  keep <- is.finite(sd_col) & sd_col > 0
  if (!any(keep)) stop("No pathway has non-zero variance.", call. = FALSE)
  Z <- scale(A[, keep, drop = FALSE])

  mask <- .same_block_mask(block)
  out <- switch(
    method,
    morans_i = {
      if (n_perm < 19L) {
        stop("method 'morans_i' needs n_perm >= 19 (permutation-based ",
             "p-values only).", call. = FALSE)
      }
      W <- .adjacency_w(coords, k = 6L, mask = mask)
      S0 <- sum(W)
      if (S0 == 0) stop("No spot pairs within blocks to build adjacency.",
                        call. = FALSE)
      stat_obs <- .stat_morans(Z, W, S0)
      null <- .perm_shared_stats(Z, coords, block, mask, n_perm,
                                 method = "morans_i", setup = list(W = W, S0 = S0))
      p_raw <- .perm_pvals(stat_obs, null)
      data.frame(stat = stat_obs, p_value = p_raw,
                 lengthscale = NA_real_, spatial_var_fraction = NA_real_,
                 dominant_scale = NA_real_)
    },
    covariance = {
      if (n_perm < 19L) {
        stop("method 'covariance' needs n_perm >= 19 (permutation-based ",
             "p-values only).", call. = FALSE)
      }
      setup <- .cov_setup(coords, block)
      oc <- .cov_studentized(Z, block, n_perm, setup)
      p_raw <- .perm_pvals(oc$stat, oc$null)
      data.frame(stat = oc$stat, p_value = p_raw,
                 lengthscale = NA_real_, spatial_var_fraction = NA_real_,
                 dominant_scale = setup$scales[oc$scale_idx])
    },
    gp = {
      if (nrow(Z) > 5000L) {
        warning("gp back end with > 5000 spots: eigendecompositions are ",
                "memory/time heavy.", call. = FALSE)
      }
      setup <- .gp_setup(coords, block)
      og <- .gp_stats(Z, setup)
      if (n_perm > 0L) {
        if (n_perm < 19L) {
          stop("method 'gp' needs n_perm >= 19 (permutation-based ",
               "p-values only).", call. = FALSE)
        }
        null <- .perm_shared_stats(Z, coords, block, mask, n_perm,
                                   method = "gp", setup = setup)
        p_raw <- .perm_pvals(og$stat, null)
      } else {
        warning("gp analytic p-values are uncalibrated for the max over ",
                "the lengthscale x rho grid; retained for diagnostics ",
                "only. Use n_perm >= 19 for inference.",
                call. = FALSE)
        p_raw <- og$p_value
      }
      data.frame(stat = og$stat, p_value = p_raw,
                 lengthscale = og$lengthscale,
                 spatial_var_fraction = og$rho,
                 dominant_scale = NA_real_)
    })

  res <- data.frame(pathway = pathways, stat = NA_real_,
                    p_value = NA_real_, stringsAsFactors = FALSE)
  res$stat[keep] <- out$stat
  res$p_value[keep] <- out$p_value
  res$p_adj_BH <- stats::p.adjust(res$p_value, method = "BH")
  res$lengthscale <- NA_real_
  res$lengthscale[keep] <- out$lengthscale
  res$spatial_var_fraction <- NA_real_
  res$spatial_var_fraction[keep] <- out$spatial_var_fraction
  res$dominant_scale <- NA_real_
  res$dominant_scale[keep] <- out$dominant_scale
  res$n_spots <- nrow(A)
  class(res) <- c("svp_results", "data.frame")
  attr(res, "method") <- method
  attr(res, "block") <- !is.null(block)
  res
}

#' @export
print.svp_results <- function(x, ...) {
  # robust to arbitrary row/column subsets (class travels with `[`)
  fdr_col <- intersect(c("p_adj_BH", "p_adj_wy"), names(x))
  sig <- if (length(fdr_col)) {
    sum(x[[fdr_col[1]]] < 0.05, na.rm = TRUE)
  } else NA_integer_
  cat("svp_results: '", attr(x, "method"), "' test on ", nrow(x),
      " pathways; ", sig, " significant at FDR 0.05\n", sep = "")
  p_col <- intersect(c("p_value", "p_perm", "p_matched"), names(x))
  ord <- if (length(p_col)) order(x[[p_col[1]]]) else seq_len(nrow(x))
  head_cols <- intersect(c("pathway", "stat", "p_value", "p_perm",
                           "p_matched", "p_adj_BH", "p_adj_wy", "lengthscale",
                           "dominant_scale"), names(x))
  if (length(head_cols) > 0L) {
    top <- x[ord, head_cols, drop = FALSE]
    class(top) <- "data.frame"  # avoid re-dispatching this print method
    print(utils::head(top, 10L))
  }
  invisible(x)
}

# ---- shared internals ------------------------------------------------------

.extract_activity <- function(activity, coords = NULL, block = NULL) {
  if (inherits(activity, "svp_activity")) {
    if (is.null(coords)) coords <- activity$coords
    if (is.null(block)) block <- activity$block
    A <- activity$activity
  } else {
    A <- tryCatch(as.matrix(activity), error = function(e) NULL)
    if (is.null(A)) {
      stop("'activity' must be a svp_activity object or a spots x pathways ",
           "matrix.", call. = FALSE)
    }
  }
  coords <- .coords_xy(coords)
  if (nrow(A) != nrow(coords)) {
    stop("activity rows (", nrow(A), ") must equal coordinate rows (",
         nrow(coords), "): activity must be spots x pathways.", call. = FALSE)
  }
  if (anyNA(A)) {
    stop("'activity' contains NA; remove the affected spots or pathways ",
         "first.", call. = FALSE)
  }
  block <- .check_block(block, nrow(A))
  list(A = A, coords = coords, block = block)
}

.stat_morans <- function(Z, W, S0) {
  # Conventional Moran's I: Z comes from scale(), so sum(z^2) = n - 1 and
  # the raw ratio is the textbook I times (n - 1)/n; rescale to the
  # textbook definition so 'stat' compares directly with ape::Moran.I,
  # Seurat, etc. (p-values are unaffected: strictly monotone transform).
  n <- nrow(Z)
  colSums(Z * as.matrix(W %*% Z)) / S0 * n / (n - 1)
}

.cov_setup <- function(coords, block = NULL) {
  mask <- .same_block_mask(block)
  l0 <- .nndist_median(coords, mask)
  l1 <- max(.coord_range(coords, block) / 2, l0 * 8)
  scales <- unique(exp(seq(log(l0), log(l1), length.out = 5L)))
  kerns <- lapply(scales, function(ls) .sqexp_kernel(coords, ls, mask))
  list(kernels = kerns, scales = scales)
}

# Raw per-scale covariance statistics T = z'Kz (pathways x scales). vapply
# drops the pathway dimension when ncol(Z) == 1; rebuild explicitly
# with nrow = 1.
.cov_stats_scales <- function(Z, setup) {
  out <- vapply(setup$kernels, function(K) colSums(Z * (K %*% Z)),
                numeric(ncol(Z)))
  if (!is.matrix(out)) out <- matrix(out, nrow = 1L)
  out
}

# Max over scales of the per-scale STANDARDIZED statistic.
# E[T] = trace(K) is scale-independent, but Var(T) = 2 trace(K^2) differs by
# orders of magnitude across scales, so the raw argmax is decided by kernel
# geometry, not by signal. This stats-only
# variant standardizes per scale by the across-pathway median/MAD (used by
# matched_enrichment(), where most columns are null-like matched
# replicates); the permutation-based entry points use .cov_studentized()
# with moments pooled over the shared permutation draws instead.
.std_covmax <- function(raw) {
  ctr <- apply(raw, 2L, stats::median)
  scl <- apply(raw, 2L, stats::mad)
  scl[!is.finite(scl) | scl <= 0] <- 1
  std <- sweep(sweep(raw, 2L, ctr, "-"), 2L, scl, "/")
  kmax <- max.col(std, ties.method = "first")
  list(stat = std[cbind(seq_len(nrow(std)), kmax)], scale_idx = kmax)
}

.stat_covmax <- function(Z, setup) {
  oc <- .std_covmax(.cov_stats_scales(Z, setup))
  list(stat = oc$stat, scale = setup$scales[oc$scale_idx])
}

# Permutation-based covariance core, shared by test_spatial_variability()
# and permutation_fdr(): per-scale null moments are pooled over the
# shared permutation draws and all pathways (null exchangeable across
# pathways for unit-variance columns); the max over scales is taken on
# studentized statistics, and every permutation draw is transformed by the
# same pooled moments, keeping the permutation comparison symmetric.
.cov_studentized <- function(Z, block, n_perm, setup) {
  n <- nrow(Z)
  p <- ncol(Z)
  n_scales <- length(setup$scales)
  null_raw <- array(NA_real_, dim = c(n_perm, p, n_scales))
  for (b in seq_len(n_perm)) {
    Zp <- Z[.permute_index(n, block), , drop = FALSE]
    null_raw[b, , ] <- .cov_stats_scales(Zp, setup)
  }
  mu <- apply(null_raw, 3L, mean)
  sdv <- apply(null_raw, 3L, stats::sd)
  sdv[!is.finite(sdv) | sdv <= 0] <- 1
  studentize <- function(M) sweep(sweep(M, 2L, mu, "-"), 2L, sdv, "/")
  so <- studentize(.cov_stats_scales(Z, setup))
  kmax <- max.col(so, ties.method = "first")
  null_std <- matrix(NA_real_, n_perm, p)
  for (b in seq_len(n_perm)) {
    sn <- studentize(null_raw[b, , ])
    null_std[b, ] <- sn[cbind(seq_len(p), max.col(sn, ties.method = "first"))]
  }
  list(stat = so[cbind(seq_len(p), kmax)], scale_idx = kmax,
       null = null_std)
}

.gp_setup <- function(coords, block = NULL) {
  mask <- .same_block_mask(block)
  ls_grid <- .ls_grid(coords, block)
  eig <- lapply(ls_grid, function(ls) {
    K <- .sqexp_kernel(coords, ls, mask)
    e <- eigen(K, symmetric = TRUE)
    list(values = e$values, vectors = e$vectors)
  })
  list(ls_grid = ls_grid, eig = eig,
       # rho = 0 makes the alternative nest the null (no LRT truncation at
       # p = 0.5); near-0 values refine the boundary region.
       rho_grid = c(0, 0.001, 0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.9))
}

.gp_stats <- function(Y, setup) {
  n <- nrow(Y)
  p <- ncol(Y)
  # Reference null: intercept-only model (rho = 0). It is identical for
  # every lengthscale, so it is fitted ONCE here instead of occupying the
  # grid (a ridge of tied rho = 0 fits made the argmax numerically
  # unstable). The grid alternative is a superset of this null, so the LRT
  # is >= 0 by construction and never truncates into a p = 0.5 mass point.
  mu0 <- colMeans(Y)
  v0 <- pmax(colMeans(Y * Y) - mu0 * mu0, 1e-12)
  ll0 <- -n / 2 * (log(2 * pi * v0) + 1)
  llbest <- ll0
  rbest <- numeric(p)         # rho = 0 when no grid point beats the null
  lbest <- rep(NA_real_, p)   # no spatial scale when the null wins
  rhos <- setup$rho_grid[setup$rho_grid > 0]
  for (li in seq_along(setup$eig)) {
    eg <- setup$eig[[li]]
    Yrot <- crossprod(eg$vectors, Y)
    oner <- drop(crossprod(eg$vectors, rep(1, n)))
    for (rho in rhos) {
      w <- 1 / (rho * eg$values + (1 - rho))
      oneV1 <- sum(w * oner * oner)
      yVone <- colSums(w * oner * Yrot)
      yVy <- colSums(w * Yrot * Yrot)
      mu <- yVone / oneV1
      v <- (yVy - 2 * mu * yVone + mu * mu * oneV1) / n
      v <- pmax(v, 1e-12)
      ll <- -n / 2 * (log(2 * pi * v) + 1)
      better <- ll > llbest
      llbest[better] <- ll[better]
      rbest[better] <- rho
      lbest[better] <- setup$ls_grid[li]
    }
  }
  lrt <- 2 * pmax(llbest - ll0, 0)
  # NOTE: the 50:50 chi-square mixture is the null of a SINGLE boundary
  # variance component, not of the max over this grid; analytic p-values
  # are diagnostic only. Inference uses permutations.
  pval <- 0.5 * stats::pchisq(lrt, 1, lower.tail = FALSE)
  list(stat = lrt, p_value = pval, rho = rbest, lengthscale = lbest,
       ll = llbest)
}

# Apply ONE shared block-aware permutation to all pathways per replicate and
# return the null statistics (B x pathways). Shared permutations are what
# preserves cross-pathway correlation for the Westfall-Young machinery.
.perm_shared_stats <- function(Z, coords, block, mask, n_perm, method,
                               setup) {
  n <- nrow(Z)
  p <- ncol(Z)
  null <- matrix(NA_real_, n_perm, p)
  for (b in seq_len(n_perm)) {
    pi_b <- .permute_index(n, block)
    Zp <- Z[pi_b, , drop = FALSE]
    null[b, ] <- switch(
      method,
      morans_i = .stat_morans(Zp, setup$W, setup$S0),
      covariance = .stat_covmax(Zp, setup)$stat,
      gp = .gp_stats(Zp, setup)$stat)
  }
  null
}

.perm_pvals <- function(stat_obs, stat_null) {
  B <- nrow(stat_null)
  ge <- sweep(stat_null, 2L, stat_obs, `>=`)
  (1 + colSums(ge)) / (B + 1)
}
