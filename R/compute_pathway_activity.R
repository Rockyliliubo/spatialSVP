#' Compute pathway activity scores across spots
#'
#' Scores pathway activity per spot from a genes x spots expression matrix
#' using one of several back ends, and returns a \code{svp_activity} object
#' that \code{\link{test_spatial_variability}} consumes directly. The scoring
#' layer is deliberately thin: the package's contribution
#' is the spatial testing on top of it.
#'
#' @param x Genes x spots matrix, or a \code{SpatialExperiment},
#'   \code{SingleCellExperiment}, or \code{Seurat} object (see
#'   \code{\link{svp_input}}).
#' @param gene_sets Gene sets: a tidy data.frame with columns
#'   \code{source} (set id), \code{target} (gene id), and optional
#'   \code{mor} (mode of regulation / weight, default 1; negative weights
#'   indicate inhibitory members); or a named list of gene-id character
#'   vectors; or a genes x sets weight matrix.
#' @param coords Spot coordinates (see \code{\link{svp_input}}); stored in
#'   the returned object so later calls need not repeat it.
#' @param assay Assay to use for container classes.
#' @param method Scoring back end: \code{"mean_z"} (gene-wise z-scores
#'   weighted by \code{mor}, normalized by \code{sum(|mor|)} so weighted
#'   resources like PROGENy keep their scale), \code{"decoupler_mlm"} /
#'   \code{"decoupler_wmean"} (decoupleR multivariate linear model /
#'   weighted mean), \code{"fgsea"} (fgsea enrichment-score per spot on
#'   z-scored genes), or \code{"gsva"} (GSVA).
#' @param preprocess \code{"logcpm_z"} (log-CPM, then gene-wise z-scoring
#'   for the z-based back ends), \code{"logcpm"}, or \code{"none"} (input is
#'   already normalized; still z-scored for z-based back ends). Note: under
#'   library-size normalization, spatially over-expressed genes claim a
#'   larger share of affected spots' libraries, so unrelated gene sets can
#'   inherit coherent (negative) spatial structure -- a compositional
#'   coupling, not a testing artifact.
#' @param residualize Residualization of the log-CPM matrix before scoring:
#'   \code{"none"} (default), \code{"libsize"} (regress each gene's log-CPM
#'   profile on log10 library size; removes the shared library-size driver
#'   of the compositional coupling described above), or
#'   \code{"pcs"} (\code{"libsize"} plus regression on the top
#'   \code{n_pcs} expression principal components, i.e. dominant expression
#'   axes). \code{"pcs"} can reduce shared spatial background, but it can also
#'   subtract genuine biology when a pathway contributes to a leading
#'   component. Report sensitivity across residualization choices. Logical
#'   \code{TRUE}/\code{FALSE} is accepted for back-compatibility
#'   (\code{TRUE} -> \code{"libsize"}, deprecated). Requires count-scale
#'   input: an error with \code{preprocess = "none"}. Library size can
#'   carry real biology (region-level RNA content), so residualization is
#'   opt-in, not the default.
#' @param n_pcs Number of expression principal components regressed out when
#'   \code{residualize = "pcs"} (default 10).
#' @param min_size Minimum number of set genes present in the data for a set
#'   to be scored.
#' @param block Optional block (section/patient) factor, one value per spot;
#'   stored in the output and used by the testing functions for
#'   block-diagonal kernels and within-block permutations.
#' @param ... Passed to decoupleR or GSVA.
#'
#' @return A \code{svp_activity} object: a list with \code{activity}
#'   (spots x pathways matrix), \code{coords}, \code{block}, \code{method},
#'   \code{preprocess}, and the filtered \code{network}.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' dim(act$activity)
#'
#' @export
compute_pathway_activity <- function(x, gene_sets, coords = NULL, assay = NULL,
                                     method = c("mean_z", "decoupler_mlm",
                                                "decoupler_wmean", "fgsea",
                                                "gsva"),
                                     preprocess = c("logcpm_z", "logcpm",
                                                    "none"),
                                     residualize = c("none", "libsize",
                                                     "pcs"),
                                     n_pcs = 10L, min_size = 5L, block = NULL,
                                     ...) {
  method <- match.arg(method)
  preprocess <- match.arg(preprocess)
  if (is.logical(residualize)) {
    if (isTRUE(residualize)) {
      warning("residualize = TRUE is deprecated; use residualize = ",
              "'libsize'.", call. = FALSE)
    }
    residualize <- if (isTRUE(residualize)) "libsize" else "none"
  }
  residualize <- match.arg(residualize, c("none", "libsize", "pcs"))
  inp <- svp_input(x, coords = coords, assay = assay)
  n <- ncol(inp$counts)
  block <- .check_block(block, n)

  net <- .as_network(gene_sets)
  sets_order <- unique(as.character(net$source))
  net <- .filter_network(net, rownames(inp$counts), min_size = min_size)
  if (nrow(net) == 0L) {
    stop("No gene set retains at least ", min_size,
         " genes present in the data.", call. = FALSE)
  }

  mat <- if (preprocess %in% c("logcpm", "logcpm_z")) {
    m <- .logcpm(inp$counts)
    if (residualize != "none") m <- .resid_libsize(m, inp$counts)
    if (residualize == "pcs") m <- .resid_pcs(m, n_pcs)
    m
  } else {
    if (residualize != "none") {
      stop("'residualize' needs count-scale input: set preprocess ",
           "to 'logcpm_z' or 'logcpm'.", call. = FALSE)
    }
    inp$counts
  }

  A <- switch(method,
    mean_z = .score_mean_z(mat, net),
    decoupler_mlm = .score_decoupler(mat, net, "mlm", ...),
    decoupler_wmean = .score_decoupler(mat, net, "wmean", ...),
    fgsea = .score_fgsea(mat, net, ...),
    gsva = .score_gsva(mat, net, ...))

  A <- as.matrix(A)
  storage.mode(A) <- "double"
  A <- A[, intersect(sets_order, colnames(A)), drop = FALSE]
  if (any(!is.finite(A))) {
    bad <- colnames(A)[colSums(!is.finite(A)) > 0]
    warning("Non-finite activity scores in: ", paste(bad, collapse = ", "),
            "; consider 'none' preprocessing or another method.",
            call. = FALSE)
    A[!is.finite(A)] <- NA
  }

  out <- list(activity = A, coords = inp$coords, block = block,
              method = method, preprocess = preprocess,
              residualize = residualize, network = net,
              n_spots = n)
  class(out) <- "svp_activity"
  out
}

#' @export
print.svp_activity <- function(x, ...) {
  resid <- x$residualize
  resid <- if (isTRUE(resid)) "libsize" else as.character(resid)
  cat("svp_activity: ", nrow(x$activity), " spots x ", ncol(x$activity),
      " pathways (method = ", x$method, ", preprocess = ", x$preprocess,
      if (!identical(resid, "none")) paste0(", resid = ", resid) else "",
      ")\n", sep = "")
  invisible(x)
}

# ---- internal scoring helpers ---------------------------------------------

.as_network <- function(gene_sets) {
  if (is.data.frame(gene_sets)) {
    nms <- tolower(colnames(gene_sets))
    src <- intersect(c("source", "set", "pathway", "geneset"), nms)
    tgt <- intersect(c("target", "gene", "symbol"), nms)
    mor <- intersect(c("mor", "weight", "score"), nms)
    if (length(src) != 1L || length(tgt) != 1L) {
      stop("gene_sets data.frame needs columns like source/target(/mor).",
           call. = FALSE)
    }
    out <- data.frame(source = as.character(gene_sets[[src[1]]]),
                      target = as.character(gene_sets[[tgt[1]]]),
                      mor = if (length(mor) == 1L) as.numeric(gene_sets[[mor[1]]]) else 1)
  } else if (is.list(gene_sets)) {
    if (is.null(names(gene_sets))) {
      stop("A gene_sets list must be named.", call. = FALSE)
    }
    out <- do.call(rbind, lapply(names(gene_sets), function(s) {
      g <- gene_sets[[s]]
      data.frame(source = s, target = as.character(g), mor = 1)
    }))
  } else if (is.matrix(gene_sets)) {
    # genes x sets weight matrix
    lst <- lapply(colnames(gene_sets), function(s) {
      w <- gene_sets[, s]
      w <- w[w != 0]
      if (length(w) == 0L) return(NULL)
      data.frame(source = s, target = names(w), mor = unname(as.numeric(w)))
    })
    lst <- Filter(Negate(is.null), lst)
    if (length(lst) == 0L) {
      stop("gene_sets matrix has no non-zero weights.", call. = FALSE)
    }
    out <- do.call(rbind, lst)
  } else {
    stop("Unsupported gene_sets format.", call. = FALSE)
  }
  out$source <- factor(out$source)
  out <- out[is.finite(out$mor), , drop = FALSE]
  out
}

.filter_network <- function(net, genes, min_size = 5L) {
  net <- net[net$target %in% genes, , drop = FALSE]
  if (nrow(net) == 0L) return(net)
  # average duplicated (source, target) rows: run-length over the sorted
  # frame, linear-time (stats::aggregate is too slow at matched-null
  # scale, >= 1e7 edges)
  ord <- order(as.character(net$source), as.character(net$target))
  src <- as.character(net$source)[ord]
  tgt <- as.character(net$target)[ord]
  mor <- net$mor[ord]
  new_run <- c(TRUE, src[-1L] != src[-length(src)] |
                 tgt[-1L] != tgt[-length(tgt)])
  grp <- cumsum(new_run)
  ag_src <- src[new_run]
  ag_tgt <- tgt[new_run]
  ag_mor <- rowsum(mor, grp)[, 1L] / tabulate(grp)
  tab <- table(ag_src)
  keep <- names(tab)[tab >= min_size]
  sel <- ag_src %in% keep
  data.frame(source = factor(ag_src[sel]), target = ag_tgt[sel],
             mor = ag_mor[sel], stringsAsFactors = FALSE)
}

.score_mean_z <- function(mat, net) {
  z <- .zscore_rows(mat)
  # Vectorized scoring: with a sparse membership matrix M (genes x sets,
  # entries mor / sum|mor| per set) all set scores are ONE crossproduct.
  # A per-set loop costs ~1 ms per set x gene, i.e. an hour on
  # matched-null-scale libraries (1e5 sets); this costs under a minute.
  sets <- unique(as.character(net$source))
  src <- as.character(net$source)
  gi <- match(as.character(net$target), rownames(z))
  ok <- !is.na(gi)
  gi <- gi[ok]
  src <- src[ok]
  mor <- net$mor[ok]
  denom <- rowsum(abs(mor), group = src)[, 1]  # rownames sorted by name
  M <- Matrix::sparseMatrix(i = gi, j = match(src, sets),
                            x = mor / denom[src],
                            dims = c(nrow(z), length(sets)),
                            dimnames = list(rownames(z), sets))
  A <- as.matrix(Matrix::crossprod(z, M))
  # Return NaN for sets with no present genes or all-zero weights.
  bad <- union(names(denom)[denom == 0], setdiff(sets, names(denom)))
  if (length(bad) > 0L) A[, bad] <- NaN
  A
}

.score_decoupler <- function(mat, net, statistic, ...) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required for method 'decoupler_*'.",
         call. = FALSE)
  }
  net2 <- data.frame(source = as.character(net$source),
                     target = net$target, mor = net$mor)
  z <- .zscore_rows(mat)
  res <- decoupleR::decouple(mat = z, network = net2, statistics = statistic,
                             minsize = 0, ...)
  res <- as.data.frame(res)
  res <- res[res$statistic == statistic, , drop = FALSE]
  if (nrow(res) == 0L) stop("decoupleR returned no rows.", call. = FALSE)
  conds <- colnames(mat)
  srcs <- unique(res$source)
  A <- matrix(NA_real_, length(conds), length(srcs),
              dimnames = list(conds, srcs))
  A[cbind(match(res$condition, conds), match(res$source, srcs))] <- res$score
  drop_cols <- colSums(is.na(A)) == nrow(A)
  if (any(drop_cols)) {
    warning("decoupleR dropped sets (its own size filtering): ",
            paste(colnames(A)[drop_cols], collapse = ", "), call. = FALSE)
    A <- A[, !drop_cols, drop = FALSE]
  }
  A
}

.score_gsva <- function(mat, net, ...) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for method 'gsva'.", call. = FALSE)
  }
  agg <- stats::aggregate(mor ~ source + target, data = net, FUN = mean)
  gsets <- split(agg$target, as.character(agg$source))
  has_param <- "gsvaParam" %in% getNamespaceExports("GSVA")
  es <- if (has_param) {
    gp <- GSVA::gsvaParam(exprData = mat, geneSets = gsets, ...)
    GSVA::gsva(gp, verbose = FALSE)
  } else {
    GSVA::gsva(mat, gsets, verbose = FALSE, ...)
  }
  A <- t(as.matrix(es))
  A[rownames(A) %in% colnames(mat), , drop = FALSE]
}

# Per-spot fgsea enrichment scores: each spot's genes are ranked by their
# (z-scored) expression and the standard running-sum ES is computed per set
# via fgsea's permutation-free calcGseaStat(). Verified identical to
# fgseaSimple() output. Rank-based: the mor weights are not used (set
# membership only). Positive ES = set members concentrated at the top of
# the spot's expression ranking.
.score_fgsea <- function(mat, net, gseaParam = 1, ...) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for method 'fgsea'.", call. = FALSE)
  }
  z <- .zscore_rows(mat)
  rn <- rownames(z)
  agg <- stats::aggregate(mor ~ source + target, data = net, FUN = mean)
  members <- split(as.character(agg$target), as.character(agg$source))
  idx <- lapply(members, function(g) {
    i <- match(g, rn)
    i[!is.na(i)]
  })
  sets <- names(idx)
  A <- matrix(NA_real_, nrow = ncol(z), ncol = length(sets),
              dimnames = list(colnames(z), sets))
  for (j in seq_len(ncol(z))) {
    zj <- z[, j]
    o <- order(zj, decreasing = TRUE)
    posmap <- integer(length(zj))
    posmap[o] <- seq_along(o)
    rz <- zj[o]
    A[j, ] <- vapply(idx, function(ii) {
      fgsea::calcGseaStat(rz, sort(posmap[ii]), gseaParam = gseaParam)
    }, numeric(1))
  }
  A
}
