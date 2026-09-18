# Internal helpers ----------------------------------------------------------

# ggplot2 aesthetic columns are NSE; silence R CMD check's global-variable NOTE
utils::globalVariables(c("x", "y", "value", "stat", "neglog10p",
                         "significant", "rank", "scale", "lab"))

`%||%` <- function(a, b) if (is.null(a)) b else a

.zscore_rows <- function(m) {
  mu <- rowMeans(m)
  sdev <- apply(m, 1L, stats::sd)
  bad <- !is.finite(sdev) | sdev == 0
  sdev[bad] <- 1
  out <- (m - mu) / sdev
  out[!is.finite(out)] <- 0
  out
}

.logcpm <- function(counts) {
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  log1p(t(t(counts) / lib) * 1e4)
}

# Regress each gene's (log-CPM, already library-residualized) profile on the
# top expression principal components -- the tissue-composition axes -- and
# return the residuals. PCs come from the eigendecomposition of the
# spots x spots Gram matrix (exactly the SVD principal scores up to scale,
# which leaves the projection unchanged), so no genes x spots transpose is
# ever materialized. Gene means are kept intact: PC scores are orthogonal
# to the constant direction because the Gram is built from gene-centered
# profiles.
.resid_pcs <- function(logcpm, n_pcs = 10L) {
  n <- ncol(logcpm)
  k <- min(as.integer(n_pcs), n - 1L)
  if (!is.finite(k) || k < 1L) return(logcpm)
  mc <- logcpm - rowMeans(logcpm)
  scores <- .pc_scores(mc, k)
  coef <- mc %*% scores                    # genes x k (scores orthonormal)
  logcpm - tcrossprod(coef, scores)
}

# Top-k principal-component spot scores of a row-centred genes x spots
# matrix. The full eigendecomposition of the spots x spots Gram matrix is
# O(n_spots^3) (~40 min at 14k spots); beyond 8000 spots switch to
# RSpectra's Lanczos iteration on the Gram OPERATOR (x -> mc' (mc x)),
# which yields the exact top-k subspace in seconds without materializing
# the Gram. Falls back to the full eigen when RSpectra is unavailable.
.pc_scores <- function(mc, k) {
  n <- ncol(mc)
  if (n > 8000L && requireNamespace("RSpectra", quietly = TRUE)) {
    matop <- function(x, args) as.vector(crossprod(mc, mc %*% x))
    RSpectra::eigs_sym(matop, k = k, n = n)$vectors
  } else {
    eigen(crossprod(mc), symmetric = TRUE)$vectors[, seq_len(k),
                                                  drop = FALSE]
  }
}

# Regress each gene's (log-CPM) profile on log10 library size and return the
# residuals: removes the shared library-size driver that makes unrelated gene
# sets inherit spatial structure under compositional normalization.
# With zero library-size variance (degenerate input) the centered matrix is
# returned; downstream z-scoring then makes residualization a no-op.
.resid_libsize <- function(logcpm, counts) {
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  xc <- log10(lib)
  xc <- xc - mean(xc)
  sxx <- sum(xc * xc)
  M <- sweep(logcpm, 1L, rowMeans(logcpm), "-")
  if (!is.finite(sxx) || sxx == 0) return(M)
  beta <- as.vector(M %*% xc) / sxx
  M - outer(beta, xc)
}

# Validate/extract two numeric coordinate columns. Columns named 'x' and 'y'
# (case-insensitive) are preferred; otherwise the first two columns are used.
.coords_xy <- function(coords) {
  if (is.null(coords)) {
    stop("'coords' is required: a data.frame or matrix with two numeric ",
         "columns (columns named 'x'/'y' are used if present, otherwise ",
         "the first two columns).", call. = FALSE)
  }
  if (is.matrix(coords)) coords <- as.data.frame(coords)
  if (!is.data.frame(coords)) {
    stop("'coords' must be a data.frame or matrix.", call. = FALSE)
  }
  if (ncol(coords) < 2L) {
    stop("'coords' needs at least two columns.", call. = FALSE)
  }
  nms <- tolower(colnames(coords))
  ix <- which(nms == "x")
  iy <- which(nms == "y")
  if (length(ix) == 1L && length(iy) == 1L) {
    out <- data.frame(x = as.numeric(coords[[ix]]),
                      y = as.numeric(coords[[iy]]),
                      row.names = rownames(coords))
  } else {
    out <- data.frame(x = as.numeric(coords[[1L]]),
                      y = as.numeric(coords[[2L]]),
                      row.names = rownames(coords))
  }
  if (anyNA(out)) stop("'coords' contains missing values.", call. = FALSE)
  out
}

.pair_dist2 <- function(coords) {
  dx <- outer(coords$x, coords$x, "-")
  dy <- outer(coords$y, coords$y, "-")
  dx * dx + dy * dy
}

.same_block_mask <- function(block) {
  if (is.null(block)) NULL else outer(block, block, "==")
}

.sqexp_kernel <- function(coords, ls, mask = NULL) {
  d2 <- .pair_dist2(coords)
  if (!is.null(mask)) d2[!mask] <- Inf
  exp(-d2 / (2 * ls * ls))
}

.nndist_median <- function(coords, mask = NULL) {
  d2 <- .pair_dist2(coords)
  diag(d2) <- Inf
  if (!is.null(mask)) d2[!mask] <- Inf
  mins <- apply(d2, 1L, min)
  vals <- mins[is.finite(mins)]
  if (length(vals) == 0L) {
    stop("Cannot compute within-block distances: every block has fewer than ",
         "two spots.", call. = FALSE)
  }
  # sqrt: .pair_dist2 stores squared distances; the lengthscale grids
  # consume distance units, so return the distance itself.
  stats::median(sqrt(vals))
}

.coord_range <- function(coords, block = NULL) {
  span <- function(cx, cy) max(diff(range(cx)), diff(range(cy)))
  if (is.null(block)) return(span(coords$x, coords$y))
  max(vapply(split(seq_len(nrow(coords)), block),
             function(i) span(coords$x[i], coords$y[i]), numeric(1)))
}

# Lengthscale candidate grid in the units of the coordinates: from the median
# nearest-neighbour distance up to half the within-block coordinate span.
.ls_grid <- function(coords, block = NULL, n = 6L) {
  mask <- .same_block_mask(block)
  l0 <- .nndist_median(coords, mask)
  l1 <- max(.coord_range(coords, block) / 2, l0 * 4)
  unique(pmax(l0, exp(seq(log(l0), log(l1), length.out = n))))
}

# Symmetrized binary k-nearest-neighbour adjacency.
.adjacency_w <- function(coords, k = 6L, mask = NULL) {
  n <- nrow(coords)
  d2 <- .pair_dist2(coords)
  diag(d2) <- Inf
  if (!is.null(mask)) d2[!mask] <- Inf
  kk <- min(k, n - 1L)
  ii <- integer(0)
  jj <- integer(0)
  for (i in seq_len(n)) {
    nn <- order(d2[i, ])[seq_len(kk)]
    nn <- nn[is.finite(d2[i, nn])]
    if (length(nn) > 0L) {
      ii <- c(ii, rep.int(i, length(nn)))
      jj <- c(jj, nn)
    }
  }
  # sparse: a dense n x n adjacency does not scale to multi-section
  # pools (1e4+ spots); the only downstream ops are sum(W) and W %*% Z
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = 1.0, dims = c(n, n))
  W <- W + Matrix::t(W)
  Matrix::drop0(1 * (W != 0))
}

# Block-aware index permutations: spots are shuffled only within blocks.
.make_perms <- function(n, block = NULL, n_perm) {
  lapply(seq_len(n_perm), function(...) .permute_index(n, block))
}

.permute_index <- function(n, block = NULL) {
  if (is.null(block)) return(sample.int(n))
  idx <- seq_len(n)
  out <- integer(n)
  for (b in unique(block)) {
    wi <- which(block == b)
    out[wi] <- idx[sample(wi)]
  }
  out
}

.check_block <- function(block, n) {
  if (is.null(block)) return(NULL)
  block <- as.factor(block)
  if (length(block) != n) {
    stop("length(block) (", length(block), ") must match the number of spots ",
         "(", n, ").", call. = FALSE)
  }
  if (anyNA(block)) stop("'block' contains NA.", call. = FALSE)
  block
}
