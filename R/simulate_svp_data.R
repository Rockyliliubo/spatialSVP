#' Simulate spatial transcriptomics data with spatially variable gene sets
#'
#' Generates a toy spatial transcriptomics experiment on a jittered grid:
#' a fraction of gene sets carries a shared spatial signal (Gaussian bump or
#' linear gradient), the rest are spatially null. Used throughout the package
#' tests and vignette for power / type-I-error calibration.
#'
#' @param n_spots Number of spots (placed on the smallest square grid that
#'   fits, then jittered).
#' @param n_genes Total number of genes; must be at least
#'   \code{n_sets * set_size}.
#' @param n_sets Total number of gene sets.
#' @param n_spatial_sets Number of sets carrying a spatial signal.
#' @param set_size Number of genes per set.
#' @param n_blocks Number of coordinate blocks (vertical strips); spots are
#'   shuffled/analyzed within blocks when the returned \code{block} is passed
#'   on. 1 = single section.
#' @param pattern Signal shape: \code{"bump"}, \code{"gradient"}, or
#'   \code{"mix"} (half bump, half gradient across the spatial sets).
#' @param effect Signal amplitude multiplier (log space).
#' @param noise Gene-level log-space noise sd.
#' @param coord_scale Multiplier applied to the grid coordinates (default 1,
#'   i.e. unit spacing). Use e.g. \code{276} to mimic Visium pixel
#'   coordinates; signal shapes and block assignments are defined on the
#'   unit grid and scale along with it, so results are directly comparable
#'   across coordinate magnitudes.
#' @param seed Optional random seed.
#'
#' @return A list with \code{counts} (genes x spots matrix),
#'   \code{coords} (data.frame), \code{gene_sets} (tidy data.frame:
#'   source/target/mor), \code{block} (factor or NULL), and \code{truth}
#'   (data.frame: set, spatial, pattern, center_x, center_y).
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' table(sim$truth$spatial)
#'
#' @export
simulate_svp_data <- function(n_spots = 900L, n_genes = 120L, n_sets = 12L,
                              n_spatial_sets = 4L, set_size = 10L,
                              n_blocks = 1L, pattern = c("bump", "gradient",
                                                         "mix"),
                              effect = 1.5, noise = 0.6, coord_scale = 1,
                              seed = NULL) {
  pattern <- match.arg(pattern)
  n_spatial_sets <- as.integer(n_spatial_sets)
  n_sets <- as.integer(n_sets)
  set_size <- as.integer(set_size)
  n_genes <- as.integer(n_genes)
  n_spots <- as.integer(n_spots)
  if (n_genes < n_sets * set_size) {
    stop("n_genes (", n_genes, ") must be >= n_sets * set_size (",
         n_sets * set_size, ").", call. = FALSE)
  }
  if (n_spatial_sets > n_sets) {
    stop("n_spatial_sets cannot exceed n_sets.", call. = FALSE)
  }
  if (n_spots < 50L) stop("n_spots must be >= 50.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  m <- ceiling(sqrt(n_spots))
  grid <- expand.grid(y = seq_len(m), x = seq_len(m))
  grid <- grid[seq_len(n_spots), ]
  coords <- data.frame(x = grid$x + stats::runif(n_spots, -0.2, 0.2),
                       y = grid$y + stats::runif(n_spots, -0.2, 0.2))
  block <- if (n_blocks == 1L) NULL else {
    as.factor(cut(coords$x, breaks = n_blocks,
                  labels = paste0("block", seq_len(n_blocks))))
  }

  n_null_sets <- n_sets - n_spatial_sets
  n_spatial_genes <- n_spatial_sets * set_size
  n_null_genes_needed <- n_null_sets * set_size

  # spatial signal per spatial set
  sig <- vector("list", n_spatial_sets)
  types <- if (pattern == "mix") {
    rep(c("bump", "gradient"), length.out = n_spatial_sets)
  } else rep(pattern, n_spatial_sets)
  centers <- data.frame(cx = stats::runif(n_spatial_sets, 1, m),
                        cy = stats::runif(n_spatial_sets, 1, m))
  for (k in seq_len(n_spatial_sets)) {
    if (types[k] == "bump") {
      sig[[k]] <- exp(-((coords$x - centers$cx[k])^2 +
                          (coords$y - centers$cy[k])^2) / (2 * (0.15 * m)^2))
    } else {
      sig[[k]] <- (coords$x - min(coords$x)) /
        max(diff(range(coords$x))) - 0.5
    }
  }

  # gene expression (log space): spatial genes followed by null genes
  logexpr <- matrix(stats::rnorm(n_genes * n_spots, 1.5, noise),
                    nrow = n_genes, ncol = n_spots)
  for (k in seq_len(n_spatial_sets)) {
    gi <- ((k - 1L) * set_size + 1L):(k * set_size)
    a <- stats::runif(set_size, 0.5, 1.5)
    logexpr[gi, ] <- logexpr[gi, ] + effect * outer(a, sig[[k]])
  }
  storage.mode(logexpr) <- "double"
  counts <- matrix(stats::rpois(n_genes * n_spots, exp(logexpr)),
                   nrow = n_genes, ncol = n_spots)
  gnames <- paste0("gene", seq_len(n_genes))
  rownames(counts) <- gnames
  colnames(counts) <- paste0("spot", seq_len(n_spots))

  # gene sets: spatial sets use dedicated spatial genes; null sets use
  # dedicated null genes (disjoint)
  # R >= 4.3.1: paste0(const, zero-length) returns length-1, so guard empties
  spat_names <- if (n_spatial_sets > 0L) {
    paste0("SPATIAL_", seq_len(n_spatial_sets))
  } else character(0)
  null_names <- if (n_null_sets > 0L) {
    paste0("NULL_", seq_len(n_null_sets))
  } else character(0)
  src <- c(spat_names, null_names)
  spatial_gene_idx <- if (n_spatial_sets > 0L) {
    unlist(lapply(seq_len(n_spatial_sets),
                  function(k) ((k - 1L) * set_size + 1L):(k * set_size)))
  } else integer(0)
  null_gene_idx <- setdiff(seq_len(n_genes), spatial_gene_idx)
  null_gene_idx <- null_gene_idx[seq_len(n_null_genes_needed)]

  genes_by_set <- c(
    lapply(seq_len(n_spatial_sets), function(k) {
      gnames[((k - 1L) * set_size + 1L):(k * set_size)]
    }),
    lapply(seq_len(n_null_sets), function(j) {
      gnames[null_gene_idx[((j - 1L) * set_size + 1L):(j * set_size)]]
    }))
  gene_sets <- do.call(rbind, Map(function(s, g) {
    data.frame(source = s, target = g, mor = 1)
  }, src, genes_by_set))

  truth <- data.frame(
    set = src,
    spatial = c(rep(TRUE, n_spatial_sets), rep(FALSE, n_null_sets)),
    pattern = c(types, rep("none", n_null_sets)),
    center_x = c(centers$cx, rep(NA_real_, n_null_sets)) * coord_scale,
    center_y = c(centers$cy, rep(NA_real_, n_null_sets)) * coord_scale)

  coords$x <- coords$x * coord_scale
  coords$y <- coords$y * coord_scale

  list(counts = counts, coords = coords, gene_sets = gene_sets,
       block = block, truth = truth)
}
