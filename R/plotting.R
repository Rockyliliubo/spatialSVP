#' Plot pathway activity in space
#'
#' Spatial scatter of one pathway's activity per spot.
#'
#' @param activity A \code{svp_activity} object or a spots x pathways matrix
#'   (then \code{coords} is required).
#' @param pathway Pathway (column) name to plot.
#' @param coords Spot coordinates; defaults to the object's stored coords.
#'
#' @return A ggplot object.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' plot_pathway_spatial(act, "SPATIAL_1")
#'
#' @export
plot_pathway_spatial <- function(activity, pathway, coords = NULL) {
  if (inherits(activity, "svp_activity")) {
    if (is.null(coords)) coords <- activity$coords
    A <- activity$activity
  } else {
    A <- as.matrix(activity)
  }
  coords <- .coords_xy(coords)
  if (!pathway %in% colnames(A)) {
    stop("Pathway not found: ", pathway, ". Available: ",
         paste(colnames(A), collapse = ", "), call. = FALSE)
  }
  df <- data.frame(x = coords$x, y = coords$y, value = A[, pathway])
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = value)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "grey95",
                                   high = "#B2182B",
                                   midpoint = stats::median(df$value)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, color = NULL, title = pathway,
                  subtitle = "pathway activity per spot") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Volcano plot of spatial variability results
#'
#' @param results A \code{svp_results} data.frame from
#'   \code{\link{test_spatial_variability}} or
#'   \code{\link{permutation_fdr}}.
#' @param fdr_col Column holding the adjusted p-values to threshold on
#'   (default \code{"p_adj_BH"}; use \code{"p_adj_wy"} for
#'   \code{permutation_fdr} output).
#' @param fdr_cutoff Significance cutoff.
#' @param p_col Column holding raw p-values (default \code{"p_value"}; use
#'   \code{"p_perm"} for \code{permutation_fdr} output).
#' @param label_top Number of top pathways to label.
#'
#' @return A ggplot object.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' res <- test_spatial_variability(act, method = "morans_i",
#'                                 n_perm = 39, seed = 1)
#' plot_volcano(res)
#'
#' @export
plot_volcano <- function(results, fdr_col = "p_adj_BH", fdr_cutoff = 0.05,
                         p_col = "p_value", label_top = 8L) {
  df <- .as_results_df(results)
  for (col in c(p_col, fdr_col, "stat")) {
    if (!col %in% colnames(df)) {
      stop("Column missing from results: ", col, call. = FALSE)
    }
  }
  df$neglog10p <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$significant <- !is.na(df[[fdr_col]]) & df[[fdr_col]] < fdr_cutoff
  df <- df[order(df$neglog10p, decreasing = TRUE), , drop = FALSE]
  df$lab <- ifelse(seq_len(nrow(df)) <= label_top, df$pathway, "")
  xlab <- paste(attr(results, "method") %||% "", "statistic")
  ggplot2::ggplot(df, ggplot2::aes(x = stat, y = neglog10p,
                                   color = significant)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "grey60", `TRUE` = "#B2182B"),
      name = paste0("FDR < ", fdr_cutoff), na.value = "grey60") +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3,
                       vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::labs(x = xlab, y = sprintf("-log10(%s)", p_col)) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Length-scale map of spatially variable pathways
#'
#' Lollipop chart of pathway-level length scales (gp back end) or dominant
#' kernel scales (covariance back end), ordered and coloured by significance.
#' Moran's-I-only results carry no scale and are rejected with a message.
#'
#' @param results A \code{svp_results} data.frame.
#' @param fdr_col Column holding adjusted p-values.
#' @param fdr_cutoff Significance cutoff.
#'
#' @return A ggplot object.
#'
#' @examples
#' sim <- simulate_svp_data(n_spots = 200, n_genes = 80, n_sets = 6,
#'                          n_spatial_sets = 2, set_size = 8, seed = 1)
#' act <- compute_pathway_activity(sim$counts, sim$gene_sets,
#'                                 coords = sim$coords, method = "mean_z")
#' res <- test_spatial_variability(act, method = "gp", seed = 1)
#' plot_lengthscale_map(res)
#'
#' @export
plot_lengthscale_map <- function(results, fdr_col = "p_adj_BH",
                                 fdr_cutoff = 0.05) {
  df <- .as_results_df(results)
  scale_col <- if ("lengthscale" %in% colnames(df) &&
                     any(is.finite(df$lengthscale))) {
    "lengthscale"
  } else if ("dominant_scale" %in% colnames(df) &&
             any(is.finite(df$dominant_scale))) {
    "dominant_scale"
  } else {
    stop("No length-scale column in results: run method = 'gp' (lengthscale) ",
         "or 'covariance' (dominant_scale) first.", call. = FALSE)
  }
  df$scale <- df[[scale_col]]
  df <- df[is.finite(df$scale), , drop = FALSE]
  df <- df[order(df$scale), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df$significant <- !is.na(df[[fdr_col]]) & df[[fdr_col]] < fdr_cutoff
  ggplot2::ggplot(df, ggplot2::aes(x = scale, y = rank,
                                   color = significant)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = scale,
                                       y = rank, yend = rank),
                          color = "grey75", na.rm = TRUE) +
    ggplot2::geom_point(size = 2.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "grey40", `TRUE` = "#B2182B"),
      name = paste0("FDR < ", fdr_cutoff), na.value = "grey40") +
    ggplot2::labs(x = if (identical(scale_col, "lengthscale")) {
      "length scale (coordinate units)"
    } else {
      "dominant kernel scale (coordinate units)"
    }, y = NULL,
      subtitle = "small scale = local patches; large scale = broad gradients") +
    ggplot2::theme_minimal(base_size = 11)
}

.as_results_df <- function(results) {
  df <- as.data.frame(results)
  if (!"pathway" %in% colnames(df) || !"stat" %in% colnames(df)) {
    stop("Results do not look like svp_results (need 'pathway' and 'stat').",
         call. = FALSE)
  }
  df
}
