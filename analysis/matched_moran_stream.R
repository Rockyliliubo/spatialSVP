# Memory-bounded evaluation of the existing Moran matched-reference statistic.
# Draw all references with matched_null() first, preserving its seed and order.
# Only score storage is chunked; gene normalization, spatial graph, reference
# ranks and multiplicity adjustment are unchanged.
matched_moran_stream <- function(counts, coords, matched, tier = "none",
                                 n_pcs = 10L, chunk_sets = 250L) {
  stopifnot(inherits(matched, "svp_matched"), tier %in% c("none", "libsize", "pcs"),
            identical(colnames(counts), rownames(coords)), chunk_sets >= 1L)
  mat <- spatialSVP:::.logcpm(counts)
  if (tier != "none") mat <- spatialSVP:::.resid_libsize(mat, counts)
  if (tier == "pcs") mat <- spatialSVP:::.resid_pcs(mat, n_pcs)
  z <- spatialSVP:::.zscore_rows(mat)
  rm(mat)
  xy <- spatialSVP:::.coords_xy(coords)
  W <- spatialSVP:::.adjacency_w(xy, k = 6L)
  S0 <- sum(W)
  stopifnot(S0 > 0)
  net <- matched$network
  set_names <- unique(as.character(net$source))
  group <- match(as.character(net$source), set_names)
  gene <- match(as.character(net$target), rownames(z))
  stopifnot(!anyNA(gene), all(is.finite(net$mor)))
  # matched_null has already filtered and combined original network edges,
  # and generated matched members without replacement.
  denom <- rowsum(abs(net$mor), group = group, reorder = FALSE)[, 1]
  denom <- setNames(denom, unique(group))
  weight <- net$mor / denom[as.character(group)]
  M <- Matrix::sparseMatrix(i = gene, j = group, x = weight,
    dims = c(nrow(z), length(set_names)), dimnames = list(rownames(z), set_names))
  stat <- setNames(rep(NA_real_, length(set_names)), set_names)
  for (start in seq.int(1L, length(set_names), by = chunk_sets)) {
    ix <- start:min(start + chunk_sets - 1L, length(set_names))
    A <- as.matrix(Matrix::crossprod(z, M[, ix, drop = FALSE]))
    sd_col <- apply(A, 2L, stats::sd)
    valid <- is.finite(sd_col) & sd_col > 0
    if (any(valid)) stat[ix[valid]] <- spatialSVP:::.stat_morans(
      scale(A[, valid, drop = FALSE]), W, S0)
  }
  output <- lapply(matched$originals, function(s) {
    refs <- stat[names(matched$parents)[matched$parents == s]]
    refs <- refs[is.finite(refs)]
    obs <- unname(stat[s])
    stopifnot(is.finite(obs), length(refs) == matched$n_match)
    med <- median(refs)
    spread <- mad(refs, center = med)
    data.frame(pathway = s, stat = obs, stat_matched_med = med,
      stat_matched_mad = spread,
      z_matched = if (spread > 0) (obs - med) / spread else NA_real_,
      n_matched = length(refs),
      p_matched = (1 + sum(refs >= obs)) / (1 + length(refs)))
  })
  output <- do.call(rbind, output)
  output$p_adj_BH <- p.adjust(output$p_matched, "BH")
  output$n_spots <- ncol(counts)
  output
}
