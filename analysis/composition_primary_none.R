# Composition sensitivity for the primary total-organization estimand.
# Existing cell-type estimates are reused; pathway scores are evaluated without
# PC regression and in memory-bounded chunks. Workers are split by environment
# variables SVP_N_CHUNKS and SVP_CHUNK_INDEX.

suppressPackageStartupMessages(library(spatialSVP))
source("analysis/section_inputs.R")

outdir <- "analysis/results/composition_primary_none_v1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- "analysis/results/cohort_residualization_v1/manifest.rds"
network_path <- "analysis/results/cohort_residualization_v1/hallmark_network.rds"
manifest <- if (file.exists(manifest_path)) readRDS(manifest_path) else svp_cohort_manifest()
network <- if (file.exists(network_path)) readRDS(network_path) else svp_hallmark()
sections <- c("GSM8452847_Pt-1A", "GSM8452865_Pt-6C", "GSM8452893_Pt-13C",
              "GSM8641021_C1_D8_ROI1_s1", "GSM8641030_C2_D11_ROI1_s1",
              "GSM8641067_C3_D12_ROI1_s1")
manifest <- manifest[match(sections, manifest$section_id), ]
stopifnot(!anyNA(manifest$section_id))

marker_modules <- list(
  acinar = c("PRSS1", "CPA1", "CTRB1", "CTRB2", "REG1A"),
  ductal_tumor = c("KRT19", "EPCAM", "MUC1", "KRT8", "KRT18"),
  fibroblast = c("COL1A1", "DCN", "COL3A1", "LUM", "COL1A2"),
  t_cell = c("CD3D", "CD8A", "IL7R", "TRAC", "CD2"),
  macrophage = c("LST1", "CD68", "AIF1", "TYROBP", "C1QA"),
  endothelial = c("PECAM1", "VWF", "CLDN5", "KDR", "FLT1"),
  b_cell = c("CD79A", "MS4A1", "CD19", "BANK1", "CD79B"),
  endocrine = c("INS", "IAPP", "GCG", "SST", "PPY"))

module_scores <- function(counts) {
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  lcpm <- log1p(t(t(counts) / lib) * 1e4)
  z <- t(scale(t(as.matrix(lcpm))))
  z[!is.finite(z)] <- 0
  ans <- sapply(marker_modules, function(g) {
    g <- intersect(g, rownames(z))
    if (!length(g)) return(rep(0, ncol(z)))
    colMeans(z[g, , drop = FALSE])
  })
  rownames(ans) <- colnames(counts)
  ans
}

score_stat <- function(A, W, S0) {
  out <- rep(NA_real_, ncol(A))
  s <- apply(A, 2L, stats::sd)
  valid <- is.finite(s) & s > 0
  if (any(valid)) {
    out[valid] <- spatialSVP:::.stat_morans(
      scale(A[, valid, drop = FALSE]), W, S0)
  }
  out
}

enrich_from_stat <- function(stat, matched) {
  ans <- lapply(matched$originals, function(pathway) {
    refs <- stat[names(matched$parents)[matched$parents == pathway]]
    refs <- refs[is.finite(refs)]
    obs <- unname(stat[pathway])
    stopifnot(is.finite(obs), length(refs) == matched$n_match)
    med <- median(refs)
    spread <- mad(refs, center = med)
    data.frame(
      pathway = pathway,
      stat = obs,
      stat_matched_med = med,
      stat_matched_mad = spread,
      z_matched = if (spread > 0) (obs - med) / spread else NA_real_,
      n_matched = length(refs),
      p_matched = (1 + sum(refs >= obs)) / (1 + length(refs)))
  })
  ans <- do.call(rbind, ans)
  ans$p_adj_BH <- p.adjust(ans$p_matched, "BH")
  ans
}

matched_moran_multi_adjust <- function(counts, coords, matched,
                                       covariates = list(), chunk_sets = 250L) {
  stopifnot(identical(colnames(counts), rownames(coords)))
  z <- spatialSVP:::.zscore_rows(spatialSVP:::.logcpm(counts))
  xy <- spatialSVP:::.coords_xy(coords)
  W <- spatialSVP:::.adjacency_w(xy, k = 6L)
  S0 <- sum(W)
  stopifnot(S0 > 0)

  qr_cov <- lapply(covariates, function(x) {
    x <- as.matrix(x)
    stopifnot(identical(rownames(x), colnames(counts)), all(is.finite(x)))
    qr(cbind(intercept = 1, x))
  })
  net <- matched$network
  set_names <- unique(as.character(net$source))
  group <- match(as.character(net$source), set_names)
  gene <- match(as.character(net$target), rownames(z))
  stopifnot(!anyNA(gene), all(is.finite(net$mor)))
  denom <- rowsum(abs(net$mor), group = group, reorder = FALSE)[, 1]
  denom <- setNames(denom, unique(group))
  weight <- net$mor / denom[as.character(group)]
  membership <- Matrix::sparseMatrix(
    i = gene, j = group, x = weight,
    dims = c(nrow(z), length(set_names)),
    dimnames = list(rownames(z), set_names))
  stats <- c(list(raw = setNames(rep(NA_real_, length(set_names)), set_names)),
             lapply(qr_cov, function(x) setNames(rep(NA_real_, length(set_names)),
                                                   set_names)))
  names(stats) <- c("raw", names(qr_cov))
  for (start in seq.int(1L, length(set_names), by = chunk_sets)) {
    ix <- start:min(start + chunk_sets - 1L, length(set_names))
    A <- as.matrix(Matrix::crossprod(z, membership[, ix, drop = FALSE]))
    stats$raw[ix] <- score_stat(A, W, S0)
    for (nm in names(qr_cov)) {
      stats[[nm]][ix] <- score_stat(qr.resid(qr_cov[[nm]], A), W, S0)
    }
  }
  lapply(stats, enrich_from_stat, matched = matched)
}

read_matrix_csv <- function(path, rowname_col = "barcode") {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(rowname_col %in% names(x), !anyDuplicated(x[[rowname_col]]))
  rn <- x[[rowname_col]]
  x[[rowname_col]] <- NULL
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  rownames(x) <- rn
  x
}

merge_results <- function(x, sid, analysis) {
  out <- do.call(rbind, lapply(names(x), function(nm) {
    d <- x[[nm]]
    d$arm <- nm
    d
  }))
  out$section_id <- sid
  out$analysis <- analysis
  rownames(out) <- NULL
  out
}

run_section <- function(i) {
  row <- manifest[i, ]
  sid <- row$section_id
  final <- file.path(outdir, paste0(sid, ".rds"))
  if (file.exists(final)) {
    message(sid, " checkpoint exists")
    return(invisible(NULL))
  }
  dat <- svp_load_section(row$h5, row$positions)
  props <- read_matrix_csv(file.path("analysis/results/deconvolution",
                                     paste0(sid, "_proportions.csv")))
  common <- colnames(dat$counts)
  stopifnot(all(common %in% rownames(props)))
  props <- props[common, , drop = FALSE]
  mods <- module_scores(dat$counts)

  mn <- matched_null(network, dat$counts, n_match = 1999L, min_size = 10L,
                     seed = row$seed)
  full <- matched_moran_multi_adjust(
    dat$counts, dat$coords, mn,
    covariates = list(marker = mods, label_transfer = props))

  cohort_checkpoint <- readRDS(file.path(
    "analysis/results/cohort_residualization_v1/none", paste0(sid, ".rds")))$results
  check <- merge(full$raw, cohort_checkpoint, by = "pathway",
                 suffixes = c(".new", ".cohort"))
  stopifnot(nrow(check) == nrow(cohort_checkpoint),
            max(abs(check$stat.new - check$stat.cohort)) < 1e-12,
            identical(check$p_matched.new, check$p_matched.cohort))
  result <- merge_results(full, sid, "full_section")

  weights <- read_matrix_csv(file.path("analysis/results/RCTD",
                                       paste0(sid, "_weights.csv")))
  kept <- intersect(colnames(dat$counts), rownames(weights))
  dsub <- list(counts = dat$counts[, kept, drop = FALSE],
               coords = dat$coords[kept, , drop = FALSE])
  weights <- weights[kept, , drop = FALSE]
  mn_rctd <- matched_null(network, dsub$counts, n_match = 1999L,
                          min_size = 10L, seed = row$seed)
  rctd <- matched_moran_multi_adjust(
    dsub$counts, dsub$coords, mn_rctd, covariates = list(RCTD = weights))
  result <- rbind(result, merge_results(rctd, sid, "RCTD_intersection"))

  malignant_col <- grep("malignant", colnames(props), ignore.case = TRUE,
                        value = TRUE)
  malignant_col <- malignant_col[!grepl("non", malignant_col,
                                        ignore.case = TRUE)]
  stopifnot(length(malignant_col) == 1L)
  keep_label <- props[, malignant_col] >= 0.5
  if (sum(keep_label) >= 500L) {
    mn_label <- matched_null(network, dat$counts[, keep_label, drop = FALSE],
                             n_match = 1999L, min_size = 10L,
                             seed = row$seed + 1000L)
    label_subset <- matched_moran_multi_adjust(
      dat$counts[, keep_label, drop = FALSE],
      dat$coords[keep_label, , drop = FALSE], mn_label)
    result <- rbind(result, merge_results(label_subset, sid,
                                          "label_malignant_subset"))
  }

  rctd_malignant <- grep("malignant", colnames(weights), ignore.case = TRUE,
                         value = TRUE)
  rctd_malignant <- rctd_malignant[!grepl("non", rctd_malignant,
                                          ignore.case = TRUE)]
  stopifnot(length(rctd_malignant) == 1L)
  keep_rctd <- weights[, rctd_malignant] >= 0.5
  if (sum(keep_rctd) >= 500L) {
    mn_rm <- matched_null(network, dsub$counts[, keep_rctd, drop = FALSE],
                          n_match = 1999L, min_size = 10L,
                          seed = row$seed + 2000L)
    rctd_subset <- matched_moran_multi_adjust(
      dsub$counts[, keep_rctd, drop = FALSE],
      dsub$coords[keep_rctd, , drop = FALSE], mn_rm)
    result <- rbind(result, merge_results(rctd_subset, sid,
                                          "RCTD_malignant_subset"))
  }
  payload <- list(design = "composition_primary_none_v1", section_id = sid,
                  seed = row$seed, n_spots = ncol(dat$counts),
                  n_label_malignant = sum(keep_label),
                  n_rctd_intersection = length(kept),
                  n_rctd_malignant = sum(keep_rctd), results = result)
  tmp <- paste0(final, ".partial")
  saveRDS(payload, tmp)
  file.rename(tmp, final)
  message(sid, " complete")
  invisible(NULL)
}

n_chunks <- as.integer(Sys.getenv("SVP_N_CHUNKS", "1"))
chunk_index <- as.integer(Sys.getenv("SVP_CHUNK_INDEX", "1"))
stopifnot(n_chunks >= 1L, chunk_index >= 1L, chunk_index <= n_chunks)
indices <- which((seq_len(nrow(manifest)) - 1L) %% n_chunks == chunk_index - 1L)
for (i in indices) {
  run_section(i)
  invisible(gc())
}
