# Independent-donor extension of the exploratory cross-tissue panel.
# One preselected 0-micrometre DLPFC section per donor from the original
# spatialLIBD metadata: 151507/Br5292, 151669/Br5595, 151673/Br8100.
suppressWarnings(suppressMessages(library(spatialSVP)))
source("analysis/section_inputs.R")
source("analysis/matched_moran_stream.R")
data_root <- Sys.getenv("SVP_DLPFC_ROOT", unset = "")
if (!nzchar(data_root) || !dir.exists(data_root)) {
  stop("Set SVP_DLPFC_ROOT to the downloaded spatialLIBD DLPFC directory.")
}
out_root <- file.path("analysis", "results", "cross_tissue_dlpfc_donors_v1")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  sample = c("151507", "151669", "151673"),
  donor = c("Br5292", "Br5595", "Br8100"),
  position_um = 0L, replicate = 1L, stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 3L, length(unique(manifest$donor)) == 3L)

hm <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens",
                                                 db = "MSigDB_Hallmark"),
    error = function(e2) msigdbr::msigdbr(species = "Homo sapiens",
                                           collection = "h")))
set_col <- intersect(c("gs_name", "pathway"), colnames(hm))[1]
gene_col <- intersect(c("gene_symbol", "hgnc_symbol"), colnames(hm))[1]
hallmark <- data.frame(source = hm[[set_col]], target = hm[[gene_col]], mor = 1)

load_section <- function(sample) {
  h5 <- file.path(data_root, paste0(sample, "_filtered_feature_bc_matrix.h5"))
  pf <- file.path(data_root, paste0(sample, "_tissue_positions_list.txt"))
  stopifnot(file.exists(h5), file.exists(pf))
  counts <- as.matrix(Seurat::Read10X_h5(h5))
  pos <- utils::read.csv(pf, header = FALSE, stringsAsFactors = FALSE)
  colnames(pos)[1:6] <- c("barcode", "in_tissue", "array_row", "array_col",
                          "pxl_row", "pxl_col")
  pos <- pos[as.integer(pos$in_tissue) == 1L, ]
  common <- intersect(colnames(counts), pos$barcode)
  counts <- counts[, common, drop = FALSE]
  coords <- data.frame(x = as.numeric(pos$pxl_col[match(common, pos$barcode)]),
                       y = as.numeric(pos$pxl_row[match(common, pos$barcode)]))
  keep <- colSums(counts > 0) >= 200
  counts <- counts[, keep, drop = FALSE]
  coords <- coords[keep, , drop = FALSE]
  rownames(coords) <- colnames(counts)
  counts <- counts[rowSums(counts > 0) >= 20, , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(coords)))
  list(counts = counts, coords = coords)
}

for (i in seq_len(nrow(manifest))) {
  sample <- manifest$sample[i]
  outdir <- file.path(out_root, sample)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  done <- file.path(outdir, "RUN_SUMMARY.txt")
  if (file.exists(done)) next
  message("START DLPFC ", sample, " donor ", manifest$donor[i], " at ", Sys.time())
  dat <- load_section(sample)
  mn_ctrl <- matched_null(hallmark, dat$counts, n_match = 2L,
                          min_size = 10L, seed = 210L + i)
  mn <- matched_null(hallmark, dat$counts, n_match = 1999L,
                     min_size = 10L, seed = 410L + i)
  tier_summaries <- list()
  for (tier in c("none", "libsize", "pcs")) {
    act <- compute_pathway_activity(dat$counts, mn_ctrl$network,
      coords = dat$coords, method = "mean_z", min_size = 10L,
      residualize = tier)
    abs_res <- permutation_fdr(act, method = "morans_i", n_perm = 999L,
                               seed = 310L + i)
    utils::write.csv(abs_res, file.path(outdir, paste0("results_absolute_", tier, ".csv")),
                     row.names = FALSE)
    is_hm <- abs_res$pathway %in% mn_ctrl$originals
    rel_res <- matched_moran_stream(dat$counts, dat$coords, mn, tier = tier)
    utils::write.csv(rel_res, file.path(outdir, paste0("results_relative_", tier, ".csv")),
                     row.names = FALSE)
    tier_summaries[[tier]] <- data.frame(sample = sample, donor = manifest$donor[i],
    tier = tier,
    spots = ncol(dat$counts), genes = nrow(dat$counts),
    hallmark_sets = length(mn$originals), random_sets = sum(!is_hm),
    absolute_hallmark_WY = sum(abs_res$p_adj_wy[is_hm] < 0.05, na.rm = TRUE),
    absolute_random_WY = sum(abs_res$p_adj_wy[!is_hm] < 0.05, na.rm = TRUE),
    relative_hallmark_BH = sum(rel_res$p_adj_BH < 0.05, na.rm = TRUE))
    message("DONE DLPFC ", sample, " ", tier, " at ", Sys.time())
  }
  summary <- do.call(rbind, tier_summaries)
  utils::write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)
  writeLines(c("Exploratory DLPFC independent-donor control",
    sprintf("sample %s; donor %s; preselected 0-micrometre replicate 1",
            sample, manifest$donor[i]),
    paste(capture.output(print(summary, row.names = FALSE)), collapse = "\n")), done)
  rm(dat, mn, mn_ctrl, act)
  invisible(gc())
}
parts <- lapply(manifest$sample, function(s)
  utils::read.csv(file.path(out_root, s, "summary.csv"), stringsAsFactors = FALSE))
all_summary <- do.call(rbind, parts)
utils::write.csv(all_summary, file.path(out_root, "donor_summary.csv"), row.names = FALSE)
stopifnot(nrow(all_summary) == 9L, length(unique(all_summary$donor)) == 3L)
cat(paste(capture.output(print(all_summary, row.names = FALSE)), collapse = "\n"), "\n")
