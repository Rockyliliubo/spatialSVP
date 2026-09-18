# Re-evaluate signed PROGENy footprints with weight-preserving references.
suppressPackageStartupMessages(library(spatialSVP))
source("analysis/section_inputs.R")
source("analysis/matched_moran_stream.R")
root <- "analysis/results/progeny_weight_preserving_v1"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
manifest_path <- "analysis/results/cohort_residualization_v1/manifest.rds"
manifest <- if (file.exists(manifest_path)) readRDS(manifest_path) else svp_cohort_manifest()
row <- manifest[manifest$section_id == "GSM8452847_Pt-1A", ]
dat <- svp_load_section(row$h5, row$positions)
network_path <- file.path(root, "progeny_network.rds")
if (file.exists(network_path)) {
  network <- readRDS(network_path)
} else {
  network <- get_progeny_sets("human", top = 250L)
  saveRDS(network, network_path)
}
stopifnot(length(unique(network$source)) == 14L, any(network$mor < 0))
mn <- matched_null(network, dat$counts, n_match = 1999L, min_size = 10L, seed = 13L)
for (s in mn$originals) {
  original <- mn$network$mor[mn$network$source == s]
  child <- names(mn$parents)[mn$parents == s][1]
  reference <- mn$network$mor[mn$network$source == child]
  if (!identical(sort(original), sort(reference)))
    stop("Install the updated spatialSVP version with signed-weight matching.")
}
res <- list()
for (tier in c("none", "libsize", "pcs")) {
  d <- matched_moran_stream(dat$counts, dat$coords, mn, tier = tier)
  d$tier <- tier
  write.csv(d, file.path(root, paste0("results_", tier, ".csv")), row.names = FALSE)
  res[[tier]] <- d
  message("PROGENy ", tier, ": ", sum(d$p_adj_BH < 0.05), "/", nrow(d))
}
all <- do.call(rbind, res)
write.csv(all, file.path(root, "results_all_tiers.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
