# Fixed none-primary/libsize-sensitivity extension of the archived pcs cohort.
# Run once with SVP_RESID_INIT=1, then optionally split disjoint workers using
# SVP_N_CHUNKS and SVP_CHUNK_INDEX. Existing outputs are validated before reuse.
suppressPackageStartupMessages(library(spatialSVP))
source("analysis/section_inputs.R")
source("analysis/matched_moran_stream.R")
outdir <- "analysis/results/cohort_residualization_v1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(outdir, "manifest.rds")
network_path <- file.path(outdir, "hallmark_network.rds")
if (Sys.getenv("SVP_RESID_INIT", "0") == "1") {
  stopifnot(!file.exists(manifest_path), !file.exists(network_path))
  saveRDS(svp_cohort_manifest(), manifest_path)
  saveRDS(svp_hallmark(), network_path)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("Initialized fixed 163-section manifest and Hallmark collection.")
  quit(status = 0L)
}
manifest <- readRDS(manifest_path)
hallmark <- readRDS(network_path)
tiers <- strsplit(Sys.getenv("SVP_RESID_TIERS", "none,libsize"), ",", fixed = TRUE)[[1]]
stopifnot(all(tiers %in% c("none", "libsize", "pcs")))
nchunks <- as.integer(Sys.getenv("SVP_N_CHUNKS", "1"))
chunk <- as.integer(Sys.getenv("SVP_CHUNK_INDEX", "1"))
stopifnot(nchunks >= 1L, chunk >= 1L, chunk <= nchunks)
selected <- which((manifest$seed_index - 1L) %% nchunks == chunk - 1L)
if (Sys.getenv("SVP_RESID_PILOT", "0") == "1")
  selected <- which(manifest$section_id == "GSM8452847_Pt-1A")
for (tier in tiers) dir.create(file.path(outdir, tier), showWarnings = FALSE)
for (i in selected) {
  sid <- manifest$section_id[i]
  paths <- setNames(file.path(outdir, tiers, paste0(sid, ".rds")), tiers)
  for (tier in tiers[file.exists(paths)]) {
    old <- readRDS(paths[tier])
    stopifnot(old$design == "residualization_v1", old$seed == manifest$seed[i],
              old$n_match == 1999L, old$tier == tier,
              identical(old$section_id, sid),
              all(old$results$n_matched == 1999L))
  }
  todo <- tiers[!file.exists(paths)]
  if (!length(todo)) next
  started <- Sys.time()
  message("START ", sid, " index ", i, "/163 at ", started)
  dat <- svp_load_section(manifest$h5[i], manifest$positions[i])
  mn <- matched_null(hallmark, dat$counts, n_match = 1999L,
                     min_size = 10L, seed = manifest$seed[i])
  message("MATCHED ", sid, " ", ncol(dat$counts), " spots; ", nrow(dat$counts),
          " genes at ", Sys.time())
  for (tier in todo) {
    res <- matched_moran_stream(dat$counts, dat$coords, mn, tier = tier)
    res$cohort <- manifest$cohort[i]
    res$patient <- manifest$patient[i]
    res$section_id <- sid
    res$n_genes <- nrow(dat$counts)
    res$tier <- tier
    item <- list(design = "residualization_v1", seed = manifest$seed[i],
      n_match = 1999L, tier = tier, section_id = sid, results = res)
    tmp <- paste0(paths[tier], ".partial")
    saveRDS(item, tmp)
    stopifnot(file.rename(tmp, paths[tier]))
    message("DONE ", sid, " ", tier, " ", sum(res$p_adj_BH < 0.05),
            "/", nrow(res), " at ", Sys.time())
  }
  rm(dat, mn, res, item)
  invisible(gc())
  message("SECTION_SECONDS ", sid, " ", as.numeric(difftime(Sys.time(), started, units = "secs")))
}
message("WORKER_COMPLETE ", chunk, "/", nchunks, " at ", Sys.time())
