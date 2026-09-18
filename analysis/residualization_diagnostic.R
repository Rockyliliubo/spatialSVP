# Quantify how fixed residualization choices remove planted pathway patterns.
# Common simulation seeds pair weak and strong signals; no tier is selected
# using the resulting pathway p-values.
suppressPackageStartupMessages(library(spatialSVP))
source("analysis/benchmark_methods.R")
outdir <- "analysis/results/residualization_diagnostic"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Ensure the benchmark used the same numerical implementations as the source.
source_env <- new.env(parent = globalenv())
sys.source("R/utils.R", envir = source_env)
sys.source("R/compute_pathway_activity.R", envir = source_env)
for (fn in c(".resid_libsize", ".resid_pcs", ".pc_scores",
             "compute_pathway_activity", ".score_mean_z")) {
  stopifnot(identical(body(get(fn, source_env)),
                     body(get(fn, asNamespace("spatialSVP")))))
}

output <- list()
for (r in seq_len(20L)) {
  for (scenario in c("coherent_weak", "coherent_strong")) {
    seed <- 820000L + r
    dat <- benchmark_simulate(scenario, seed)
    logexpr <- log1p(sweep(dat$counts, 2, colSums(dat$counts), "/") * 1e4)
    lib <- log10(colSums(dat$counts))
    # Independent QR calculation validates the package's regression geometry.
    qlib <- qr(cbind(1, lib))
    after_lib <- t(qr.resid(qlib, t(logexpr)))
    stopifnot(max(abs(after_lib - spatialSVP:::.resid_libsize(
      logexpr, dat$counts))) < 1e-9)
    centered <- after_lib - rowMeans(after_lib)
    pcs <- svd(centered, nu = 0L, nv = 10L)$v
    after_pc <- t(qr.resid(qr(cbind(1, pcs)), t(after_lib)))
    stopifnot(max(abs(after_pc - spatialSVP:::.resid_pcs(after_lib, 10L))) < 1e-8)
    mn <- matched_null(dat$sets, dat$counts, n_match = 999L,
                       min_size = 10L, seed = seed + 100000L)
    for (tier in c("none", "libsize", "pcs1", "pcs3", "pcs10")) {
      k <- if (grepl("^pcs", tier)) as.integer(sub("pcs", "", tier)) else 0L
      rz <- if (k > 0L) "pcs" else tier
      act <- compute_pathway_activity(dat$counts, mn$network,
        coords = dat$coords, residualize = rz, n_pcs = k, min_size = 10)
      en <- matched_enrichment(act, mn)
      en$replicate <- r
      en$scenario <- scenario
      en$tier <- tier
      en$planted <- unname(dat$planted[en$pathway])
      en$signal_R2_removed <- NA_real_
      en$signal_R2_library <- NA_real_
      for (j in 1:3) {
        center <- c(4, 8, 12)[j]
        signal <- as.numeric(scale(exp(-((dat$coords[, 1] - center)^2 +
          (dat$coords[, 2] - 8)^2) / 12)))
        signal_after_lib <- qr.resid(qlib, signal)
        design <- if (tier == "none") matrix(1, nrow(dat$coords), 1) else
          if (k == 0L) cbind(1, lib) else cbind(1, lib, pcs[, seq_len(k), drop = FALSE])
        ix <- match(names(dat$sets)[j], en$pathway)
        en$signal_R2_removed[ix] <- 1 - sum(qr.resid(qr(design), signal)^2) / sum(signal^2)
        en$signal_R2_library[ix] <- 1 - sum(signal_after_lib^2) / sum(signal^2)
      }
      output[[length(output) + 1L]] <- en
    }
  }
  message("Paired diagnostic replicate ", r, "/20")
}
x <- do.call(rbind, output)
write.csv(x, file.path(outdir, "paired_results.csv"), row.names = FALSE)
pl <- x[x$planted, ]
summary <- aggregate(cbind(detection = as.numeric(pl$p_adj_BH < 0.05),
  signal_R2_removed = pl$signal_R2_removed,
  signal_R2_library = pl$signal_R2_library),
  by = pl[c("scenario", "tier")], FUN = mean)
write.csv(summary, file.path(outdir, "paired_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
