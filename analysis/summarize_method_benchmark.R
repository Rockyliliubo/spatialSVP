# Validate and summarize the saved method-benchmark replicates.
# Does not treat structured-background absolute-null calls as false positives.
suppressWarnings(suppressMessages(library(spatialSVP)))

indir <- Sys.getenv("SVP_BENCH_OUT", "analysis/results/method_benchmark")
files <- list.files(indir, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("No benchmark replicate files under ", indir)
obj <- lapply(files, readRDS)
stopifnot(all(vapply(obj, function(z) identical(z$design_version, "1.0"), logical(1))),
          length(unique(vapply(obj, `[[`, integer(1), "n_reference"))) == 1L)
x <- do.call(rbind, lapply(obj, `[[`, "results"))
rownames(x) <- NULL

scenarios <- c("global_null", "structured_background", "coherent_weak",
               "coherent_strong", "opposing_members", "correlated_nonspatial")
methods <- c("SPARKX_ACAT", "SPARKX_fgsea", "STpathway_dCov",
             "spatialSVP_matched_none", "spatialSVP_matched_pcs",
             "spatialSVP_permutation_none", "spatialSVP_permutation_pcs")
key <- x[c("scenario", "rep", "method", "pathway")]
if (anyDuplicated(key)) stop("Duplicated scenario/rep/method/pathway rows")
expected_reps <- as.integer(Sys.getenv("SVP_BENCH_REPS", "100"))
stopifnot(setequal(unique(x$rep), seq_len(expected_reps)))
expected <- expand.grid(scenario = scenarios, rep = seq_len(expected_reps),
                        method = methods, pathway = sprintf("set%02d", 1:12),
                        stringsAsFactors = FALSE)
key$observed <- TRUE
missing <- merge(expected, key, all.x = TRUE, by = names(expected))
missing <- missing[is.na(missing$observed), ]
if (nrow(missing)) stop("Incomplete benchmark: ", nrow(missing), " rows missing")
stopifnot(all(is.finite(x$p_value)), all(x$p_value >= 0 & x$p_value <= 1),
          all(is.finite(x$q_value)), all(x$q_value >= 0 & x$q_value <= 1),
          all(is.finite(x$elapsed_seconds)), all(x$elapsed_seconds >= 0))

split_rep <- split(x, interaction(x$scenario, x$method, x$rep, drop = TRUE))
rep_metrics <- do.call(rbind, lapply(split_rep, function(d) {
  sig <- d$q_value < 0.05
  has_signal <- any(d$planted)
  competitive <- grepl("^competitive", d$endpoint[1])
  null_eligible <- if (d$scenario[1] %in% c("global_null", "correlated_nonspatial")) {
    rep(TRUE, nrow(d))
  } else if (competitive) {
    !d$planted
  } else {
    rep(FALSE, nrow(d))
  }
  data.frame(scenario = d$scenario[1], method = d$method[1], rep = d$rep[1],
    any_rejection = any(sig), rejection_fraction = mean(sig),
    planted_detection = if (has_signal) mean(sig[d$planted]) else NA_real_,
    null_rejection_fraction = if (any(null_eligible)) mean(sig[null_eligible]) else NA_real_,
    relative_fdp = if (has_signal && competitive && sum(sig) > 0)
      sum(sig & !d$planted) / sum(sig) else if (has_signal && competitive) 0 else NA_real_,
    elapsed_seconds = d$elapsed_seconds[1], stringsAsFactors = FALSE)
}))
rownames(rep_metrics) <- NULL

mean_ci <- function(v, seed) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(mean = NA, lo = NA, hi = NA))
  set.seed(seed)
  boot <- replicate(2000, mean(sample(v, length(v), replace = TRUE)))
  c(mean = mean(v), lo = unname(quantile(boot, 0.025)),
    hi = unname(quantile(boot, 0.975)))
}
wilson_ci <- function(v) {
  n <- length(v)
  p <- mean(v)
  z <- qnorm(0.975)
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(mean = p, lo = max(0, center - half), hi = min(1, center + half))
}
groups <- split(rep_metrics, interaction(rep_metrics$scenario, rep_metrics$method,
                                         drop = TRUE))
summary <- do.call(rbind, Map(function(d, i) {
  vals <- lapply(c("any_rejection", "rejection_fraction", "planted_detection",
                   "null_rejection_fraction", "relative_fdp"), function(nm)
                     mean_ci(as.numeric(d[[nm]]), 180000L + i))
  names(vals) <- c("any_rejection", "rejection_fraction", "planted_detection",
                   "null_rejection_fraction", "relative_fdp")
  vals$any_rejection <- wilson_ci(as.numeric(d$any_rejection))
  data.frame(scenario = d$scenario[1], method = d$method[1], n_replicates = nrow(d),
    any_rejection = vals$any_rejection["mean"],
    any_rejection_lo = vals$any_rejection["lo"],
    any_rejection_hi = vals$any_rejection["hi"],
    rejection_fraction = vals$rejection_fraction["mean"],
    rejection_fraction_lo = vals$rejection_fraction["lo"],
    rejection_fraction_hi = vals$rejection_fraction["hi"],
    planted_detection = vals$planted_detection["mean"],
    planted_detection_lo = vals$planted_detection["lo"],
    planted_detection_hi = vals$planted_detection["hi"],
    null_rejection_fraction = vals$null_rejection_fraction["mean"],
    null_rejection_fraction_lo = vals$null_rejection_fraction["lo"],
    null_rejection_fraction_hi = vals$null_rejection_fraction["hi"],
    relative_fdp = vals$relative_fdp["mean"],
    relative_fdp_lo = vals$relative_fdp["lo"],
    relative_fdp_hi = vals$relative_fdp["hi"],
    median_seconds = median(d$elapsed_seconds), stringsAsFactors = FALSE)
}, groups, seq_along(groups)))
rownames(summary) <- NULL
summary <- summary[order(match(summary$scenario, scenarios),
                         match(summary$method, methods)), ]
utils::write.csv(x, file.path(indir, "benchmark_long.csv"), row.names = FALSE)
utils::write.csv(rep_metrics, file.path(indir, "benchmark_replicate_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(summary, file.path(indir, "benchmark_summary.csv"), row.names = FALSE)

global <- summary[summary$scenario == "global_null",
                  c("method", "n_replicates", "any_rejection",
                    "any_rejection_lo", "any_rejection_hi",
                    "rejection_fraction", "rejection_fraction_lo",
                    "rejection_fraction_hi", "median_seconds")]
signal <- summary[summary$scenario %in% c("coherent_weak", "coherent_strong",
                                          "opposing_members"),
                  c("scenario", "method", "planted_detection",
                    "planted_detection_lo", "planted_detection_hi",
                    "null_rejection_fraction", "relative_fdp")]
writeLines(c(
  sprintf("Validated %d replicates and %d method-pathway rows; n_reference=%d.",
          length(files), nrow(x), obj[[1]]$n_reference),
  "Any-rejection intervals use Wilson binomial intervals over independent replicates.",
  "Other intervals bootstrap independent replicates; intervals can degenerate at all-zero/all-one outcomes.",
  "Calls in structured_background are descriptive because its absolute null is false.",
  "", "Global-null calibration:",
  paste(capture.output(print(global, row.names = FALSE, digits = 3)), collapse = "\n"),
  "", "Planted-set detection (absolute methods are not assigned false-discovery rates when shared spatial background makes their null false):",
  paste(capture.output(print(signal, row.names = FALSE, digits = 3)), collapse = "\n")
), file.path(indir, "RUN_SUMMARY.txt"))
cat(readLines(file.path(indir, "RUN_SUMMARY.txt")), sep = "\n")
