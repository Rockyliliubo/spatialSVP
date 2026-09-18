# Run from the package root. Replicates are saved independently for deterministic reruns.
# SVP_BENCH_REPS=100; SVP_BENCH_WORKERS=1; SVP_BENCH_WORKER=1.
# Use SVP_BENCH_OUT for a separate smoke-test directory.
suppressPackageStartupMessages(library(spatialSVP))
source("analysis/benchmark_methods.R")
n_rep <- as.integer(Sys.getenv("SVP_BENCH_REPS", "100"))
n_reference <- as.integer(Sys.getenv("SVP_BENCH_REFERENCE", "999"))
workers <- as.integer(Sys.getenv("SVP_BENCH_WORKERS", "1"))
worker <- as.integer(Sys.getenv("SVP_BENCH_WORKER", "1"))
outdir <- Sys.getenv("SVP_BENCH_OUT", "analysis/results/method_benchmark")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
scenarios <- c("global_null", "structured_background", "coherent_weak",
               "coherent_strong", "opposing_members", "correlated_nonspatial")
jobs <- expand.grid(rep = seq_len(n_rep), scenario = scenarios, stringsAsFactors = FALSE)
jobs$job <- seq_len(nrow(jobs))
jobs <- jobs[(jobs$job - 1L) %% workers == worker - 1L, ]
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, ]
  path <- file.path(outdir, sprintf("%s_%03d.rds", j$scenario, j$rep))
  if (file.exists(path)) {
    old <- readRDS(path)
    stopifnot(identical(old$n_reference, n_reference), old$design_version == "1.0")
    next
  }
  seed <- 170000L + match(j$scenario, scenarios) * 1000L + j$rep
  message(sprintf("START %s replicate %d at %s", j$scenario, j$rep, Sys.time()))
  dat <- benchmark_simulate(j$scenario, seed)
  ans <- benchmark_run(dat, n_reference, seed + 100000L)
  ans$scenario <- j$scenario
  ans$rep <- j$rep
  ans$seed <- seed
  ans$planted <- unname(dat$planted[ans$pathway])
  saveRDS(list(design_version = "1.0", n_reference = n_reference, results = ans,
               session = capture.output(sessionInfo())), path)
  message(sprintf("DONE %s replicate %d at %s", j$scenario, j$rep, Sys.time()))
}
