# Comparator adapters and reproducible count-data simulations.
# SPARK-X: official SPARK package, mixture kernels and unadjusted gene P values.
# STpathway: the energy::dcov.test call used by the authors' published wrapper,
# with user-supplied gene sets replacing GO discovery. See benchmark_protocol.md.

benchmark_simulate <- function(scenario, seed, n_spots = 225L, n_genes = 1500L) {
  stopifnot(scenario %in% c("global_null", "structured_background", "coherent_weak",
                          "coherent_strong", "opposing_members", "correlated_nonspatial"))
  set.seed(seed)
  side <- ceiling(sqrt(n_spots))
  xy <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))[seq_len(n_spots), ]
  xy <- xy + matrix(runif(n_spots * 2, -0.15, 0.15), ncol = 2)
  genes <- sprintf("g%04d", seq_len(n_genes))
  spots <- sprintf("s%04d", seq_len(n_spots))
  rownames(xy) <- spots
  sizes <- rep(c(12L, 24L, 48L), 4)
  ids <- sample(genes, sum(sizes), replace = FALSE)
  sets <- split(ids, rep(sprintf("set%02d", seq_along(sizes)), sizes))
  baseline <- runif(n_genes, log(1), log(15))
  depth <- rnorm(n_spots, 0, 0.2)
  eta <- outer(baseline, rep(1, n_spots)) + outer(rep(1, n_genes), depth)
  bg <- as.numeric(scale(sin(xy[, 1] / 3) + cos(xy[, 2] / 4)))
  if (scenario %in% c("structured_background", "coherent_weak", "coherent_strong", "opposing_members")) {
    eta <- eta + outer(rnorm(n_genes, 0, 0.35), bg)
  }
  planted <- rep(FALSE, length(sets))
  names(planted) <- names(sets)
  if (scenario %in% c("coherent_weak", "coherent_strong", "opposing_members")) {
    planted[1:3] <- TRUE
    amplitude <- if (scenario == "coherent_weak") 0.35 else 0.8
    for (j in 1:3) {
      center <- c(4, 8, 12)[j]
      signal <- as.numeric(scale(exp(-((xy[, 1] - center)^2 + (xy[, 2] - 8)^2) / 12)))
      member <- match(sets[[j]], genes)
      weight <- if (scenario == "opposing_members") rep(c(-1, 1), length.out = length(member)) else rep(1, length(member))
      eta[member, ] <- eta[member, ] + amplitude * outer(weight, signal)
    }
  }
  if (scenario == "correlated_nonspatial") {
    for (j in seq_along(sets)) {
      member <- match(sets[[j]], genes)
      eta[member, ] <- eta[member, ] + outer(rep(0.8, length(member)), rnorm(n_spots))
    }
  }
  counts <- matrix(rnbinom(length(eta), mu = exp(eta), size = 10), nrow = n_genes,
                   dimnames = list(genes, spots))
  list(counts = counts, coords = xy, sets = sets, planted = planted,
       scenario = scenario, seed = seed)
}

benchmark_run <- function(dat, n_reference = 999L, seed = 1L, stpathway = TRUE,
                          on_method = NULL, progress = FALSE) {
  stopifnot(requireNamespace("SPARK", quietly = TRUE),
            requireNamespace("energy", quietly = TRUE),
            requireNamespace("fgsea", quietly = TRUE))
  counts <- as.matrix(dat$counts)
  coords <- as.matrix(dat$coords)
  stopifnot(identical(colnames(counts), rownames(coords)), all(is.finite(counts)),
            all(counts >= 0), all(colSums(counts) > 0))
  sets <- lapply(dat$sets, function(s) intersect(unique(s), rownames(counts)))
  stopifnot(all(lengths(sets) >= 10L))
  results <- list()
  add <- function(method, endpoint, pathway, p, started) {
    stopifnot(length(pathway) == length(sets), setequal(pathway, names(sets)),
              all(is.finite(p)), all(p >= 0 & p <= 1))
    results[[method]] <<- data.frame(method = method, endpoint = endpoint,
      pathway = pathway, p_value = unname(p), q_value = p.adjust(p, "BH"),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    if (is.function(on_method)) on_method(results[[method]])
    if (progress) message("Completed ", method, " at ", Sys.time())
  }
  tic <- Sys.time()
  sx <- SPARK::sparkx(counts, coords, numCores = 1, option = "mixture", verbose = FALSE)
  gene_p <- setNames(sx$res_mtest$combinedPval, rownames(sx$res_mtest))
  stopifnot(setequal(names(gene_p), rownames(counts)), all(is.finite(gene_p)))
  acat <- vapply(sets, function(s) SPARK::ACAT(gene_p[s]), numeric(1))
  add("SPARKX_ACAT", "absolute", names(acat), acat, tic)
  tic <- Sys.time()
  ranks <- -log10(pmax(gene_p, .Machine$double.xmin))
  set.seed(seed + 31L)
  fg <- fgsea::fgseaMultilevel(pathways = sets, stats = ranks, minSize = 10,
                               maxSize = max(lengths(sets)), scoreType = "pos",
                               eps = 0, nproc = 1)
  add("SPARKX_fgsea", "competitive_gene_evidence", fg$pathway, fg$pval, tic)
  if (stpathway) {
    tic <- Sys.time()
    logexpr <- log1p(sweep(counts, 2, colSums(counts), "/") * 1e4)
    dc <- vapply(seq_along(sets), function(j) {
      if (progress) message("STpathway dCov ", j, "/", length(sets),
                            " ", names(sets)[j], " at ", Sys.time())
      set.seed(seed + 71L)
      energy::dcov.test(t(logexpr[sets[[j]], , drop = FALSE]), coords,
                       R = n_reference)$p.value
    }, numeric(1))
    add("STpathway_dCov", "absolute", names(sets), dc, tic)
  }
  tic <- Sys.time()
  mn <- spatialSVP::matched_null(sets, counts, n_match = n_reference,
                                 min_size = 10, seed = seed)
  match_seconds <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  for (rz in c("none", "pcs")) {
    tic <- Sys.time() - match_seconds
    act <- spatialSVP::compute_pathway_activity(counts, mn$network, coords = coords,
             method = "mean_z", min_size = 10, residualize = rz)
    en <- spatialSVP::matched_enrichment(act, mn, method = "morans_i")
    add(paste0("spatialSVP_matched_", rz), "competitive_score", en$pathway, en$p_matched, tic)
    tic <- Sys.time()
    act$activity <- act$activity[, names(sets), drop = FALSE]
    pp <- spatialSVP::permutation_fdr(act, method = "morans_i",
                                     n_perm = n_reference, seed = seed)
    add(paste0("spatialSVP_permutation_", rz), "absolute_score", pp$pathway, pp$p_perm, tic)
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out$n_genes <- lengths(sets)[out$pathway]
  out
}
