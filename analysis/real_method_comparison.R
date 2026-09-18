# Direct comparison on the Pt-1A section using the official SPARK-X package,
# the published STpathway dCov kernel, and spatialSVP. Run from package root.
suppressWarnings(suppressMessages({ library(spatialSVP) }))
source("analysis/benchmark_methods.R")
data_root <- Sys.getenv("SVP_DATA_ROOT", unset = "")
if (!nzchar(data_root) || !dir.exists(data_root)) {
  stop("Set SVP_DATA_ROOT to the downloaded GSE274557 directory.")
}
sample_id <- Sys.getenv("SVP_SAMPLE", "GSM8452847_Pt-1A")
n_reference <- as.integer(Sys.getenv("SVP_REAL_REFERENCE", "999"))
outdir <- file.path("analysis", "results", "direct_method_comparison", sample_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
h5 <- file.path(data_root, paste0(sample_id, "_filtered_feature_bc_matrix.h5"))
spatial_tar <- file.path(data_root, paste0(sample_id, "_spatial.tar.gz"))
stopifnot(file.exists(h5), file.exists(spatial_tar))
counts <- as.matrix(Seurat::Read10X_h5(h5))
td <- file.path(tempdir(), paste0("svp_compare_", sample_id))
dir.create(td, recursive = TRUE, showWarnings = FALSE)
utils::untar(spatial_tar, exdir = td)
pos_file <- list.files(td, pattern = "tissue_positions", recursive = TRUE,
                       full.names = TRUE)[1]
pos <- utils::read.csv(pos_file, header = FALSE, stringsAsFactors = FALSE)
if (pos[1, 1] == "barcode") pos <- pos[-1, ]
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

hm <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens",
                                                 db = "MSigDB_Hallmark"),
    error = function(e2) msigdbr::msigdbr(species = "Homo sapiens",
                                           collection = "h")))
set_col <- intersect(c("gs_name", "pathway"), colnames(hm))[1]
gene_col <- intersect(c("gene_symbol", "hgnc_symbol"), colnames(hm))[1]
sets <- split(hm[[gene_col]], hm[[set_col]])
sets <- lapply(sets, function(s) intersect(unique(s), rownames(counts)))
sets <- sets[lengths(sets) >= 10]
stopifnot(length(sets) == 50L, identical(colnames(counts), rownames(coords)))
message(sample_id, ": ", nrow(counts), " genes x ", ncol(counts),
        " spots; n_reference=", n_reference)
dat <- list(counts = counts, coords = as.matrix(coords), sets = sets,
            planted = setNames(rep(FALSE, length(sets)), names(sets)))
saveRDS(list(sample_id = sample_id, spots = ncol(counts), genes = nrow(counts),
             pathways = lengths(sets), n_reference = n_reference),
        file.path(outdir, "input_manifest.rds"))
warnings_seen <- character()
res <- withCallingHandlers(
  benchmark_run(dat, n_reference = n_reference, seed = 191717L,
    stpathway = TRUE, progress = TRUE,
    on_method = function(d) saveRDS(d, file.path(outdir,
      paste0("checkpoint_", d$method[1], ".rds")))),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    message("Recorded warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
writeLines(unique(warnings_seen), file.path(outdir, "warnings.txt"))
res$q_value <- ave(res$p_value, res$method,
                   FUN = function(p) p.adjust(p, method = "BH"))
utils::write.csv(res, file.path(outdir, "method_results.csv"), row.names = FALSE)

wide <- reshape(res[c("pathway", "method", "p_value", "q_value")],
                idvar = "pathway", timevar = "method", direction = "wide")
utils::write.csv(wide, file.path(outdir, "method_results_wide.csv"), row.names = FALSE)
calls <- aggregate(q_value < 0.05 ~ method + endpoint, res, sum)
names(calls)[3] <- "n_BH_lt_0_05"
calls$n_pathways <- 50L
utils::write.csv(calls, file.path(outdir, "call_counts.csv"), row.names = FALSE)
writeLines(c(sprintf("Direct method comparison: %s", sample_id),
  sprintf("%d genes x %d spots; 50 Hallmark pathways; %d references/permutations",
          nrow(counts), ncol(counts), n_reference),
  "Official SPARK-X uses raw counts and mixture kernels.",
  "STpathway uses its published energy::dcov.test kernel on logCP10k expression; the authors' GO discovery wrapper is not used.",
  "Methods target different null hypotheses; call counts are not an accuracy ranking.", "",
  paste(capture.output(print(calls, row.names = FALSE)), collapse = "\n")),
  file.path(outdir, "RUN_SUMMARY.txt"))
cat(readLines(file.path(outdir, "RUN_SUMMARY.txt")), sep = "\n")
