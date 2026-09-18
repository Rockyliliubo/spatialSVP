# Shared section loading and fixed Hallmark membership for cohort analyses.
svp_hallmark <- function() {
  hm <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
    error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens",
      db = "MSigDB_Hallmark"), error = function(e2)
        msigdbr::msigdbr(species = "Homo sapiens", collection = "h")))
  src <- intersect(c("gs_name", "pathway"), names(hm))[1]
  gene <- intersect(c("gene_symbol", "hgnc_symbol"), names(hm))[1]
  net <- unique(data.frame(source = hm[[src]], target = hm[[gene]], mor = 1))
  stopifnot(length(unique(net$source)) == 50L)
  net
}

svp_load_section <- function(h5, positions) {
  raw <- Seurat::Read10X_h5(h5)
  if (is.list(raw)) raw <- raw[["Gene Expression"]]
  stopifnot(!is.null(raw))
  pos <- read.csv(positions, header = FALSE, stringsAsFactors = FALSE)
  if (pos[1, 1] == "barcode") pos <- pos[-1, , drop = FALSE]
  stopifnot(ncol(pos) >= 6L, !anyDuplicated(pos[, 1]))
  names(pos)[1:6] <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
  pos <- pos[as.integer(pos$in_tissue) == 1L, , drop = FALSE]
  common <- intersect(colnames(raw), pos$barcode)
  stopifnot(length(common) > 0L)
  raw <- raw[, common, drop = FALSE]
  keep <- Matrix::colSums(raw > 0) >= 200L
  raw <- raw[, keep, drop = FALSE]
  stopifnot(ncol(raw) >= 2L)
  raw <- raw[Matrix::rowSums(raw > 0) >= 20L, , drop = FALSE]
  counts <- as.matrix(raw)
  ix <- match(colnames(counts), pos$barcode)
  coords <- data.frame(x = as.numeric(pos$pxl_col[ix]), y = as.numeric(pos$pxl_row[ix]),
                       row.names = colnames(counts))
  stopifnot(all(is.finite(as.matrix(coords))), !anyDuplicated(rownames(counts)),
            identical(colnames(counts), rownames(coords)), all(colSums(counts) > 0))
  list(counts = counts, coords = coords)
}

svp_cohort_manifest <- function() {
  source("analysis/patient_id_utils.R", local = TRUE)
  root82 <- Sys.getenv("SVP_DATA_ROOT_GSE282302", unset = "")
  root45 <- Sys.getenv("SVP_DATA_ROOT_GSE274557", unset = "")
  if (!nzchar(root82) || !dir.exists(root82)) {
    stop("Set SVP_DATA_ROOT_GSE282302 to the downloaded GSE282302 directory.")
  }
  if (!nzchar(root45) || !dir.exists(root45)) {
    stop("Set SVP_DATA_ROOT_GSE274557 to the downloaded GSE274557 directory.")
  }
  h82 <- list.files(root82, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
  h45 <- list.files(root45, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
  h45 <- h45[!grepl("PDX", basename(h45))]
  sid <- function(h) sub("_filtered_feature_bc_matrix\\.h5$", "", basename(h))
  m82 <- data.frame(cohort = "GSE282302", h5 = h82, section_id = sid(h82),
                    patient = parse_gse282302_patient(sid(h82)))
  m45 <- data.frame(cohort = "GSE274557", h5 = h45, section_id = sid(h45),
                    patient = parse_gse274557_patient(sid(h45)))
  m <- rbind(m82, m45)
  m$seed_group <- ifelse(m$cohort == "GSE282302",
    sub("^GSM[0-9]+_([A-Z0-9]+)_.*$", "\\1", m$section_id), m$patient)
  m <- m[order(m$cohort, m$seed_group, m$section_id), ]
  m$seed_index <- seq_len(nrow(m))
  m$seed <- 22000L + m$seed_index
  m$positions <- vapply(seq_len(nrow(m)), function(i) {
    if (m$cohort[i] == "GSE282302") {
      opts <- file.path(root82, paste0(m$section_id[i],
        c("_tissue_positions_list.csv.gz", "_tissue_positions.csv.gz")))
      opts[file.exists(opts)][1]
    } else {
      list.files(file.path(root45, paste0(m$section_id[i], "_spatial")),
        pattern = "tissue_positions.*\\.csv$", recursive = TRUE, full.names = TRUE)[1]
    }
  }, character(1))
  rownames(m) <- NULL
  assert_manifest_counts(m)
  stopifnot(!anyDuplicated(m$section_id), all(file.exists(m$h5)),
            !anyNA(m$positions), all(file.exists(m$positions)))
  m
}
