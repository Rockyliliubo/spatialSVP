root <- "analysis/results/composition_primary_none_v1"
sections <- c("GSM8452847_Pt-1A", "GSM8452865_Pt-6C", "GSM8452893_Pt-13C",
              "GSM8641021_C1_D8_ROI1_s1", "GSM8641030_C2_D11_ROI1_s1",
              "GSM8641067_C3_D12_ROI1_s1")
paths <- file.path(root, paste0(sections, ".rds"))
stopifnot(all(file.exists(paths)))
objects <- lapply(paths, readRDS)
stopifnot(all(vapply(objects, function(x) {
  identical(x$design, "composition_primary_none_v1") && x$section_id %in% sections
}, logical(1))))
results <- do.call(rbind, lapply(objects, `[[`, "results"))
stopifnot(!anyDuplicated(results[c("section_id", "pathway", "analysis", "arm")]),
          all(results$n_matched == 1999L),
          all(results$p_matched >= 1 / 2000 & results$p_matched <= 1))

pair_adjustment <- function(analysis, adjusted_arm, label) {
  raw <- results[results$analysis == analysis & results$arm == "raw",
                 c("section_id", "pathway", "p_matched", "p_adj_BH")]
  adj <- results[results$analysis == analysis & results$arm == adjusted_arm,
                 c("section_id", "pathway", "p_matched", "p_adj_BH")]
  pair <- merge(raw, adj, by = c("section_id", "pathway"),
                suffixes = c("_raw", "_adjusted"))
  raw_sig <- pair$p_adj_BH_raw < 0.05
  adjusted_sig <- pair$p_adj_BH_adjusted < 0.05
  data.frame(adjustment = label, comparisons = nrow(pair),
             raw_calls = sum(raw_sig), adjusted_calls = sum(adjusted_sig),
             retained_raw_calls = sum(raw_sig & adjusted_sig),
             lost_raw_calls = sum(raw_sig & !adjusted_sig),
             gained_calls = sum(!raw_sig & adjusted_sig))
}

adjustment <- rbind(
  pair_adjustment("full_section", "marker", "marker_modules"),
  pair_adjustment("full_section", "label_transfer", "label_transfer"),
  pair_adjustment("RCTD_intersection", "RCTD", "RCTD"))
adjustment$retention_fraction <- adjustment$retained_raw_calls /
  pmax(adjustment$raw_calls, 1L)

subset_names <- c("label_malignant_subset", "RCTD_malignant_subset")
subset <- do.call(rbind, lapply(subset_names, function(nm) {
  d <- results[results$analysis == nm & results$arm == "raw", ]
  if (!nrow(d)) return(NULL)
  do.call(rbind, lapply(split(d, d$section_id), function(x) {
    data.frame(analysis = nm, section_id = x$section_id[1],
               pathways = nrow(x), calls = sum(x$p_adj_BH < 0.05))
  }))
}))

sample_summary <- do.call(rbind, lapply(objects, function(x) data.frame(
  section_id = x$section_id, n_spots = x$n_spots,
  n_label_malignant = x$n_label_malignant,
  n_rctd_intersection = x$n_rctd_intersection,
  n_rctd_malignant = x$n_rctd_malignant)))

write.csv(results, file.path(root, "results_long.csv"), row.names = FALSE)
write.csv(adjustment, file.path(root, "adjustment_summary.csv"), row.names = FALSE)
write.csv(subset, file.path(root, "malignant_subset_summary.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(root, "sample_summary.csv"), row.names = FALSE)
summary_text <- c(
  "Primary total-organization composition sensitivity: six sections.",
  "All pathway and matched-reference scores use no expression-PC regression.",
  capture.output(print(adjustment, row.names = FALSE)), "",
  "Malignant-enriched subset results:",
  capture.output(print(subset, row.names = FALSE)))
writeLines(summary_text, file.path(root, "RUN_SUMMARY.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")
