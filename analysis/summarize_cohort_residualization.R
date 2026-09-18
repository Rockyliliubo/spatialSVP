# Aggregate only complete tiers and compare with the archived pcs results.
suppressPackageStartupMessages(library(spatialSVP))
root <- "analysis/results/cohort_residualization_v1"
manifest <- readRDS(file.path(root, "manifest.rds"))
old <- read.csv("analysis/results/cohort_WP1_v0.3.0/cohort_long_table.csv",
                 stringsAsFactors = FALSE)
stopifnot(nrow(old) == 8104L, length(unique(old$section_id)) == 163L)
cols <- c("cohort", "patient", "section_id", "pathway", "n_spots", "n_genes",
          "stat", "stat_matched_med", "n_matched", "p_matched", "p_adj_BH")
read_tier <- function(tier) {
  paths <- file.path(root, tier, paste0(manifest$section_id, ".rds"))
  if (!all(file.exists(paths))) stop("Incomplete tier ", tier, ": ",
                                    sum(file.exists(paths)), "/163 sections")
  obj <- lapply(paths, readRDS)
  stopifnot(all(vapply(obj, function(d) identical(d$design, "residualization_v1") &&
    d$n_match == 1999L && identical(d$tier, tier), logical(1))))
  x <- do.call(rbind, lapply(obj, function(d) d$results))
  stopifnot(nrow(x) == 8104L, !anyDuplicated(x[c("section_id", "pathway")]),
    all(x$n_matched == 1999L), all(is.finite(x$p_adj_BH)))
  check <- merge(x, old, by = c("section_id", "pathway"), suffixes = c(".new", ".old"))
  stopifnot(nrow(check) == 8104L,
    identical(check$patient.new, check$patient.old),
    all(check$n_spots.new == check$n_spots.old),
    all(check$n_genes.new == check$n_genes.old))
  write.csv(x, file.path(root, paste0("cohort_", tier, ".csv")), row.names = FALSE)
  x[cols]
}
all <- rbind(transform(read_tier("none"), tier = "none"),
             transform(read_tier("libsize"), tier = "libsize"),
             transform(old[cols], tier = "pcs"))
all$sig <- all$p_adj_BH < 0.05
all$rank_effect <- 1 - all$p_matched
patient <- aggregate(cbind(sig, rank_effect) ~ tier + cohort + patient + pathway,
                     data = all, FUN = mean)
names(patient)[names(patient) == "sig"] <- "fraction_sections"
patient$positive <- patient$fraction_sections >= 0.5
recurrence <- aggregate(positive ~ tier + cohort + pathway, patient, mean)
counts <- aggregate(positive ~ tier + pathway, patient, sum)
names(counts)[3] <- "positive_patients"
robust <- aggregate(positive ~ tier + pathway, recurrence, function(x) all(x >= 0.8))
names(robust)[3] <- "recurrent_80pct_both_cohorts"
counts <- merge(counts, robust)
write.csv(patient, file.path(root, "patient_effects_all_tiers.csv"), row.names = FALSE)
write.csv(recurrence, file.path(root, "recurrence_by_cohort.csv"), row.names = FALSE)
write.csv(counts, file.path(root, "recurrence_summary.csv"), row.names = FALSE)

# Compare clinical tests using exactly the same patient-level rank endpoint.
clin <- read.csv("analysis/results/WP2_clinical_meta_v0.3.0/clinical_table_S1.csv")
stopifnot(nrow(clin) == 13L, !anyDuplicated(clin$patient))
clinical <- merge(patient[patient$cohort == "GSE274557", ], clin, by = "patient")
stopifnot(nrow(clinical) == 3L * 50L * 13L)
groups <- split(clinical, interaction(clinical$tier, clinical$pathway, drop = TRUE))
tests <- do.call(rbind, lapply(groups, function(d) {
  varying <- length(unique(d$rank_effect)) > 1L
  sp <- if (varying) suppressWarnings(cor.test(d$rank_effect, d$os_days,
    method = "spearman", exact = FALSE)) else NULL
  wt <- if (varying) suppressWarnings(wilcox.test(rank_effect ~ treated, data = d)) else NULL
  data.frame(tier = d$tier[1], pathway = d$pathway[1], n_patients = nrow(d),
    rho_os = if (varying) unname(sp$estimate) else NA_real_,
    p_os = if (varying) sp$p.value else 1,
    median_untreated = median(d$rank_effect[!d$treated]),
    median_treated = median(d$rank_effect[d$treated]),
    p_treatment = if (varying) wt$p.value else 1,
    constant_effect = !varying)
}))
tests$q_os <- ave(tests$p_os, tests$tier, FUN = function(p) p.adjust(p, "BH"))
tests$q_treatment <- ave(tests$p_treatment, tests$tier, FUN = function(p) p.adjust(p, "BH"))
write.csv(tests, file.path(root, "clinical_associations_all_tiers.csv"), row.names = FALSE)
write.csv(clinical, file.path(root, "clinical_patient_effects.csv"), row.names = FALSE)
write.csv(all, file.path(root, "cohort_all_tiers.csv"), row.names = FALSE)
summary <- aggregate(recurrent_80pct_both_cohorts ~ tier, counts, sum)
writeLines(c("Complete cohort residualization comparison: 163 sections, 52 patients.",
  "Primary tier: none. Fixed sensitivity tiers: libsize and archived ten-PC pcs.",
  "1,999 matched references; section and patient eligibility unchanged.",
  paste(capture.output(print(summary, row.names = FALSE)), collapse = "\n"),
  "The legacy binomial section-level heterogeneity P values are not used to establish independence of sections.",
  "Clinical associations are exploratory and adjusted separately within each tier."),
  file.path(root, "RUN_SUMMARY.txt"))
print(summary, row.names = FALSE)
