# Patient identifier parsing for the two PDAC Visium cohorts.
#
# These functions deliberately fail on unmatched identifiers. Silent fallback
# would merge biological replicates and invalidate patient-level inference.

.parse_patient_id <- function(section_id, pattern, replacement, cohort,
                              strict = TRUE) {
  section_id <- as.character(section_id)
  matched <- grepl(pattern, section_id, perl = TRUE)
  patient <- rep(NA_character_, length(section_id))
  patient[matched] <- sub(pattern, replacement, section_id[matched],
                          perl = TRUE)
  if (strict && any(!matched)) {
    bad <- unique(section_id[!matched])
    stop(sprintf("Unmatched %s section identifier(s): %s", cohort,
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  patient
}

parse_gse282302_patient <- function(section_id, strict = TRUE) {
  .parse_patient_id(
    section_id,
    pattern = "^GSM[0-9]+_(C[0-9]+_D[0-9]+)(?:_|$).*$",
    replacement = "\\1",
    cohort = "GSE282302",
    strict = strict
  )
}

parse_gse274557_patient <- function(section_id, strict = TRUE) {
  .parse_patient_id(
    section_id,
    pattern = "^GSM[0-9]+_Pt-?([0-9]+)(?:[A-Za-z].*|_.*)?$",
    replacement = "Pt\\1",
    cohort = "GSE274557",
    strict = strict
  )
}

parse_patient_by_cohort <- function(cohort, section_id, strict = TRUE) {
  cohort <- as.character(cohort)
  section_id <- as.character(section_id)
  if (length(cohort) != length(section_id)) {
    stop("cohort and section_id must have equal length", call. = FALSE)
  }
  known <- cohort %in% c("GSE282302", "GSE274557")
  if (strict && any(!known)) {
    stop(sprintf("Unknown cohort(s): %s",
                 paste(unique(cohort[!known]), collapse = ", ")),
         call. = FALSE)
  }
  out <- rep(NA_character_, length(section_id))
  ix82 <- cohort == "GSE282302"
  ix45 <- cohort == "GSE274557"
  out[ix82] <- parse_gse282302_patient(section_id[ix82], strict = strict)
  out[ix45] <- parse_gse274557_patient(section_id[ix45], strict = strict)
  out
}

assert_manifest_counts <- function(manifest,
                                   expected_sections = c(GSE282302 = 108L,
                                                         GSE274557 = 55L),
                                   expected_patients = c(GSE282302 = 39L,
                                                         GSE274557 = 13L)) {
  required <- c("cohort", "patient", "section_id")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing)) {
    stop("Manifest is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(manifest[, required]) || any(manifest$patient == "")) {
    stop("Manifest contains unmatched or empty identifiers", call. = FALSE)
  }
  if (anyDuplicated(manifest$section_id)) {
    stop("Manifest contains duplicated section identifiers", call. = FALSE)
  }
  observed_sections <- vapply(names(expected_sections), function(x) {
    length(unique(manifest$section_id[manifest$cohort == x]))
  }, integer(1))
  observed_patients <- vapply(names(expected_patients), function(x) {
    length(unique(manifest$patient[manifest$cohort == x]))
  }, integer(1))
  if (!identical(unname(observed_sections), unname(as.integer(expected_sections)))) {
    stop(sprintf("Section-count mismatch: observed %s; expected %s",
                 paste(observed_sections, collapse = ","),
                 paste(expected_sections, collapse = ",")), call. = FALSE)
  }
  if (!identical(unname(observed_patients), unname(as.integer(expected_patients)))) {
    stop(sprintf("Patient-count mismatch: observed %s; expected %s",
                 paste(observed_patients, collapse = ","),
                 paste(expected_patients, collapse = ",")), call. = FALSE)
  }
  invisible(list(sections = observed_sections, patients = observed_patients))
}
