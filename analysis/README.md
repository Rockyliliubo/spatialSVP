# Reproducibility analyses

The scripts in this directory reproduce the simulation benchmark and the
manuscript-facing analyses. They are run from the package root after installing
the package and its suggested dependencies.

Raw expression matrices are not redistributed. Set these environment variables
to the corresponding downloaded data directories:

- `SVP_DATA_ROOT_GSE274557`: GSE274557 Visium files;
- `SVP_DATA_ROOT_GSE282302`: GSE282302 Visium files;
- `SVP_DLPFC_ROOT`: spatialLIBD HumanPilot DLPFC matrices and coordinates;
- `SVP_SNRNA_REF`: GSE202051 single-nucleus reference, for deconvolution.

The main workflow is:

1. `benchmark_simulation.R` and `summarize_method_benchmark.R` evaluate the six
   simulation scenarios described in `benchmark_protocol.md`.
2. `real_method_comparison.R` runs the same-input Pt-1A comparison.
3. `cohort_residualization.R` evaluates the two pancreatic cohorts; independent
   worker splits can be set with `SVP_N_CHUNKS` and `SVP_CHUNK_INDEX`.
4. `summarize_cohort_residualization.R` constructs patient and cohort summaries.
5. `dlpfc_donor_validation.R` runs the three-donor DLPFC extension.
6. `composition_primary_none.R` reuses saved label-transfer and RCTD estimates
   for the six-section composition sensitivity analysis.
7. `progeny_weight_validation.R` evaluates signed PROGENy footprints with
   weight-preserving matched references.

Large analyses save one result file per replicate or section. Existing files
are validated before they are reused, so interrupted runs can resume with the
same seeds. `figures/make_figure1_publication.py` and
`figures/make_figures2_6_publication.py` generate the six main figures;
Figures 2–6 use the saved result tables. See `figures/README.md` for dependencies.
