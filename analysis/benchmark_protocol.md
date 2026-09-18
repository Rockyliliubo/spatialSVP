# Method-comparison protocol

The six simulation scenarios, sample sizes and methods were fixed before the
formal simulation. This document also records interpretation and reporting
clarifications made during validation of the outputs. The comparison
separates two questions that are often conflated in spatial pathway analysis:

1. **Absolute association:** is a gene or pathway associated with location?
2. **Competitive association:** is the pathway more spatially organized than comparable genes or gene sets?

Methods are interpreted against the question they test. A higher call count under one null is not treated as evidence that the method is more accurate under the other.

## Methods

- **SPARK-X + ACAT:** raw counts are tested gene by gene with the official `SPARK::sparkx` mixture-kernel implementation. Gene-level p-values within each set are combined with `SPARK::ACAT`. This is an absolute pathway endpoint: one or more member genes may drive the result.
- **SPARK-X + fgsea:** the same SPARK-X gene-level evidence is ranked and aggregated with `fgseaMultilevel`. This is a competitive gene-evidence endpoint.
- **STpathway dCov:** the published Brownian distance-covariance kernel is applied to log-normalized pathway expression and spatial coordinates with `energy::dcov.test`. The comparison uses prespecified Hallmark or simulated sets rather than the authors' GO retrieval wrapper. This is an absolute multivariate endpoint.
- **spatialSVP matched null:** the pathway activity statistic is compared with size- and expression-matched random gene sets. This is a competitive score-level endpoint.
- **spatialSVP permutation null:** spot labels are permuted with shared permutations. The direct benchmark applies BH to the resulting unadjusted pathway p-values, matching the multiplicity rule used for the other methods. The package also returns Westfall–Young single-step maxT p-values, which are used in the separate original real-tissue reference-set experiments. This is an absolute score-level endpoint.

The comparator adapters are in `analysis/benchmark_methods.R`. Official SPARK-X code is used directly. The STpathway adapter reproduces the statistical kernel in the authors' public wrapper; it is not presented as a reimplementation of the complete Python workflow.

## Simulation design

Each of 100 replicates contains 225 spots, 1,500 genes and 12 non-overlapping gene sets (four each of size 12, 24 and 48). Counts follow a negative-binomial model with gene-specific abundance and spot-specific depth. Six scenarios are evaluated:

- `global_null`: expression is independent of coordinates;
- `structured_background`: many genes share a spatial background, but no set is planted as exceptional;
- `coherent_weak` and `coherent_strong`: three sets receive concordant localized effects of two amplitudes;
- `opposing_members`: three sets contain spatially patterned genes with effects in opposite directions, so an aggregate mean can cancel;
- `correlated_nonspatial`: genes within each set share non-spatial variation, testing sensitivity to within-set dependence.

The absolute null is genuinely true only in `global_null` and `correlated_nonspatial`. In scenarios with a shared spatial background, non-planted sets are eligible controls for competitive endpoints but not false positives for absolute endpoints. Power therefore means detection of the three planted sets, whereas false-discovery summaries are reported only where the corresponding null is defined.

The matched-null, coordinate-permutation and distance-covariance tests use
999 references or permutations; fgsea uses its adaptive multilevel algorithm.
Pathway-level p-values are adjusted by Benjamini–Hochberg within each method
and replicate. Any-rejection rates use Wilson intervals over 100 independent
replicates. Other metrics use bootstrap intervals resampling replicates;
an all-zero or all-one outcome can yield a degenerate bootstrap interval.
The tests sharing non-spatial modules probe absence of spatial association,
but the constructed correlation differs between annotated sets and random
references, so this single scenario does not establish general competitive
exchangeability. SPARK-X plus fgsea is a baseline adapter, not a claim about
every gene-level aggregation procedure.

Runtime is descriptive. SPARK-X gene testing is charged to the ACAT row;
the fgsea row measures only its additional aggregation cost. Matched-set
generation is included in each matched spatialSVP row, whereas permutation
rows exclude upstream scoring. These are component timings, not comparable
end-to-end costs.

## Real-section comparison

The direct comparison uses the Pt-1A section and all 50 Hallmark pathways after the same gene and spot filters used by the main analysis. Input transformations follow the native requirements of each method: raw counts for SPARK-X and log counts per 10,000 for the STpathway kernel. The analysis reports method-specific adjusted p-values, call counts and pathway-level concordance. These outputs illustrate how the tested hypotheses differ; they are not treated as ground-truth accuracy labels.
