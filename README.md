# spatialSVP

**Spatially variable pathway activity for spatial transcriptomics.**

Pathway activity scoring tools (decoupleR, GSVA, PaaSc, ...) tell you *how active*
each pathway is in each spot; gene-level spatial variability tools
(SpatialDE, SPARK-X, nnSVG) tell you *which genes* are spatially structured.
**spatialSVP** tests whether a **pathway activity map is more spatially
organized than comparable gene-set maps**. This differs from
combining gene-level spatial evidence and from direct multivariate tests such
as STpathway; the methods answer related but non-identical questions.

Two null hypotheses, on purpose:

- **Absolute null** — "this score vector has no spatial structure at all."
  On structured tissue this can be false for many ordinary gene sets, so the
  absolute test is a shared-permutation **Westfall-Young maxT** companion
  analysis rather than evidence that an annotated pathway is exceptional.
- **Relative (matched) null** — "is this pathway *more* spatially organized
  than random gene sets of the same size and the same per-gene expression
  distribution?" This is the competitive score test: `matched_null()`
  draws the matched random sets, `matched_enrichment()` turns their spatial
  statistics into an empirical one-sided p-value per pathway. No coordinate
  permutations — cost scales with the number of gene sets.

Per-pathway length scales (Gaussian-process back end) and a three-tier
residualization switch (`none / libsize / pcs`) complete the workflow.

## Scope

- STpathway (Brief Bioinform 2024) directly tests dependence between
  multivariate pathway expression and coordinates. spatialSVP instead tests
  a derived activity map against matched gene-set scores.
- SPathDB (NAR 2025) distributes precomputed spatial pathway activity for
  1.69M spots; PaaSc (PLoS Comput Biol 2025) focuses on pathway scoring.
- Running decoupleR on Visium still leaves the user to hand-wire spatial
  testing downstream; that glue is exactly this package.
- Pathway scores from one matrix can be correlated through shared genes.
  Alongside BH-adjusted p-values, spatialSVP provides
  **gene-set-level permutation with Westfall-Young maxT**: one shared,
  block-aware permutation per replicate keeps the empirical cross-pathway
  correlation inside the joint null.

## Installation

```r
# from GitHub
install.packages("remotes")
remotes::install_github("Rockyliliubo/spatialSVP")
```

Requires R >= 4.3. Suggested packages enable extra back ends: `decoupleR`,
`GSVA`, `fgsea`, `Seurat`/`SeuratObject`, `SpatialExperiment`, `msigdbr`.

## Quick start

```r
library(spatialSVP)

# 1. draw size/expression-matched random gene sets (the relative null)
mn <- matched_null(gene_sets, counts, n_match = 50, seed = 1)

# 2. score pathways + matched sets per spot. The primary example keeps total
#    organization; residualize = "libsize" or "pcs" defines separate
#    residual-organization sensitivity analyses.
act <- compute_pathway_activity(counts, mn$network, coords = coords,
                                 method = "mean_z", residualize = "none")

# 3. matched-null enrichment: empirical one-sided p per pathway
enr <- matched_enrichment(act, mn, method = "morans_i")
enr[order(enr$p_matched), ]

# 4. companion analyses
fdr <- permutation_fdr(act, method = "morans_i", n_perm = 999)  # absolute null, WY maxT
res <- test_spatial_variability(act, method = "gp")             # + length scales

# 5. visualize
plot_pathway_spatial(act, "HALLMARK_GLYCOLYSIS")
plot_volcano(enr, p_col = "p_matched", fdr_col = "p_adj_BH")
plot_lengthscale_map(res)
```

Multi-section or multi-patient designs: pass `block` and all kernels become
block-diagonal while permutations shuffle spots only within blocks.

## API overview

| Function | Purpose |
|---|---|
| `svp_input()` | normalize matrix / SpatialExperiment / Seurat input |
| `compute_pathway_activity()` | score pathways per spot (mean_z, decoupler mlm/wmean, fgsea, GSVA); three-tier residualization |
| `matched_null()` / `matched_enrichment()` | size/expression-matched random sets and the relative-null enrichment test |
| `test_spatial_variability()` | Moran's I / multi-scale covariance / Gaussian-process LRT back ends |
| `permutation_fdr()` | shared-permutation Westfall-Young maxT family-wise error control (block-aware) |
| `get_progeny_sets()` | PROGENy footprint gene sets via decoupleR at runtime |
| `plot_pathway_spatial()` / `plot_volcano()` / `plot_lengthscale_map()` | visualization |
| `simulate_svp_data()` | count-data simulation for methods teaching/testing |

## Methods notes

- Tests are one-sided by design: the alternative is *positive* spatial
  autocorrelation at scales at or above the spot spacing; anti-spatial
  (checkerboard) patterns are deliberately not called significant.
- The `gp` back end is a SpatialDE-inspired Gaussian-process LRT over a
  lengthscale x mixing-proportion grid. p-values are permutation-based by
  default (`n_perm = 199`); the analytic mixture-chi-square p remains as an
  explicit diagnostic (`n_perm = 0`) because it is uncalibrated for the grid
  maximum. Reported `lengthscale` / `dominant_scale`
  values are in the units of the supplied coordinates (pixels for raw
  Visium arrays); `dominant_scale` from the covariance back end is a coarse
  descriptor — sharp scale estimates come from gp `lengthscale`.
- Permutation back ends enforce `n_perm >= 19`; we never fabricate analytic
  nulls where none is defensible.
- **Matched-null semantics**: a small `p_matched` means "more spatially
  organized than random sets of the same size and expression level" — not
  "biologically important". The matched sets control for size and per-gene
  mean expression; weighted pathways also preserve the parent's signed weight
  multiset within expression bins. Within-set correlation is not matched.
- **Compositional caveat**: under library-size normalization, spatially
  over-expressed genes claim a larger share of affected spots' libraries, so
  unrelated sets can inherit coherent negative spatial structure. The
  three-tier `residualize` switch removes the library-size driver
  (`"libsize"`) and additionally the top expression PCs (`"pcs"`); the
  stronger tiers can reduce shared background but can also subtract genuine
  biology. In the accompanying benchmark, ten PCs removed most of a strong
  planted pattern. Treat tiers as different estimands and report them
  separately rather than choosing the smallest p-value.

## Key references

- Svensson V, Teichmann SA, Stegle O. SpatialDE. *Nat Methods*
  2018;15(5):343-346. PMID 29553579.
- Zhu J, Sun S, Zhou X. SPARK-X. *Genome Biol* 2021;22(1):184. PMID 34154649.
- Weber LM, Saha A, Datta A, Hansen KD, Hicks SC. nnSVG. *Nat Commun*
  2023;14(1):4059. PMID 37429865.
- Badia-I-Mompel P, et al. decoupleR. *Bioinform Adv* 2022;2(1):vbac016.
  PMID 36699385.
- Schubert M, et al. PROGENy. *Nat Commun* 2018;9(1):20. PMID 29295995.
- Li F, et al. SPathDB. *Nucleic Acids Res* 2025;53(D1):D1205-D1214.
  PMID 39546631.
- Li S, Wang R, Liu S, Li SC. SPIDER (spatially variable ligand-receptor
  interactions). *Nat Commun* 2025;16(1):7784. PMID 40841363.
- Liao X, et al. PaaSc (pathway activity scoring for single-cell and
  spatial data). *PLoS Comput Biol* 2025;21(11):e1013666. PMID 41212919.
- Tian L, Xiao J, Yu T. STpathway. *Brief Bioinform* 2024;25(6):bbae543.
  DOI 10.1093/bib/bbae543.

## License

MIT (c) 2026 Liubo Li.

## Citation

If you use spatialSVP, cite the archived software release:

Li L. spatialSVP 0.3.0: matched-null testing for spatially variable pathway
activity. Zenodo. 2026. https://doi.org/10.5281/zenodo.22831113

The project repository is https://github.com/Rockyliliubo/spatialSVP. The
concept DOI for all versions is https://doi.org/10.5281/zenodo.22831112.
