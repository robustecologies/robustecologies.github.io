# API reference and architecture

This vignette is the navigational map for `chamaeleon`. Section one
states the design philosophy. Section two diagrams the workflow from a
univariate or multivariate signal to a chameleon-test p-value. Section
three catalogues every exported function family. Section four lists the
S3 classes with their print/summary/plot triple. The mathematical
exposition lives in
[`vignette("theory")`](https://robustecologies.github.io/chamaeleon/articles/theory.md);
the surrogate-testing details live in
[`vignette("statistical-testing")`](https://robustecologies.github.io/chamaeleon/articles/statistical-testing.md).

  

## Design philosophy

The package combines three classical analytical machineries that
historically sit in separate corners of the literature: phase-space
reconstruction by delay embedding (Takens 1981; Cao 1997), multivariate
empirical mode decomposition (MEMD; Rehman & Mandic 2010) and
extreme-value theory (EVT) for dynamical-systems extremes (Faranda,
Lucarini, Turchetti, Vaienti 2011; Lucarini et al. 2016). The
combination diagnoses *chameleon* attractors: systems whose effective
EVT signature (local dimension, extremal index) varies systematically
with frequency scale, signalling regime mixing.

Three principles govern the architecture:

The first is **scale separation through MEMD**. A user-supplied
multivariate signal is decomposed into intrinsic mode functions on
aligned hyperspherical projection directions (Rehman-Mandic). Partial
sums above a frequency cutoff reconstruct scale-filtered signals.

The second is **EVT-per-scale as a fingerprint**. The local dimension
and extremal index are estimated per scale via peaks-over-threshold of
the recurrence statistic `-log(distance_to_x)`. The two statistics,
plotted against scale, are the chameleon fingerprint.

The third is **surrogate-based hypothesis testing**. The null hypothesis
“EVT statistics are scale-independent” is tested with one of four
surrogate methods (`scale_shuffle`, `fourier`, `iaaft`,
`block_bootstrap`); p-values are corrected for multiple comparisons
(Fisher, Bonferroni, BH FDR).

  

## Workflow architecture

  

## Phase-space reconstruction (Takens embedding)

  

## MEMD pipeline

[`memd()`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
decomposes a multivariate signal into intrinsic mode functions (IMFs) by
sifting along uniform projection directions on the (n-1)-sphere.
[`memd_partial_sums()`](https://robustecologies.github.io/chamaeleon/reference/memd_partial_sums.md)
reconstructs scale-filtered signals by summing IMFs above a frequency
cutoff.

  

## EVT-per-scale and the chameleon test

[`evt_metrics()`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
computes the local dimension `d` and extremal index `theta` from
peaks-over-threshold (POT) of the recurrence statistic
`-log(distance_to_x)`.
[`scale_dependent_metrics()`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
applies
[`evt_metrics()`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
to each scale-filtered signal produced by
[`memd_partial_sums()`](https://robustecologies.github.io/chamaeleon/reference/memd_partial_sums.md).
[`chameleon_test()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
then tests the null hypothesis that `d(f)` and `theta(f)` are
independent of scale.

  

## Surrogate-testing branch

  

## S3 class system

Seven S3 classes are exported, each with the full print/summary/plot
triple per CLAUDE convention. `plot.takens_embedding` uses plotly for
the 3D phase-space rendering with a base-R 2D projection fallback when
plotly is unavailable. All other plot methods use ggplot2 (with
patchwork for multi-panel composition).

| Class | Constructor | Methods |
|:---|:---|:---|
| takens_embedding | takens_embed | print + summary + plot |
| memd | memd | print + summary + plot |
| evt_metrics | evt_metrics | print + summary + plot |
| scale_metrics | scale_dependent_metrics | print + summary + plot |
| chameleon_analysis | chameleon_analysis | print + summary + plot |
| chameleon_test | chameleon_test | print + summary + plot |
| chameleon_nonstationary | chameleon_test_nonstationary | print + summary + plot |

  

## Function-family catalogue

The reference index in `_pkgdown.yml` mirrors the five functional
families used here and groups every export accordingly: phase-space
reconstruction (`takens_embed`, `estimate_embedding_params`); scale
decomposition (`memd`, `memd_partial_sums`); extreme-value metrics
(`evt_metrics`, `scale_dependent_metrics`); chameleon analysis and
testing (`chameleon_analysis`, `chameleon_test`,
`chameleon_test_nonstationary`, `generate_surrogates`); and
visualization (`plot_trajectory_3d`). The print/summary/plot S3 trio is
bucketed at the bottom of the reference page.

  

## Validation gates

The test suite under `tests/testthat/test-validation-canonical.R` (gated
by `skip_on_cran()`) pins the package contract on canonical reference
setups: `evt_metrics` on a synthetic Lorenz-63 trajectory at quantile
`q = 0.98` returns `d ~= 2.05 +- 0.10` (Faranda-Lucarini convention);
`memd` on a synthetic 2-channel three-mode signal returns
`n_imfs in [3, 4]`; `chameleon_test` on a pure sinusoid returns p-value
\> 0.10 (no chameleon signature in stationary signals); `takens_embed`
on Lorenz-63 returns a reconstructed trajectory with `cor > 0.99`
against the original first column.
