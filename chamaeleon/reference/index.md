# Package index

## Package overview

Top-level documentation entry.

- [`chamaeleon`](https://robustecologies.github.io/chamaeleon/reference/chamaeleon-package.md)
  [`chamaeleon-package`](https://robustecologies.github.io/chamaeleon/reference/chamaeleon-package.md)
  : chamaeleon: Multiscale analysis of chameleon behaviour in stochastic
  strange attractors

## Phase-space reconstruction

Takens delay embedding and the Cao false-nearest-neighbours dimension
estimator.

- [`takens_embed()`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
  : Takens time-delay embedding
- [`estimate_embedding_params()`](https://robustecologies.github.io/chamaeleon/reference/estimate_embedding_params.md)
  : Estimate embedding parameters

## Scale decomposition (MEMD)

Multivariate empirical mode decomposition (Rehman & Mandic 2010) and
partial-sum reconstruction.

- [`memd()`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
  : Multivariate empirical mode decomposition
- [`memd_partial_sums()`](https://robustecologies.github.io/chamaeleon/reference/memd_partial_sums.md)
  : Compute partial sums of MIMFs for scale-dependent analysis

## Extreme-value metrics

Local dimension and extremal index from peaks-over-threshold of the
recurrence statistic.

- [`evt_metrics()`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
  : Compute EVT-based instantaneous dynamical metrics
- [`scale_dependent_metrics()`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
  : Compute scale-dependent EVT metrics

## Chameleon analysis and testing

End-to-end pipeline plus surrogate-based hypothesis tests.

- [`chameleon_analysis()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
  : Complete chameleon attractor analysis
- [`chameleon_test()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
  : Statistical test for chameleon attractor behavior
- [`chameleon_test_nonstationary()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test_nonstationary.md)
  : Test for non-stationary chameleon behavior
- [`plot(`*`<chameleon_surrogates>`*`)`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  [`generate_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  : Plot method for chameleon_surrogates objects

## Visualization

Dedicated 3D trajectory renderer (plotly).

- [`plot_trajectory_3d()`](https://robustecologies.github.io/chamaeleon/reference/plot_trajectory_3d.md)
  : Plot 3D attractor trajectory as line

## S3 methods

Print, summary and plot methods for the seven exported S3 classes.

- [`print(`*`<chameleon_analysis>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_analysis.md)
  : Print method for chameleon_analysis objects
- [`print(`*`<chameleon_nonstationary>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_nonstationary.md)
  : Print method for chameleon_nonstationary objects
- [`print(`*`<chameleon_surrogates>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_surrogates.md)
  : Print method for chameleon_surrogates objects
- [`print(`*`<chameleon_test>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_test.md)
  : Print method for chameleon_test objects
- [`print(`*`<evt_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.evt_metrics.md)
  : Print method for evt_metrics objects
- [`print(`*`<memd>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.memd.md)
  : Print method for memd objects
- [`print(`*`<scale_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.scale_metrics.md)
  : Print method for scale_metrics objects
- [`print(`*`<summary.chameleon_analysis>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.chameleon_analysis.md)
  : Print method for summary.chameleon_analysis objects
- [`print(`*`<summary.chameleon_surrogates>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.chameleon_surrogates.md)
  : Print method for summary.chameleon_surrogates objects
- [`print(`*`<summary.evt_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.evt_metrics.md)
  : Print method for summary.evt_metrics objects
- [`print(`*`<summary.memd>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.memd.md)
  : Print method for summary.memd objects
- [`print(`*`<summary.scale_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.scale_metrics.md)
  : Print method for summary.scale_metrics objects
- [`print(`*`<summary.takens_embedding>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.summary.takens_embedding.md)
  : Print method for summary.takens_embedding objects
- [`print(`*`<takens_embedding>`*`)`](https://robustecologies.github.io/chamaeleon/reference/print.takens_embedding.md)
  : Print method for takens_embedding objects
- [`summary(`*`<chameleon_analysis>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_analysis.md)
  : Summary method for chameleon_analysis objects
- [`summary(`*`<chameleon_nonstationary>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_nonstationary.md)
  : Summary method for chameleon_nonstationary objects
- [`summary(`*`<chameleon_surrogates>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_surrogates.md)
  : Summary method for chameleon_surrogates objects
- [`summary(`*`<chameleon_test>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_test.md)
  : Summary method for chameleon_test objects
- [`summary(`*`<evt_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.evt_metrics.md)
  : Summary method for evt_metrics objects
- [`summary(`*`<memd>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.memd.md)
  : Summary method for memd objects
- [`summary(`*`<scale_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.scale_metrics.md)
  : Summary method for scale_metrics objects
- [`summary(`*`<takens_embedding>`*`)`](https://robustecologies.github.io/chamaeleon/reference/summary.takens_embedding.md)
  : Summary method for takens_embedding objects
- [`plot(`*`<chameleon_surrogates>`*`)`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  [`generate_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  : Plot method for chameleon_surrogates objects
- [`plot(`*`<chameleon_analysis>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_analysis.md)
  : Plot chameleon analysis summary
- [`plot(`*`<chameleon_nonstationary>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_nonstationary.md)
  : Plot method for chameleon_nonstationary objects
- [`plot(`*`<chameleon_test>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_test.md)
  : Plot method for chameleon_test objects
- [`plot(`*`<evt_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.evt_metrics.md)
  : Plot EVT metrics
- [`plot(`*`<memd>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.memd.md)
  : Plot MEMD decomposition
- [`plot(`*`<scale_metrics>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.scale_metrics.md)
  : Plot scale-dependent metrics
- [`plot(`*`<takens_embedding>`*`)`](https://robustecologies.github.io/chamaeleon/reference/plot.takens_embedding.md)
  : Plot 3D attractor from Takens embedding
- [`plot_trajectory_3d()`](https://robustecologies.github.io/chamaeleon/reference/plot_trajectory_3d.md)
  : Plot 3D attractor trajectory as line
