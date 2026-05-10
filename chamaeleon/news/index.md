# Changelog

## chamaeleon 0.3.0

### CRAN-readiness consolidation

- New S3 trio for `chameleon_surrogates`:
  [`print.chameleon_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_surrogates.md)
  reports method, surrogate count, input type and per-surrogate
  dimensionality;
  [`summary.chameleon_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_surrogates.md)
  returns marginal statistics across the ensemble (per-scale mean and
  standard deviation of dimension and persistence for
  `input_type == "scale_metrics"`; per-element mean and standard
  deviation for `input_type == "series"`). Closes the asymmetry under
  which
  [`plot.chameleon_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  was the only registered method.

- [`takens_embed()`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
  gains a `verbose` argument (default `FALSE`). The informational
  messages emitted when `tau` or `m` are auto-estimated are now printed
  only when `verbose = TRUE`. Existing call sites that supplied both
  arguments are unaffected.

  

## chamaeleon 0.2.2

### API symmetry for `generate_surrogates()`

- [`generate_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  now accepts a numeric time-series vector or a `chameleon_analysis`
  object under the default `method = "scale_shuffle"`. In those cases
  the scale-resolved metrics are extracted internally via
  [`chameleon_analysis()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md),
  mirroring the auto-extraction already performed by
  [`chameleon_test()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md).
  Previously the call errored with
  `"scale_shuffle requires scale_metrics object"`. A new `...` argument
  forwards extra parameters to
  [`chameleon_analysis()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
  when this internal extraction is triggered. The diagnostic for
  `method = "scale_shuffle"` on a numeric matrix has been clarified.

  

## chamaeleon 0.2.1

### Visualization migrated to ggplot2

- `plot.evt_metrics`, `plot.chameleon_analysis`, `plot.scale_metrics`
  (static) and the `plot_2d_projection` fallback no longer use
  `graphics::` calls. All four now return ggplot2 + patchwork
  compositions with mandatory subtitle and caption. Return value is now
  the composed object (previously `invisible(NULL)`).

- [`generate_surrogates()`](https://robustecologies.github.io/chamaeleon/reference/generate_surrogates.md)
  now returns an object of class `"chameleon_surrogates"`. A
  `plot.chameleon_surrogates` method is provided: for time-series
  surrogates it overlays density curves; for scale-metrics surrogates it
  shows a surrogate envelope (5th-95th percentile ribbon) against the
  observed metrics.

  

## chamaeleon 0.2.0

### Algorithmic corrections and validation gates

- [`memd()`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
  now uses the Rehman-Mandic monotonicity criterion (total extrema
  across all channels at most 1) to terminate sifting. The previous
  heuristic flagged a residual as monotonic when each channel had more
  than three extrema, which could halt decomposition before all
  meaningful modes were extracted. Behaviour change: `n_imfs` may differ
  from previous releases on signals with high-frequency components.

- [`chameleon_test()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
  now errors (rather than warns) when permutation-invariant statistics
  (`D_range`, `theta_range`, `D_cv`, `theta_cv`) are paired with
  `surrogate_method = "scale_shuffle"`. Those combinations always yield
  p-value = 1.0 and were silent traps for users; the change makes the
  API contract explicit. Use trend-based statistics (`D_trend`,
  `theta_trend`) or gradient statistics (`D_max_gradient`,
  `theta_max_gradient`) with `scale_shuffle`, or change to `fourier`,
  `iaaft` or `block_bootstrap`.

- New canonical-validation test suite
  (`tests/testthat/test-validation-canonical.R`, `skip_on_cran()`) pins
  the EVT dimension on a synthetic Lorenz-63 trajectory at `q = 0.98` to
  within `[1.5, 2.7]` (Faranda-Lucarini convention; published value
  ~2.05), the MEMD mode count on a synthetic 2-channel three-mode signal
  to `[3, 4]`, the chameleon-test p-value on a stationary sinusoid to
  `> 0.10`, and the Takens reconstruction first-column correlation on
  Lorenz-63 to `> 0.99`.

  

## chamaeleon 0.1.0

Initial CRAN release implementing the chameleon attractor methodology
from Alberti et al. (2023) for detecting scale-dependent dynamical
behavior in non-stationary time series.

### Core functionality

- [`chameleon_analysis()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md):
  Complete workflow for chameleon attractor detection, integrating
  Takens embedding, MEMD decomposition, and EVT metrics into a single
  pipeline. Automatically estimates embedding parameters when not
  provided, handles multivariate phase-space reconstruction, and
  computes scale-dependent instantaneous dimension and persistence.
  Returns an S3 object with comprehensive
  [`print()`](https://rdrr.io/r/base/print.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.

- [`takens_embed()`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md):
  Time-delay embedding for phase-space reconstruction from scalar time
  series. Implements Takens’ theorem to unfold the attractor geometry
  using delay coordinates. Automatic parameter estimation via false
  nearest neighbors (FNN) for embedding dimension and autocorrelation
  function (ACF) crossing for time delay. The FNN algorithm identifies
  the minimum dimension where false neighbors (artifacts of projection)
  become negligible.

- [`memd()`](https://robustecologies.github.io/chamaeleon/reference/memd.md):
  Multivariate empirical mode decomposition following the algorithm of
  Rehman & Mandic (2010). Decomposes the embedded trajectory into
  intrinsic mode functions (IMFs) at different timescales without
  requiring predefined basis functions. Implemented in C++/Rcpp with
  OpenMP parallelization for the sifting process across projection
  directions. Returns IMFs ordered by decreasing mean frequency,
  enabling scale-dependent analysis.

- [`evt_metrics()`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md):
  Extreme value theory metrics via generalized Pareto distribution (GPD)
  fitting. For each reference state on the attractor, computes the local
  dimension d (number of active degrees of freedom) and inverse
  persistence theta (extremal index measuring trajectory stability).
  Uses the Pickands-Balkema-de Haan theorem relating recurrence
  statistics to GPD parameters. Threshold selection via quantile-based
  approach with configurable percentile.

- [`scale_dependent_metrics()`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md):
  Computation of D(t,f) and theta(t,f) across frequency scales using
  MEMD partial sums. For each frequency cutoff, reconstructs the
  filtered trajectory and computes EVT metrics, revealing how effective
  dimensionality and persistence vary with observation scale. This is
  the core diagnostic for chameleon behavior: non-constant metrics
  indicate scale-dependent dynamics.

### Statistical testing (original contribution)

The statistical testing framework is an **original contribution** of
this package, extending the Alberti et al. (2023) methodology with
formal hypothesis testing capabilities not present in the original
paper.

- [`chameleon_test()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md):
  Surrogate-based hypothesis testing for chameleon behavior. Tests the
  null hypothesis that D(f) and theta(f) are constant across scales
  using permutation-based inference. Computes multiple test statistics:

  - Trend statistics (Spearman correlation with frequency)
  - Gradient statistics (maximum absolute derivative)
  - Range statistics (max - min across scales)

  Supports multiple correction methods for multiple comparisons
  (Bonferroni, Holm, Benjamini-Hochberg) and reports effect sizes
  (Cohen’s d) alongside p-values for substantive interpretation.

- [`chameleon_test_nonstationary()`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test_nonstationary.md):
  Sliding window analysis for detecting temporal variation in chameleon
  character. Divides the time series into overlapping windows, computes
  a chameleon index for each window, and applies Mann-Kendall trend
  tests to detect systematic changes over time. Useful for identifying
  regime transitions or non-stationary dynamics where the
  scale-dependent character itself evolves.

### Surrogate methods

Four surrogate generation methods are implemented for constructing null
distributions under the hypothesis of scale-invariant dynamics:

- Fourier phase randomization: Randomizes phases while preserving the
  power spectrum exactly. Destroys nonlinear structure while maintaining
  linear correlations.

- IAAFT (Iterative Amplitude Adjusted Fourier Transform): Preserves both
  the power spectrum and the amplitude distribution of the original
  series. More conservative than simple phase randomization for data
  with non-Gaussian marginals.

- Block bootstrap: Resamples contiguous blocks of the time series,
  preserving local temporal structure within blocks while destroying
  long-range dependencies. Block length selected automatically based on
  series length.

- Scale shuffle: Permutes the association between frequency scales and
  metric values. Tests whether the observed pattern of scale-dependence
  is distinguishable from random assignment of metrics to scales.

### Documentation

- Vignette “Introduction to chameleon attractor theory” provides
  theoretical foundations including Takens embedding theorem, MEMD
  algorithm, and EVT framework. Includes worked examples with synthetic
  dynamical systems (Lorenz, Rossler, Van der Pol, Duffing, Henon)
  demonstrating the complete analysis workflow.

- Vignette “Statistical testing for chameleon behavior” covers
  hypothesis testing methodology, surrogate generation, effect size
  interpretation, and practical guidelines, with a systematic comparison
  across synthetic test cases.

### Technical implementation

- C++17 implementation using RcppArmadillo for efficient linear algebra
  operations in MEMD sifting and distance calculations.

- OpenMP parallelization for computationally intensive operations: MEMD
  sifting across projection directions and surrogate generation for
  statistical tests.

- Automatic core detection with configurable parallelism. Respects the
  `_R_CHECK_LIMIT_CORES_` environment variable for CRAN compliance
  during package checking.

- Memory-efficient distance calculations for EVT metrics using streaming
  algorithms where possible.

### References

Alberti T, Daviaud F, Donner RV, Dubrulle B, Faranda D, Lucarini V
(2023). “Chameleon attractors in turbulent flows.” Chaos, Solitons and
Fractals, 168, 113195. <doi:10.1016/j.chaos.2023.113195>
