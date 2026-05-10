# chamaeleon: Multiscale analysis of chameleon behaviour in stochastic strange attractors

Implementation of multiscale times-series analysis combining Takens
embedding, multivariate empirical mode decomposition (MEMD), and extreme
value theory (EVT) for non-stationary dynamical systems. Computes
scale-dependent instantaneous dimension and inverse persistence to
characterize chameleon attractors, a kind of stochastic strange
attractors whose geometric and topological properties vary across time
and scales. Core methodology based on Alberti et al. (2023)
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)
. Provides an original surrogate-based statistical testing framework for
formal hypothesis testing of chameleon behavior against scale-invariant
null models, including multiple test statistics, effect size estimation,
and non-stationary analysis.

Implementation of multiscale dynamical systems analysis combining Takens
embedding, multivariate empirical mode decomposition (MEMD), and extreme
value theory (EVT) based metrics. Computes scale-dependent instantaneous
dimension and inverse persistence to characterize chameleon attractors,
whose geometric and topological properties vary across time and scales.

## Main functions

- [`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md):
  Complete analysis workflow

- [`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md):
  Phase-space reconstruction

- [`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md):
  Multivariate empirical mode decomposition

- [`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md):
  EVT-based dynamical metrics

- [`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md):
  Scale-dependent D(t,f) and theta(t,f)

## Visualization

All main objects have S3 plot methods:

- `plot(embedded)`: 3D phase-space visualization for takens_embedding
  objects

- `plot(memd_obj)`: MEMD decomposition modes, spectra, and
  reconstruction

- `plot(scale_metrics)`: Scale-dependent metrics plots

- `plot(evt_metrics)`: EVT metrics time series and distributions

- `plot(result)`: Multi-panel analysis summary for chameleon_analysis
  objects

- [`plot_trajectory_3d`](https://robustecologies.github.io/chamaeleon/reference/plot_trajectory_3d.md):
  3D trajectory line plot

## Theory

See the THEORY.md file in the package directory for detailed
mathematical documentation, or run
`vignette("theory", package = "chameleon")`.

## References

Alberti T, Daviaud F, Donner RV, Dubrulle B, Faranda D, Lucarini V
(2023). Chameleon attractors in turbulent flows. Chaos, Solitons and
Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

## See also

Useful links:

- <https://github.com/robustecologies/chamaeleon>

- <https://robustecologies.github.io/chamaeleon>

- Report bugs at <https://github.com/robustecologies/chamaeleon/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))

Authors:

- Pablo Almaraz <pablo.almaraz@csic.es>
  ([ORCID](https://orcid.org/0000-0003-1416-2695))
