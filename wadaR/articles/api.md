# API reference and architecture

This vignette is the navigational map for `wadaR`. It diagrams the four
subsystems (basin computation, attractor specification, Wada detection,
bifurcation analysis), enumerates the function families with cross-links
to the reference index, and lists every S3 class with its
print/summary/plot triple. The mathematical foundations of Wada basins
live in the companion vignette
[`vignette("wada-basins")`](https://robustecologies.github.io/wadaR/articles/wada-basins.md).

  

## Design philosophy

`wadaR` separates four concerns that historically were entangled in
single ad-hoc scripts: how a basin matrix is computed, how attractors
are specified, how the Wada property is tested, and how parameter sweeps
generate bifurcation diagrams. Each is its own family of exported
functions; an exported S3 class is the contract between them.

The package follows a **pipeline contract**: a system specification
(`compiled_system`, `forced_damped_pendulum`, `henon_heiles_system`,
`newton_fractal_system`, `multispecies_competition`) produces an
integrator-ready object;
[`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
integrates trajectories from a grid of initial conditions and labels
each cell by its asymptotic attractor; the resulting `wada_basins`
object is consumed by one of three Wada-detection methods
(`wada_grid_method`, `wada_merging_method`, `wada_straddle_method`) or
by
[`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md)
for an information-theoretic boundary diagnostic; bifurcation sweeps
repeat the pipeline across a parameter range; visualization renders
ggplot2 (2D) or plotly (3D) outputs.

C++ kernels accelerate the inner loops (per-cell trajectory integration,
boundary refinement, Hausdorff distance). OpenMP parallelises across
grid cells. R-side helpers compose the user-facing API.

  

## Architecture overview

  

## Interoperability with janos

[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
accepts two mutually exclusive routes for the dynamics. The native route
takes a string of raw C++ in the `cpp_dynamics` argument, splices it
into the basin template under wadaR’s `state[i]` and bare-name parameter
conventions, and compiles through
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html). The
interoperability route takes a
[`janos::model_spec`](https://robustecologies.github.io/janos/reference/model_spec.html)
object in the `model` argument and adapts it transparently: the internal
helper `adapt_model_spec()` calls the exported hook
[`janos::model_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/model_spec_rhs_cpp.html)
(Suggests), which returns the per-state right-hand sides as C++ strings
under janos’s `y[i]` and `parms[i]` conventions, then rewrites those
references to wadaR’s idioms before splicing. The compiled object is
numerically identical along both routes; the only difference is how the
user expressed the equations in the first place.

Only deterministic ODE and discrete-map `model_spec` objects are
accepted: stochastic differential equations, delay differential
equations, jump-diffusion, partial differential equations, piecewise
deterministic Markov processes, reaction-diffusion master equations and
continuous-time Markov chains are rejected because their right-hand
sides do not reduce to a single per-state derivative or successor
expression. Models built with non-trivial `positive_states` clamps are
also rejected, since positivity is enforced inside the integrator loop
and is not transferable through a per-state expression vector. The
argument `dim`, `type`, and `params` of
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
are inferred from the `model_spec` when not supplied; explicit values
still win when given.

  

## Attractor catalogue

The four attractor constructors describe what the integrator should
recognise as a steady-state outcome. They are consumed by
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
and by the basin-classification kernel.

  

## Wada-detection method tree

The three Wada methods each test the same property, with different
empirical strengths.

  

## Bifurcation sweep

[`bifurcation_basins()`](https://robustecologies.github.io/wadaR/reference/bifurcation_basins.md)
re-runs the basin pipeline across a parameter range and aggregates the
per-step Wada/non-Wada label, basin entropy, and W_m proportions into a
single `bifurcation_result` object that supports
[`animate()`](https://robustecologies.github.io/wadaR/reference/animate.md),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html).

  

## 3D visualization stack

For 3D systems (e.g. Henon-Heiles in the 4D phase space sectioned to
3D),
[`compute_basins_3d()`](https://robustecologies.github.io/wadaR/reference/compute_basins_3d.md)
produces a `basin_result_3d` object rendered through plotly. The
[`slice_3d_basins()`](https://robustecologies.github.io/wadaR/reference/slice_3d_basins.md)
helper extracts a 2D slice at a fixed third coordinate;
[`animate_3d_rotation()`](https://robustecologies.github.io/wadaR/reference/animate_3d_rotation.md)
produces a rotating-camera HTML widget.

  

## S3 class system

Twelve S3 classes are exported, each with the full print/summary/plot
triple per CLAUDE convention. Constructors and consumers are
bidirectionally cross-referenced in the roxygen `@seealso` blocks.

| Class | Constructor | Methods |
|:---|:---|:---|
| wada_basins | compute_basins / compute_newton_basins / compute_competition_basins | print + summary + plot |
| wada_grid_result | wada_grid_method | print + summary + plot |
| wada_merging_result | wada_merging_method | print + summary + plot |
| wada_straddle_result | wada_straddle_method | print + summary + plot |
| wada_analysis | detect_wada | print + summary + plot |
| basin_entropy_result | basin_entropy | print + summary + plot |
| merged_basins | merge_basins | print + summary + plot |
| compiled_system | compiled_system | print + summary + plot |
| attractor_detection | detect_attractors | print + summary + plot |
| bifurcation_result | bifurcation_basins | print + summary + plot |
| basin_result_3d | compute_basins_3d | print + summary + plot |
| multispecies_competition | multispecies_competition | print + summary + plot |

  

## Function-family catalogue

The reference index in `_pkgdown.yml` mirrors the eight functional
families used here and groups every export accordingly: attractor
specification, basin computation, Wada detection methods, basin entropy
and boundary tools, canonical and custom systems, bifurcation analysis,
3D visualization, interactive dashboard. The print/summary/plot S3 trio
is bucketed at the bottom of the reference page.

  

## Validation gates

Per the package contract, the test suite ships canonical regression
tests that exercise the Wada methods on systems with known
classifications (Newton fractal for at : Wada; forced damped pendulum at
: Wada; same pendulum at : non-Wada control). These tests live under
`tests/testthat/test-validation-canonical.R` and are gated by
`skip_on_cran()`. They define the empirical contract that any future
change to `wada_grid_method`, `wada_merging_method` or `basin_entropy`
must preserve.
