# Changelog

## SolidAngleR 0.6.0

  

### Breaking changes

- `omega()` removed. The same multivariate-normal orthant integration is
  now exposed as the `"mvn"` branch of the dispatcher; call
  `compute_solid_angle(V, method = "mvn")` instead. The new branch
  returns the scalar normalized solid angle directly; the per-dimension
  geometric mean previously available as `omega(V)$omega_scaled` is
  recovered as `compute_solid_angle(V, method = "mvn")^(1 / ncol(V))`.

  

### New features

- [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
  gains `method = "mvn"`, dispatching to
  [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
  over the positive orthant via the Genz-Bretz quasi-Monte Carlo
  algorithm. Emits the same condition-number and high-dimension warnings
  as the deleted `omega()`.

  

### Documentation

- Complete roxygen rewrite across every exported symbol: explicit
  `@title` in sentence case, expanded `@description` and `@details` with
  LaTeX equations, verified `@references` with DOI/ISBN, runnable
  `@examples`, and bidirectional `@seealso` links between the S3
  constructor
  [`diagnose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
  and its `print`, `summary`, `plot` methods.
- Internal helpers consistently marked `@noRd`; section headers across
  `R/` normalised to the canonical triple-line convention.
- `_pkgdown.yml` reference index reorganised into ten functional groups
  (main interface, low-dimensional formulas, series methods, geometric
  and decomposition methods, Monte Carlo and sampling, cone diagnostics,
  cone structure utilities, geometry primitives, visualisation, S3
  methods).
- Vignettes: heading-spacing convention enforced across all seven
  articles; chunk evaluation in vignette `a.theoretical-analysis.Rmd`
  set uniformly to `eval = FALSE`. 2D plots in vignettes `b`, `c`, `d`,
  `f`, `g`, and the histogram, density, validation, performance, vMF and
  rejection-comparison panels of `e`, converted to ggplot2 with subtitle
  and grey caption per the package plot-default rule. Two diagrammatic
  visualisations in `e` (Givens rotation arrow plot and Lambertian
  ray-tracing side view) are intentionally retained in base graphics,
  where the arrow / normal annotations are clearer than the ggplot
  equivalent.

  

## SolidAngleR 0.5.1

  

### Bug fixes

- [`tridiagonal_series()`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md):
  default convergence tolerance relaxed from `.Machine$double.eps`
  (about 2.2e-16) to `1e-10`, eliminating spurious non-convergence on
  well-posed tridiagonal cones whose signed residuals sit at
  double-precision noise level. Callers that need the previous strict
  behaviour must pass `tol = .Machine$double.eps` explicitly.
- Numerical-stability corrections in the closed-form, series and
  geometric backends (`formulas_2d.R`, `formulas_3d.R`,
  `geometric_methods.R`, `series_hypergeometric.R`,
  `series_tridiagonal.R`) and in the compiled sampling layer
  (`cone_sampling.cpp`, `cone_sampling.R`), removing silent `NaN`
  outputs at near-degenerate configurations.

  

## SolidAngleR 0.5.0

  

### New features

- Class `cone_diagnosis` gains its full S3 triad: new methods
  [`summary.cone_diagnosis()`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md),
  [`print.summary.cone_diagnosis()`](https://robustecologies.github.io/SolidAngleR/reference/print.summary.cone_diagnosis.md)
  and
  [`plot.cone_diagnosis()`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md).
  [`diagnose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
  returns a first-class diagnostic object with the full print, summary
  and plot lifecycle.
- [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
  dispatch path extended: when `method = "auto"` the internal router now
  consults a `cone_diagnosis` object instead of recomputing the spectral
  and tridiagonal tests on every call.

  

## SolidAngleR 0.4.0

  

### Changes

- Package passes `R CMD check --as-cran` with zero errors, zero warnings
  and a single system-level NOTE (non-portable compiler flag inherited
  from the host toolchain).
- Internal utility layer consolidated around the canonical helpers
  `get_colored_symbol()`, `require_package()` and `require_packages()`;
  deprecated internal helpers removed. No exported signature has
  changed.

  

## SolidAngleR 0.3.0

  

### New features

- New export
  [`shinySolidAngleR()`](https://robustecologies.github.io/SolidAngleR/reference/shinySolidAngleR.md):
  launches a `shinydashboard` front-end in the user’s browser with one
  tab per computational backend (closed forms, hypergeometric series,
  tridiagonal series, decomposition, geometric methods, `omega()`
  wrapper, cone sampler, diagnostics).

  

## SolidAngleR 0.2.0

  

### New features

- New export
  [`hypergeometric_series_nd()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md):
  general Ribando series for arbitrary dimension n greater than or equal
  to 4, with recursive weak-composition enumeration of the multi-index
  lattice. Removes the earlier n = 2, 3 restriction of
  [`hypergeometric_series()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
  and integrates transparently into
  [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md).
- New export
  [`plot_cone_3d()`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md):
  interactive plotly rendering of simplicial cones, showing generator
  rays, triangular faces and the spherical cap cut by the cone.
- [`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
  now dispatches to a C++ backend through Rcpp and RcppArmadillo. The
  O(n) Givens-rotation pipeline replaces the pure-R inner loop and
  delivers order-of-magnitude speed-ups on batch sampling. The R call
  signature is unchanged.

  

### Dependencies

- `Rcpp` and `RcppArmadillo` promoted from `Suggests` to `Imports` and
  `LinkingTo`; the package now requires a C++17 toolchain at install
  time.

  

### Bug fixes

- [`hypergeometric_series_nd()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md):
  corrected the multi-index recursion that produced incorrect values for
  n greater than or equal to 4 at intermediate truncation depths.

  

## SolidAngleR 0.1.0

  

First public release. Forty-three exported functions covering the full
stack of solid-angle methods for polyhedral cones in arbitrary
dimension. The closed forms for n = 2 and n = 3 are exposed through
[`solid_angle_2d()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md),
[`solid_angle_2d_inner()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d_inner.md),
[`solid_angle_3d()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md),
[`solid_angle_3d_det()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_det.md),
[`solid_angle_3d_from_rays()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_from_rays.md),
[`solid_angle_polyhedral()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md)
and
[`spherical_triangle_area()`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md).
The Aomoto-Ribando hypergeometric series for positive-definite cones in
n = 2, 3 is available as
[`hypergeometric_series()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md),
and the Fitisone-Zhou tridiagonal specialisation as
[`tridiagonal_series()`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md).
The recursive Brion-Vergne decomposition, which lifts the
positive-definiteness restriction, is exposed through
[`decompose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/decompose_cone.md)
and
[`solid_angle_decomposition()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md).
Mazonka’s geometric formulas for right circular, segmental, surface and
intersecting cones are wrapped by
[`solid_angle_cone()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone.md),
[`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md),
[`solid_angle_conical_surface()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_conical_surface.md),
[`solid_angle_intersecting_cones()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
and
[`plotIntersectingCones()`](https://robustecologies.github.io/SolidAngleR/reference/plotIntersectingCones.md).
Feasibility-domain computations via multivariate-normal integration are
delivered by `omega()`, and the Monte Carlo estimator by
[`solid_angle_monte_carlo()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_monte_carlo.md).
The pure-R uniform-on-cone sampler stack includes
[`generate_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md),
[`generate_hollow_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md),
[`generate_point_on_sphere()`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md),
[`theta_to_omega()`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md),
[`omega_to_theta()`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md),
[`generate_planar_angle_inverse()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md),
[`generate_planar_angle_rejection()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md),
[`rotate_from_canonical()`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md),
[`rejection_cost()`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
and
[`verify_cone_uniformity()`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md).
Matrix and cone utilities are provided as
[`compute_associated_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md),
[`compute_dot_product_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md),
[`is_positive_definite()`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md),
[`is_tridiagonal()`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md),
[`is_linearly_independent()`](https://robustecologies.github.io/SolidAngleR/reference/is_linearly_independent.md),
[`normalize_vectors()`](https://robustecologies.github.io/SolidAngleR/reference/normalize_vectors.md),
[`create_tridiagonal_cone()`](https://robustecologies.github.io/SolidAngleR/reference/create_tridiagonal_cone.md),
[`cross_product_3d()`](https://robustecologies.github.io/SolidAngleR/reference/cross_product_3d.md),
[`angle_between()`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)
and
[`lhuilier_angle()`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md).
Extreme-ray extraction is handled by
[`generate_spanning_trees()`](https://robustecologies.github.io/SolidAngleR/reference/generate_spanning_trees.md).
The main dispatcher is
[`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md),
with batch form
[`compute_solid_angles()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angles.md),
and the diagnostic constructor
[`diagnose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
with its
[`print.cone_diagnosis()`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md)
method.

  

### Pre-release development history

Three correctness fixes predate the first public commit and are recorded
here for provenance. First, the operator-precedence bug
`gamma(1 + x / 2)` was replaced by the correct `gamma((1 + x) / 2)`
throughout `hypergeometric_series_2d()`, `hypergeometric_series_3d()`
and
[`tridiagonal_series()`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md);
the incorrect form produced approximately tenfold errors in solid-angle
computations. Second, every geometric backend was renormalised to return
values in \[0, 1\]:
[`solid_angle_cone()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone.md),
[`solid_angle_polyhedral()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md),
[`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
and
[`solid_angle_intersecting_cones()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
previously returned steradians and are now consistent with the series
and Monte Carlo backends. Third, spurious `names` attributes were
stripped from
[`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
and
[`solid_angle_decomposition()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
return values.
