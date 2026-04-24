# Package index

## Package overview

Package-level documentation, scope, and entry points.

- [`SolidAngleR`](https://robustecologies.github.io/SolidAngleR/reference/SolidAngleR-package.md)
  [`SolidAngleR-package`](https://robustecologies.github.io/SolidAngleR/reference/SolidAngleR-package.md)
  : SolidAngleR: Computing normalized solid angles of polyhedral cones

## Main interface

High-level dispatcher and batch wrapper that select the most appropriate
method based on dimension, conditioning, and structure of the cone.

- [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
  : Compute the normalized solid angle of a polyhedral cone
- [`compute_solid_angles()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angles.md)
  : Vectorised solid-angle computation over a list of cones

## Low-dimensional analytical formulas

Closed-form computations in two and three dimensions via Van Oosterom-
Strackee, L’Huilier, and spherical-triangle area formulas.

- [`solid_angle_2d()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
  : Normalized solid angle of a planar cone via signed area
- [`solid_angle_2d_inner()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d_inner.md)
  : Normalized solid angle of a planar cone via inner product
- [`solid_angle_3d()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
  : Normalized solid angle of a three-dimensional simplicial cone
- [`solid_angle_3d_det()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_det.md)
  : Normalized solid angle in R^3 from a generator matrix
- [`solid_angle_3d_from_rays()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_from_rays.md)
  : Calculate solid angle of a 3D cone from extreme rays
- [`spherical_triangle_area()`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md)
  : Area of the spherical triangle on the unit sphere
- [`lhuilier_angle()`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)
  : Spherical excess of a triangle via L'Huilier's theorem

## Series methods

Ribando hypergeometric series for positive-definite cones and the
Fitisone-Zhou tridiagonal reduction for tridiagonal cones.

- [`hypergeometric_series()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
  : Solid angle via Ribando's multivariable hypergeometric series
- [`hypergeometric_series_nd()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md)
  : General n-dimensional Ribando series for solid angle
- [`tridiagonal_series()`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
  : Solid angle via the tridiagonal Fitisone-Zhou series

## Geometric and decomposition methods

Mazonka-style geometric algorithms for conical surfaces and spherical-
cap intersections, plus the general decomposition of arbitrary cones.

- [`solid_angle_cone()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone.md)
  : Normalized solid angle of a right circular cone
- [`solid_angle_polyhedral()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md)
  : Normalized solid angle of a polyhedral cone via fan triangulation
- [`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
  : Normalized solid angle of a circular cone cut by a plane
- [`solid_angle_intersecting_cones()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
  : Normalized solid angle of the intersection of two circular cones
- [`solid_angle_conical_surface()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_conical_surface.md)
  : Normalized solid angle of a conical surface defined by a curve
- [`solid_angle_decomposition()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
  : Solid angle of a simplicial cone via signed decomposition
- [`decompose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/decompose_cone.md)
  : Signed decomposition of a simplicial cone into PD or low-dimensional
  pieces

## Monte Carlo and sampling

Rcpp-accelerated uniform sampling on spheres and spherical caps, Monte
Carlo solid-angle estimation, and rejection-cost diagnostics.

- [`solid_angle_monte_carlo()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_monte_carlo.md)
  : Generic Monte Carlo estimator of a solid angle
- [`generate_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
  : Uniform random point on a spherical cap
- [`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
  : Vectorised uniform sampling on a spherical cap
- [`generate_hollow_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md)
  : Uniform random point on a hollow spherical cap
- [`generate_point_on_sphere()`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md)
  : Uniform random point on the n-dimensional unit sphere
- [`generate_planar_angle_inverse()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)
  : Random planar angle by inverse-transform sampling
- [`generate_planar_angle_rejection()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
  : Random planar angle by one-dimensional rejection sampling
- [`rejection_cost()`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
  : Expected sample cost of naive sphere-rejection sampling on a cap
- [`verify_cone_uniformity()`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
  : Goodness-of-fit tests for cone-sample uniformity

## Cone diagnostics

S3 constructor returning a structured diagnostic object with
conditioning, stability classification, and method-selection guidance.

- [`diagnose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
  : Diagnose a polyhedral cone and recommend a computation method

## Cone structure utilities

Algebraic predicates and constructors used by the dispatcher to decide
which method applies to a given cone.

- [`compute_associated_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md)
  : Associated matrix M_n(C) of a simplicial cone
- [`compute_dot_product_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_dot_product_matrix.md)
  : Pairwise dot-product matrix V^T V
- [`is_positive_definite()`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md)
  : Test whether a symmetric matrix is positive definite
- [`is_linearly_independent()`](https://robustecologies.github.io/SolidAngleR/reference/is_linearly_independent.md)
  : Test whether the columns of a matrix are linearly independent
- [`is_tridiagonal()`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
  : Test whether a matrix has tridiagonal structure
- [`create_tridiagonal_cone()`](https://robustecologies.github.io/SolidAngleR/reference/create_tridiagonal_cone.md)
  : Construct a tridiagonal simplicial cone from consecutive angles
- [`normalize_vectors()`](https://robustecologies.github.io/SolidAngleR/reference/normalize_vectors.md)
  : Column-wise unit normalisation of a matrix

## Geometry primitives

Low-level vector and angular utilities reused across families.

- [`angle_between()`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)
  : Angle between two vectors via the dot product
- [`cross_product_3d()`](https://robustecologies.github.io/SolidAngleR/reference/cross_product_3d.md)
  : Cross product of two three-dimensional vectors
- [`rotate_from_canonical()`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md)
  : Givens rotation from canonical axis to an arbitrary direction
- [`theta_to_omega()`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md)
  : Planar angle to solid-angle fraction in n dimensions
- [`omega_to_theta()`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md)
  : Solid-angle fraction to planar angle in n dimensions
- [`generate_spanning_trees()`](https://robustecologies.github.io/SolidAngleR/reference/generate_spanning_trees.md)
  : Enumerate all directed spanning trees of a complete digraph

## Visualization and interactive tools

3D plotly visualisations of cones and intersecting spherical caps, plus
a Shiny dashboard for exploration.

- [`plot_cone_3d()`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md)
  : Interactive 3D plotly visualisation of a simplicial cone
- [`plotIntersectingCones()`](https://robustecologies.github.io/SolidAngleR/reference/plotIntersectingCones.md)
  : Interactive 3D visualisation of two intersecting cones
- [`shinySolidAngleR()`](https://robustecologies.github.io/SolidAngleR/reference/shinySolidAngleR.md)
  : Launch the SolidAngleR Shiny dashboard

## S3 methods

Print, summary, and plot methods for the cone_diagnosis class. Users
normally do not call these directly; listed here for completeness and
cross-reference.

- [`print(`*`<cone_diagnosis>`*`)`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md)
  : Print a cone diagnosis report
- [`print(`*`<summary.cone_diagnosis>`*`)`](https://robustecologies.github.io/SolidAngleR/reference/print.summary.cone_diagnosis.md)
  : Print an extended cone diagnostic summary
- [`summary(`*`<cone_diagnosis>`*`)`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md)
  : Extended diagnostic summary for a polyhedral cone
- [`plot(`*`<cone_diagnosis>`*`)`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
  : Plot the eigenvalue spectrum of a diagnosed cone
- [`plotIntersectingCones()`](https://robustecologies.github.io/SolidAngleR/reference/plotIntersectingCones.md)
  : Interactive 3D visualisation of two intersecting cones
- [`plot_cone_3d()`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md)
  : Interactive 3D plotly visualisation of a simplicial cone
