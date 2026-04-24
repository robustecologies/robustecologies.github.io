# Technical architecture and API reference

## Overview

The `SolidAngleR` package computes the normalized solid angle measure of
polyhedral (simplicial) cones in \\\mathbb{R}^n\\ for arbitrary
dimension \\n \geq 2\\. The computational core implements five distinct
algorithmic families, each targeting a specific class of cone geometry:
closed-form analytical formulas for low dimensions, multivariable
hypergeometric series for cones with positive definite associated
matrices, a tridiagonal specialization of the hypergeometric series for
cones whose Gram matrix has banded structure, a recursive decomposition
method that handles the fully general case, and multivariate normal
integration via Quasi-Monte Carlo. A vectorized C++ backend (via Rcpp
and RcppArmadillo) accelerates uniform cone sampling in arbitrary
dimensions. This document describes the mathematical foundations,
numerical algorithms, software architecture, and complete API surface of
the package.

Layered architecture of SolidAngleR. Each cluster is a subsystem; arrows
show the flow from user-facing inputs through the dispatcher and
computational backends down to the C++ sampler.

  

## Mathematical framework

  

### The solid angle measure

Given \\n\\ linearly independent vectors \\v_1, \ldots, v_n \in
\mathbb{R}^n\\, the simplicial cone \\C = \text{cone}(v_1, \ldots, v_n)
= \\\sum\_{i=1}^n \lambda_i v_i : \lambda_i \geq 0\\\\ traces a region
on the unit sphere \\S^{n-1}\\. The normalized solid angle measure
\\\omega(C) \in \[0, 1\]\\ is the fraction of \\S^{n-1}\\ covered by
\\C\\:

\\ \omega(C) = \frac{\text{vol}\_{n-1}(C \cap
S^{n-1})}{\text{vol}\_{n-1}(S^{n-1})} \\

For the positive orthant \\C = \\x \in \mathbb{R}^n : x_i \geq 0\\\\, we
have \\\omega(C) = 1/2^n\\. All methods in the package return this
normalized measure.

  

### The associated matrix \\M_n(C)\\

Central to the package is the associated matrix \\M_n(C)\\ (Definition
1.1 in Fitisone & Zhou [\[1\]](#ref1)), constructed from the unit
generators \\\hat{v}\_i = v_i / \\v_i\\\\:

\\ M\_{ij} = \begin{cases} 1 & \text{if } i = j \\ -\|\hat{v}\_i \cdot
\hat{v}\_j\| & \text{if } i \neq j \end{cases} \\

The positive definiteness of \\M_n(C)\\ determines whether the
hypergeometric series converges (Theorem 1.5 in [\[1\]](#ref1)). The
package computes \\M_n(C)\\ via
[`compute_associated_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md)
and tests positive definiteness by eigendecomposition with tolerance
\\\epsilon = 10^{-10}\\.

  

## Algorithmic methods

  

### Closed formulas (dimensions 2 and 3)

**Two dimensions.** For \\v_1, v_2 \in \mathbb{R}^2\\, the solid angle
is the planar angle between the rays, normalized by \\2\pi\\:

\\ \omega = \frac{1}{2\pi} \arctan2\\\left(\|\det\[v_1, v_2\]\|,\\ v_1
\cdot v_2\right) \\

The implementation uses `atan2` for correct quadrant handling and
returns 0 for collinear or opposite vectors. Colinear detection uses a
tolerance of \\10^{-15}\\ on the determinant magnitude.

**Implementation:** `R/formulas_2d.R` exports `solid_angle_2d(v1, v2)`.

**Three dimensions.** For \\v_1, v_2, v_3 \in \mathbb{R}^3\\, the
Euler-Lagrange formula [\[2\]](#ref2) computes the spherical excess:

\\ E = 2\\\arctan2\\\left(\left\|v_1 \cdot (v_2 \times v_3)\right\|,\\
1 + v_1 \cdot v_2 + v_1 \cdot v_3 + v_2 \cdot v_3\right) \\

and \\\omega = E / (4\pi)\\. The numerator is the absolute value of the
scalar triple product; the denominator handles the case where \\1 + \sum
\text{dot products} \< 0\\ through the `atan2` branch cut. If the triple
product is smaller than \\10^{-15}\\, the cone is degenerate and the
function returns 0.

**Implementation:** `R/formulas_3d.R` exports
`solid_angle_3d(v1, v2, v3)`, `solid_angle_3d_det(V)` (matrix
interface), and `spherical_triangle_area(v1, v2, v3)` (returns \\E\\ in
steradians rather than the normalized measure).

  

### Hypergeometric series (Theorem 1.5)

For \\n \geq 4\\ with \\M_n(C)\\ positive definite, Ribando’s
multivariable hypergeometric series [\[3\]](#ref3) computes the solid
angle as a power series in the pairwise dot products \\\alpha\_{ij} =
\hat{v}\_i \cdot \hat{v}\_j\\:

\\ \omega(C) = \frac{\|\det V\|}{(4\pi)^{n/2}} \sum\_{\mathbf{a} \in
\mathbb{N}\_0^N} \left\[ \frac{(-2)^{\|\mathbf{a}\|}}{\mathbf{a}!}
\prod\_{i=1}^n \Gamma\\\left(\frac{1 + \deg_i(\mathbf{a})}{2}\right)
\right\] \boldsymbol{\alpha}^{\mathbf{a}} \\

where \\N = \binom{n}{2}\\ is the number of distinct pairs, \\\mathbf{a}
= (a\_{12}, a\_{13}, \ldots, a\_{n-1,n})\\ is a multi-index of
exponents, \\\|\mathbf{a}\| = \sum a\_{jk}\\ is the total degree,
\\\mathbf{a}! = \prod a\_{jk}!\\, \\\deg_i(\mathbf{a}) = \sum\_{k \neq
i} a\_{ik}\\ is the degree of vertex \\i\\ in the multigraph defined by
\\\mathbf{a}\\, and \\\boldsymbol{\alpha}^{\mathbf{a}} = \prod
\alpha\_{jk}^{a\_{jk}}\\.

The series converges absolutely if and only if \\M_n(C)\\ is positive
definite. Convergence speed depends on the eigenvalue spread: nearly
orthogonal cones (\\\alpha\_{ij} \approx 0\\) converge in a few terms,
while cones with large dot products require many terms even when
\\M_n(C)\\ is positive definite.

**Algorithmic strategy.** The implementation enumerates multi-indices by
total degree \\d = 0, 1, 2, \ldots\\, generating all weak compositions
of \\d\\ into \\N\\ parts via recursive backtracking. For each degree,
the degree contribution is accumulated, and convergence is declared when
the absolute change between consecutive degrees falls below `tol`.
Specialized implementations for \\n = 2\\ (single-variable series in
\\\alpha\_{12}\\) and \\n = 3\\ (three-variable series in
\\\alpha\_{12}, \alpha\_{13}, \alpha\_{23}\\) avoid the overhead of the
general partition enumeration.

**Complexity.** The number of weak compositions of \\d\\ into \\N\\
parts is \\\binom{d + N - 1}{N - 1}\\, which grows combinatorially with
both \\d\\ and \\N\\. For \\n \geq 5\\ this becomes expensive, though
convergence typically occurs at low \\d\\ for well-conditioned cones.

**Implementation:** `R/series_hypergeometric.R` exports
`hypergeometric_series(V, max_terms, tol, check_pd)` as the dispatcher,
`hypergeometric_series_nd(V, max_terms, tol)` as the general
\\n\\-dimensional engine, and internal helpers
`hypergeometric_series_2d()`, `hypergeometric_series_3d()`,
`compute_series_term_3d()`.

  

### Tridiagonal series (Theorem 4.1)

When the Gram matrix \\V^T V\\ is tridiagonal (i.e., \\\hat{v}\_i \cdot
\hat{v}\_j = 0\\ for \\\|i - j\| \> 1\\), the hypergeometric series
simplifies dramatically. The summation reduces from \\N = \binom{n}{2}\\
variables to \\n - 1\\ consecutive dot products \\\beta_i = \hat{v}\_i
\cdot \hat{v}\_{i+1}\\:

\\ \omega(C) = \frac{\|\det V\|}{(4\pi)^{n/2}} \sum\_{\mathbf{b} \in
\mathbb{N}\_0^{n-1}} \left\[ \frac{(-2)^{\|\mathbf{b}\|}}{\mathbf{b}!}\\
\Gamma\\\left(\frac{1+b_1}{2}\right) \prod\_{i=2}^{n-1}
\Gamma\\\left(\frac{1+b\_{i-1}+b_i}{2}\right)
\Gamma\\\left(\frac{1+b\_{n-1}}{2}\right) \right\]
\boldsymbol{\beta}^{\mathbf{b}} \\

By Theorem 4.1 of [\[1\]](#ref1), tridiagonality of \\V^T V\\ implies
that \\M_n(C)\\ is automatically positive definite, so this series
always converges. The reduced dimensionality (\\n - 1\\ instead of
\\\binom{n}{2}\\) yields a substantial computational advantage,
particularly in high dimensions.

**Implementation:** `R/series_tridiagonal.R` exports
`tridiagonal_series(V, max_terms, tol)` and internal
`compute_tridiagonal_term(b, beta)`. The multi-index enumeration uses
`generate_partitions(n, k)` from `R/utils.R`.

  

### Decomposition method (Theorem 3.3)

When \\M_n(C)\\ is not positive definite, the hypergeometric series does
not converge. The decomposition method, based on Brion-Vergne theory
[\[1\]](#ref1), recursively partitions the cone into sub-cones, each of
which either has positive definite associated matrix or has affine
dimension less than \\n\\ (contributing zero to the solid angle).

**Algorithm.** Given generators \\V = \[v_1, \ldots, v_n\]\\:

1.  If \\M_n(C)\\ is already positive definite, return \\V\\ as a single
    terminal cone.
2.  Compute \\\delta_i = \text{sign}(\hat{v}\_i \cdot \hat{v}\_n)\\ for
    \\i = 1, \ldots, n-1\\.
3.  If all \\\delta_i = 0\\ (\\\hat{v}\_n\\ orthogonal to all other
    generators), project onto the \\(n-1)\\-dimensional subspace
    orthogonal to \\\hat{v}\_n\\ and recurse.
4.  Otherwise, for each \\i\\ with \\\delta_i \neq 0\\, construct a
    decomposed cone \\C_i\\ with generators:

\\ u\_{i,k} = \begin{cases} v_k & \text{if } \delta_k = 0 \text{ or } k
= i \\ \epsilon\_{ik}\left(v_k - \frac{v_k \cdot v_n}{v_i \cdot v_n}
v_i\right) & \text{otherwise} \end{cases} \\

where \\\epsilon\_{ik}\\ is a sign factor from Equation 15 of
[\[1\]](#ref1). Each \\C_i\\ carries a sign \\s_i = \pm 1\\ computed
from the delta sequence. The solid angle of the original cone is
\\\omega(C) = \sum_i s_i \cdot \omega(C_i)\\.

5.  Recurse on any sub-cone whose associated matrix is still not
    positive definite.

By Corollary 3.4, the decomposition produces at most \\(n-1)!\\ terminal
cones. The worst-case complexity is therefore \\O(2^n)\\, though
well-conditioned cones typically decompose into far fewer pieces.

**Implementation:** `R/decomposition.R` exports
`decompose_cone(V, method)`,
`solid_angle_decomposition(V, max_terms, tol)`, and internal helpers
`decompose_by_line(V)`, `compute_epsilon(i, k, delta)`.

  

### Multivariate normal integration (Theorem 2.2)

Ribando’s approach [\[3\]](#ref3) reformulates the solid angle as the
probability that a multivariate normal random variable falls in the
positive orthant under a transformed coordinate system:

\\ \omega(C) = \int\_{\mathbb{R}\_{\geq 0}^n} \phi_n(\mathbf{x};\\
\mathbf{0},\\ \Sigma)\\ d\mathbf{x} \\

where \\\Sigma = (V^T V)^{-1}\\ and \\\phi_n\\ is the \\n\\-dimensional
Gaussian density. The algorithm computes \\B = V^T V\\, inverts to
obtain \\\Sigma\\, then evaluates the integral via
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html),
which uses the Genz-Bretz Quasi-Monte Carlo algorithm [\[4\]](#ref4)
with Halton sequences and importance sampling.

This method works for any cone where \\B\\ is positive definite,
provides results exact up to numerical integration tolerance, and avoids
the combinatorial explosion of the hypergeometric series. However, the
matrix inversion becomes ill-conditioned for large \\n\\ or nearly
dependent generators; the function warns when \\\kappa(B) \> 10^{10}\\
and is unreliable for \\n \> 20\\.

**Complexity:** \\O(n^3)\\ for matrix inversion plus \\O(M \times n^2)\\
for the Quasi-Monte Carlo integration, where \\M\\ is the number of
adaptive quadrature points.

**Implementation:** the dispatcher branch
`compute_solid_angle(V, method = "mvn")` (in `R/main_interface.R`)
computes \\B = V^\top V\\, recovers \\\Sigma = B^{-1}\\, and forwards
the integral to
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html).
The function returns the scalar normalized solid angle; the
per-dimension geometric mean is recovered as `omega^(1 / ncol(V))` when
needed.

  

### Geometric methods for special structures

The package provides specialized formulas for geometrically defined
cones in \\\mathbb{R}^3\\ that bypass the general-purpose machinery.

**Right circular cones.** The Mazonka formula [\[5\]](#ref5) (Equation
9) gives the normalized solid angle of a cone with apex half-angle
\\\theta\\:

\\ \omega = \frac{2\pi(1 - \cos\theta)}{4\pi} = \frac{1 - \cos\theta}{2}
\\

**Implementation:** `solid_angle_cone(theta)` in
`R/geometric_methods.R`.

**Polyhedral cones (spherical polygons).** The Van Oosterom-Strackee
formula [\[2\]](#ref2) computes the solid angle of a spherical triangle
defined by three unit vectors \\v_1, v_2, v_3\\:

\\ \tan\\\left(\frac{\Omega}{2}\right) = \frac{\|v_1 \cdot (v_2 \times
v_3)\|}{1 + v_1 \cdot v_2 + v_1 \cdot v_3 + v_2 \cdot v_3} \\

For spherical polygons with \\n \> 3\\ vertices, the function decomposes
into triangular fans from the first vertex and sums the contributions.

**Implementation:** `solid_angle_polyhedral(vertices)` in
`R/geometric_methods.R`.

**Cone segments.** A cone segment is the portion of a circular cone cut
by a plane at angle \\\gamma\\ to the axis. The Mazonka formula
[\[5\]](#ref5) (Equations 35, 39, 41) computes auxiliary angles:

\\ \phi = \arccos\\\left(\frac{\cos\gamma}{\tan\theta}\right), \qquad
\beta = \arccos\\\left(\frac{\sin\gamma}{\sin\theta}\right) \\

and the solid angle \\\Omega = 2(\phi - \cos\theta \cdot \beta)\\,
normalized by \\4\pi\\. When the analytic expression yields a negative
or invalid value, the function falls back to numerical integration of
the spherical-cap intersection via
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) with 2000
subdivisions.

**Implementation:** `solid_angle_cone_segment(theta, gamma)` in
`R/geometric_methods.R`.

**Intersecting cones.** For two circular cones with half-angles
\\\theta_1\\, \\\theta_2\\ and axis separation \\\alpha\\, the
intersection solid angle is computed by numerical integration of the
azimuthally admissible arc length:

\\ \Omega = \int_0^{\theta_1} L(\theta)\\\sin\theta\\d\theta \\

where \\L(\theta) = 2\arccos C(\theta)\\ with \\C(\theta) =
(\cos\theta_2 - \cos\alpha\cos\theta)/(\sin\alpha\sin\theta)\\, clamped
to \\\[-1, 1\]\\. Special cases (no overlap, containment, co-directed
axes, opposite axes, hemisphere intersections) are handled analytically
before falling back to integration.

**Implementation:**
`solid_angle_intersecting_cones(theta1, theta2, alpha, method, fallback, mc_samples, mc_seed)`
in `R/geometric_methods.R`.

**L’Huilier’s theorem.** The spherical excess of a triangle with
arc-length sides \\a, b, c\\ is:

\\ E =
4\arctan\sqrt{\tan\\\left(\frac{s}{2}\right)\tan\\\left(\frac{s-a}{2}\right)\tan\\\left(\frac{s-b}{2}\right)\tan\\\left(\frac{s-c}{2}\right)}
\\

where \\s = (a + b + c)/2\\. This provides a numerically stable
alternative to Girard’s theorem \\E = \alpha + \beta + \gamma - \pi\\.

**Implementation:** `lhuilier_angle(a, b, c)` in `R/extreme_rays.R`.

  

## Uniform cone sampling

The package implements the \\O(n)\\ algorithm of Arun & Venkatapathi
[\[6\]](#ref6) for generating points uniformly distributed on a
spherical cap (the intersection of a cone with \\S^{n-1}\\) in arbitrary
dimension.

  

### Angle-to-solid-angle mapping

The function \\\Theta: \[0, \pi\] \to \[0, 1\]\\ maps a planar
half-angle \\\theta\\ to the fraction of the sphere enclosed by the
corresponding cone:

\\ \Theta(\theta) = \begin{cases} \frac{1}{2}\\I\\\left(\sin^2\theta;\\
\frac{n-1}{2},\\ \frac{1}{2}\right) & \theta \in \left\[0,
\frac{\pi}{2}\right\] \\ 1 - \frac{1}{2}\\I\\\left(\sin^2\theta;\\
\frac{n-1}{2},\\ \frac{1}{2}\right) & \theta \in \left(\frac{\pi}{2},
\pi\right\] \end{cases} \\

where \\I(x; \alpha, \beta)\\ is the regularized incomplete beta
function, computed via [`pbeta()`](https://rdrr.io/r/stats/Beta.html).
The inverse \\\Theta^{-1}\\ uses
[`qbeta()`](https://rdrr.io/r/stats/Beta.html). The piecewise definition
avoids catastrophic cancellation near \\\theta = \pi/2\\.

**Implementation:** `theta_to_omega(theta, n)` and
`omega_to_theta(omega, n)` in `R/cone_sampling.R`.

  

### The O(n) sampling algorithm

Given a cone axis \\\hat{\mu} \in S^{n-1}\\ and half-angle \\\theta_0\\,
a uniformly distributed point on the spherical cap is generated in four
steps:

**Step 1. Generate random planar angle.** Sample \\\theta \in \[0,
\theta_0\]\\ from the distribution \\F\_\theta(\theta) = \Theta(\theta)
/ \Theta(\theta_0)\\. Two methods are available:

- *Inverse transform sampling:* generate \\U \sim \text{Uniform}(0,
  \Theta(\theta_0))\\ and return \\\theta = \Theta^{-1}(U)\\. Fast but
  vulnerable to floating-point underflow for large \\n\\ and small
  \\\theta_0\\.

- *1D rejection sampling:* generate \\U \sim \text{Uniform}(0, 1)\\ and
  \\\theta \sim \text{Uniform}(0, \theta_0)\\; accept if \\\log U \<
  (n-2)\[\log\sin\theta - \log\sin\min(\theta_0, \pi/2)\]\\. The log
  transform prevents underflow. Expected number of rejections is
  \\O(n)\\, maintaining overall \\O(n)\\ complexity.

**Step 2. Generate perpendicular component.** Sample \\\hat{x}\_\perp\\
uniformly on \\S^{n-2}\\ using the Marsaglia method [\[7\]](#ref7): draw
\\z_1, \ldots, z\_{n-1} \sim N(0, 1)\\ and normalize. The do-while loop
rejects samples with \\\\\mathbf{z}\\ \< 10^{-15}\\.

**Step 3. Construct canonical vector.** Form \\\hat{x} = (\sin\theta
\cdot \hat{x}\_\perp,\\ \cos\theta)\\ aligned with the \\n\\-th
canonical basis vector \\\hat{e}\_n\\.

**Step 4. Rotate to target orientation.** Apply a Givens rotation in the
plane containing \\\hat{e}\_n\\ and \\\hat{\mu}\\:

\\ \hat{y} = \hat{x} + P(G - I_2)P^T\hat{x} \\

where \\P = \[\hat{e}\_n,\\ (\hat{\mu} - \mu_n\hat{e}\_n)/\\\hat{\mu} -
\mu_n\hat{e}\_n\\\]\\ is the \\n \times 2\\ orthonormal basis for the
rotation plane and \\G\\ is the \\2 \times 2\\ Givens matrix:

\\ G = \begin{bmatrix} \mu_n & -\sqrt{1 - \mu_n^2} \\ \sqrt{1 - \mu_n^2}
& \mu_n \end{bmatrix} \\

This rotation costs \\O(n)\\ operations (two dot products and two vector
additions), unlike a general \\n\\-dimensional rotation matrix which
would cost \\O(n^2)\\.

**Edge cases.** When \\\hat{\mu}\\ is nearly parallel to \\\hat{e}\_n\\
(\\\|1 - \mu_n^2\| \< 10^{-10}\\), the second column of \\P\\ is set to
zero, and the rotation reduces to the identity (parallel) or a
reflection (antiparallel). The implementation guards against negative
arguments to `sqrt` via `max(0, 1 - \mu_n^2)`.

**Hollow cones.** The algorithm extends to annular regions \\\theta_1
\leq \theta \leq \theta_2\\ by sampling \\U \sim
\text{Uniform}(\Theta(\theta_1), \Theta(\theta_2))\\ and applying
\\\theta = \Theta^{-1}(U)\\.

**Implementation (R):** `R/cone_sampling.R` exports
[`generate_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md),
[`generate_hollow_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md),
[`generate_planar_angle_inverse()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md),
[`generate_planar_angle_rejection()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md),
[`rotate_from_canonical()`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md),
[`generate_point_on_sphere()`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md),
[`rejection_cost()`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md),
[`verify_cone_uniformity()`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md).

  

## C++ computational backend

  

### Architecture

The C++ layer consists of a single source file `src/cone_sampling.cpp`
(124 lines) implementing the vectorized cone sampling algorithm via
RcppArmadillo. The design mirrors the R implementation exactly but
operates in a tight loop without R-level overhead, providing significant
speedup for large sample counts.

The package links against RcppArmadillo
(`// [[Rcpp::depends(RcppArmadillo)]]`) for linear algebra operations.
No OpenMP parallelization is used; the performance gain comes entirely
from eliminating the R interpreter overhead in the sampling loop.

  

### C++ helper functions (not exported)

**`theta_to_omega_cpp(double theta, int n)`** and
**`omega_to_theta_cpp(double omega, int n)`** mirror the R functions
exactly, calling `R::pbeta()` and `R::qbeta()` from the R math library
for the regularized incomplete beta function.

**`generate_point_on_sphere_cpp(int n)`** generates a uniform random
point on \\S^{n-1}\\ using `arma::randn(n)` for Gaussian sampling and
`arma::norm(x, 2)` for normalization. The do-while loop guards against
near-zero norms.

**`rotate_from_canonical_cpp(arma::vec x, arma::vec mu_hat)`** performs
the \\O(n)\\ Givens rotation using Armadillo matrix operations.
Constructs the \\n \times 2\\ projection matrix \\P\\ and the \\2 \times
2\\ Givens matrix \\G\\, then computes \\\hat{y} = \hat{x} + P(G -
I_2)P^T\hat{x}\\ via `P.t() * x` and matrix-vector products. Guards
against \\\|1 - \mu_n^2\| \leq 10^{-10}\\ (axis-aligned case) and uses
`std::max(0.0, ...)` for the square root argument.

  

### Exported function

``` cpp
// [[Rcpp::export]]
NumericMatrix generate_cone_samples_cpp(int n_samples, NumericVector mu_hat_r, double theta0)
```

This is the only function exported to R. It converts the R
`NumericVector` to an `arma::vec`, precomputes \\\Omega_0 =
\Theta(\theta_0)\\, then loops over `n_samples` iterations. Each
iteration performs inverse transform sampling (\\U \sim
\text{Uniform}(0, \Omega_0)\\, \\\theta = \Theta^{-1}(U)\\), generates a
perpendicular component on \\S^{n-2}\\, constructs the canonical vector,
applies the Givens rotation, and stores the result column-wise in an \\n
\times n\_{\text{samples}}\\ `NumericMatrix`.

The R wrapper
[`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
attempts the C++ backend first and falls back to a pure R loop if the
compiled code is unavailable.

C++ sampling pipeline. Each blue node is a per-sample stage; the green
nodes are wrappers, and the red node is the precomputed cap-area
constant. The pipeline runs N times in a tight loop inside
generate_cone_samples_cpp.

  

### Compilation

The package uses the standard Rcpp compilation pipeline. No `Makevars`
file is needed beyond the defaults provided by
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html),
since no OpenMP or custom compiler flags are required. The `NAMESPACE`
directive `useDynLib(SolidAngleR, .registration = TRUE)` registers the
compiled routines at load time.

  

## Software architecture

The package follows a two-layer design: an R orchestration layer handles
method selection, input validation, diagnostics, and user interaction,
while the C++ backend accelerates the single computationally intensive
operation (vectorized cone sampling).

  

### R layer

    R/
      SolidAngleR-package.R        Package-level documentation
      main_interface.R             compute_solid_angle (incl. method = "mvn"),
                                   compute_solid_angles, diagnose_cone, S3 triad
      utils.R                      Console output, dependency checking, matrix utilities
      formulas_2d.R                solid_angle_2d, solid_angle_2d_inner
      formulas_3d.R                solid_angle_3d, solid_angle_3d_det, spherical_triangle_area
      series_hypergeometric.R      hypergeometric_series, hypergeometric_series_nd
      series_tridiagonal.R         tridiagonal_series
      decomposition.R              decompose_cone, solid_angle_decomposition
      geometric_methods.R          solid_angle_cone, solid_angle_polyhedral,
                                   solid_angle_cone_segment, solid_angle_intersecting_cones,
                                   solid_angle_monte_carlo
      extreme_rays.R               angle_between, lhuilier_angle, spanning trees
      cone_sampling.R              generate_cone_sample(s), rotate_from_canonical,
                                   theta_to_omega, omega_to_theta, verification
      plot_cone_3d.R               3D visualization (plotly)
      plotIntersectingCones.R      Interactive cone intersection visualization
      shinySolidAngleR.R           Shiny dashboard
      RcppExports.R                Auto-generated C++ interface

The standalone `omega.R` source from versions \<= 0.5.1 was removed in
0.6.0; the multivariate-normal orthant integral is now exposed as
`compute_solid_angle(V, method = "mvn")`, dispatched through an internal
helper in `main_interface.R` that calls
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html).

  

### C++ layer

    src/
      cone_sampling.cpp            O(n) sampling algorithm (RcppArmadillo)
      RcppExports.cpp              Auto-generated Rcpp registration

  

### Method selection logic

The main entry point `compute_solid_angle(V, method = "auto")` first
normalises the columns of \\V\\ and tests linear independence; if either
step fails the function returns `NA` with a warning. With
`method = "auto"`, the dispatcher consults a `cone_diagnosis` object and
routes the call to the cheapest convergent backend. The full decision
tree is shown below.

Dispatcher decision tree of compute_solid_angle(V, method = ‘auto’).
Diamonds are predicates; rectangles are computational backends. The
‘mvn’ branch is selected explicitly by the user; the ‘auto’ path never
returns it.

For forced method selection (`method = "formula"`, `"series"`,
`"tridiagonal"`, `"decomposition"`, `"mvn"`), validation checks are
applied (the formula method errors for \\n \> 3\\; the series method
errors if \\M_n(C)\\ is not positive definite; the `mvn` branch warns
when \\\kappa(V^\top V) \> 10^{10}\\ and refuses \\n \> 20\\).

  

### S3 class: `cone_diagnosis`

The `diagnose_cone(V)` function returns an S3 object of class
`"cone_diagnosis"` with the following structure:

| Component | Type | Description |
|:---|:---|:---|
| `dimension` | `integer` | Ambient dimension \\n\\ |
| `linearly_independent` | `logical` | Whether generators are linearly independent |
| `associated_matrix_PD` | `logical` | Whether \\M_n(C)\\ is positive definite |
| `is_tridiagonal` | `logical` | Whether \\V^T V\\ is tridiagonal |
| `eigenvalues` | `numeric(n)` | Eigenvalues of \\M_n(C)\\ |
| `min_eigenvalue` | `numeric` | Smallest eigenvalue |
| `suggested_method` | `character` | Recommended computation method |

Three S3 methods provide the standard inspection interface. The
[`print()`](https://rdrr.io/r/base/print.html) method gives a one-screen
diagnostic report with Unicode indicators;
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`"summary.cone_diagnosis"` object with condition number \\\kappa =
\max\|\lambda_i\| / \min\|\lambda_i\|\\, eigenvalue spread, stability
classification (excellent / good / moderate / poor / critical based on
\\\kappa\\ thresholds at 10, 100, \\10^6\\, \\10^{10}\\), and a prose
rationale for the method recommendation;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) renders an
eigenvalue ggplot panel coloured by sign (blue for positive, red for
negative), with a dashed reference line at zero.

S3 lifecycle of the cone_diagnosis class. The constructor
diagnose_cone() returns a list-based object that flows through the print
/ summary / plot trio and feeds the dispatcher’s auto-selection logic.

  

## API reference

  

### Core interface

``` r

compute_solid_angle(V, method = "auto", max_terms = 1000, tol = 1e-10, normalize = TRUE)
```

The primary entry point. `V` is an \\n \times n\\ matrix whose columns
generate the simplicial cone. Returns a numeric value in \\\[0, 1\]\\.

``` r

compute_solid_angles(cone_list, ...)
```

Vectorized computation over a list of cone matrices via `sapply`.

``` r

diagnose_cone(V)
```

Returns a `cone_diagnosis` S3 object with structural analysis and method
recommendation.

  

### Series methods

| Function | Signature | Convergence condition |
|:---|:---|:---|
| `hypergeometric_series` | `(V, max_terms, tol, check_pd)` | \\M_n(C)\\ positive definite |
| `hypergeometric_series_nd` | `(V, max_terms, tol)` | \\M_n(C)\\ positive definite |
| `tridiagonal_series` | `(V, max_terms, tol)` | \\V^T V\\ tridiagonal (auto-PD) |

All series methods return `list(solid_angle, n_terms, converged, ...)`.

  

### Decomposition

| Function | Signature | Output |
|:---|:---|:---|
| `decompose_cone` | `(V, method = "line")` | `list(cones, signs, n_cones, all_pd)` |
| `solid_angle_decomposition` | `(V, max_terms, tol)` | Numeric solid angle |

  

### Geometric methods

| Function | Signature | Domain |
|:---|:---|:---|
| `solid_angle_cone` | `(theta)` | Right circular cone, \\\theta \in (0, \pi)\\ |
| `solid_angle_polyhedral` | `(vertices)` | Spherical polygon, \\n \times 3\\ matrix |
| `solid_angle_cone_segment` | `(theta, gamma)` | Cone segment cut by plane |
| `solid_angle_intersecting_cones` | `(theta1, theta2, alpha, method, fallback, mc_samples, mc_seed)` | Two-cone intersection |
| `solid_angle_monte_carlo` | `(region_test, n_samples, seed)` | Arbitrary region via MC |
| `spherical_triangle_area` | `(v1, v2, v3)` | Spherical excess in steradians |
| `lhuilier_angle` | `(a, b, c)` | L’Huilier’s theorem |

  

### Multivariate normal integration

``` r

compute_solid_angle(V, method = "mvn")
```

Returns the scalar normalized solid angle in \\\[0, 1\]\\. The
per-dimension geometric mean \\\omega^{1/n}\\ is recovered inline as
`result^(1 / ncol(V))`. Requires `mvtnorm`.

  

### Cone sampling

| Function | Signature | Complexity |
|:---|:---|:---|
| `generate_point_on_sphere` | `(n)` | \\O(n)\\ |
| `generate_cone_sample` | `(mu_hat, theta0, method)` | \\O(n)\\ |
| `generate_cone_samples` | `(n_samples, mu_hat, theta0, method)` | \\O(n \cdot n\_{\text{samples}})\\ (C++) |
| `generate_hollow_cone_sample` | `(mu_hat, theta1, theta2, method)` | \\O(n)\\ |
| `theta_to_omega` | `(theta, n)` | \\O(1)\\ |
| `omega_to_theta` | `(omega, n)` | \\O(1)\\ |
| `generate_planar_angle_inverse` | `(theta0, n)` | \\O(1)\\ |
| `generate_planar_angle_rejection` | `(theta0, n)` | \\O(n)\\ expected |
| `rotate_from_canonical` | `(x, mu_hat)` | \\O(n)\\ |
| `verify_cone_uniformity` | `(samples, mu_hat, theta0, n_bins)` | \\O(n\_{\text{samples}})\\ |
| `rejection_cost` | `(theta0, n)` | \\O(1)\\ |

  

### Matrix and cone utilities

| Function | Signature | Purpose |
|:---|:---|:---|
| `compute_associated_matrix` | `(V)` | \\M_n(C)\\ from Definition 1.1 |
| `compute_dot_product_matrix` | `(V)` | \\V^T V\\ |
| `is_positive_definite` | `(M, tol)` | Eigenvalue-based PD check |
| `is_linearly_independent` | `(V, tol)` | QR rank test |
| `is_tridiagonal` | `(M, tol)` | Off-band element test |
| `normalize_vectors` | `(V)` | Column-wise unit normalization |
| `cross_product_3d` | `(a, b)` | 3D cross product |
| `create_tridiagonal_cone` | `(angles)` | Cholesky-based tridiagonal cone constructor |

  

### Visualization and interaction

``` r

plotIntersectingCones(theta1, theta2, alpha, ...)   # Interactive 3D (plotly)
shinySolidAngleR(launch.browser = TRUE)             # Shiny dashboard
```

  

## Numerical considerations

**Dot product clamping.** All functions that compute `acos(x)` clamp the
argument to \\\[-1, 1\]\\ to prevent `NaN` from floating-point
overshoot. This arises in
[`solid_angle_3d()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md),
[`angle_between()`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md),
[`lhuilier_angle()`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md),
and the cone sampling verification.

**Degenerate cone detection.** A cone is considered degenerate (zero
solid angle) when the scalar triple product is below \\10^{-15}\\ (3D
formulas), the rank of \\V\\ is less than \\n\\ (decomposition), or the
norm of a generated random vector is below \\10^{-15}\\ (sphere
sampling).

**Condition number monitoring.** The
`compute_solid_angle(method = "mvn")` branch computes \\\kappa(B)\\
exactly via `kappa(B, exact = TRUE)` and warns when \\\kappa \>
10^{10}\\. The `summary.cone_diagnosis` method classifies stability
using thresholds at \\\kappa \in \\10, 100, 10^6, 10^{10}\\\\.

**Convergence tolerance.** All series methods default to `tol = 1e-10`
for convergence detection. Convergence is assessed per-degree: the
series terminates when the absolute contribution of an entire degree
level (sum of all terms at that total degree) falls below `tol`, with a
minimum of 3 degrees to avoid premature termination.

**Numerical integration fallback.** The
[`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
function uses
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) with
`rel.tol = 1e-10`, `abs.tol = 1e-10`, and `subdivisions = 2000L` when
the analytic Mazonka formula yields negative or non-finite results.

  

## Performance characteristics

| Method | Complexity | Typical use case |
|:---|:---|:---|
| Closed formulas (2D/3D) | \\O(1)\\ | Low-dimensional cones |
| Hypergeometric series | \\O(M \cdot \binom{d+N-1}{N-1})\\ | PD cones, \\n \leq 6\\ |
| Tridiagonal series | \\O(M \cdot \binom{d+n-2}{n-2})\\ | Tridiagonal \\V^T V\\, any \\n\\ |
| Decomposition | \\O((n-1)! \cdot T\_{\text{series}})\\ worst case | Non-PD cones |
| MVN integration | \\O(n^3 + M \cdot n^2)\\ | Any PD cone, \\n \leq 20\\ |
| Cone sampling (C++) | \\O(n \cdot n\_{\text{samples}})\\ | Uniform sampling, any \\n\\ |
| Monte Carlo | \\O(n\_{\text{samples}})\\ | Verification, complex regions |

  

## References

**\[1\]** Fitisone, A., & Zhou, Y. (2023). Solid angle measure of
polyhedral cones. *arXiv:2304.11102* (math.CO).
<https://arxiv.org/abs/2304.11102>

**\[2\]** Van Oosterom, A., & Strackee, J. (1983). The solid angle of a
plane triangle. *IEEE Transactions on Biomedical Engineering*,
BME-30(2), 125-126. DOI:
[10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

**\[3\]** Ribando, J. M. (2006). Measuring solid angles beyond dimension
three. *Discrete & Computational Geometry*, 36(3), 479-487. DOI:
[10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

**\[4\]** Genz, A., & Bretz, F. (2009). *Computation of multivariate
normal and t probabilities*. Lecture Notes in Statistics, Vol. 195.
Springer-Verlag. DOI:
[10.1007/978-3-642-01689-9](https://doi.org/10.1007/978-3-642-01689-9)

**\[5\]** Mazonka, O. (2012). Solid angle of conical surfaces,
polyhedral cones, and intersecting spherical caps. *arXiv:1205.1396v2*.
DOI:
[10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

**\[6\]** Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for
generating uniform random vectors in n-dimensional cones. *Sankhya A:
The Indian Journal of Statistics*, 87(2), 327-348. DOI:
[10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

**\[7\]** Marsaglia, G. (1972). Choosing a point from the surface of a
sphere. *Annals of Mathematical Statistics*, 43(2), 645-646. DOI:
[10.1214/aoms/1177692644](https://doi.org/10.1214/aoms/1177692644)

**\[8\]** Aomoto, K. (1977). Analytic structure of Schlafli function.
*Nagoya Mathematical Journal*, 68, 1-16.
<https://projecteuclid.org/euclid.nmj/1118786429>

**\[9\]** Song, C., Rohr, R. P., & Saavedra, S. (2018). A guideline to
study the feasibility domain of multi-trophic and changing ecological
communities. *Journal of Theoretical Biology*, 450, 30-36. DOI:
[10.1016/j.jtbi.2018.04.030](https://doi.org/10.1016/j.jtbi.2018.04.030)
