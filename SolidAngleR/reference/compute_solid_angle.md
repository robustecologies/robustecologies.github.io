# Compute the normalized solid angle of a polyhedral cone

High-level dispatcher that returns the normalized solid angle measure of
a polyhedral cone in arbitrary ambient dimension. The function inspects
the dimension and the spectral structure of the associated matrix and
routes the computation to the most appropriate analytical, series,
decomposition, or quasi-Monte Carlo backend.

## Usage

``` r
compute_solid_angle(
  V,
  method = "auto",
  max_terms = 1000,
  tol = 1e-10,
  normalize = TRUE
)
```

## Arguments

- V:

  A numeric matrix whose columns are the cone generators. For a
  simplicial cone the matrix is square (\\n \times n\\) with linearly
  independent columns.

- method:

  Character scalar selecting the backend. Allowed values:

  - `"auto"`: automatic selection based on dimension and the spectrum of
    the associated matrix (default).

  - `"formula"`: closed-form analytical formula. Available for \\n = 2\\
    (planar arc length) and \\n = 3\\ (Van Oosterom- Strackee
    tetrahedron formula).

  - `"series"`: Ribando's hypergeometric series. Requires the associated
    matrix \\M_n(C)\\ to be positive definite.

  - `"tridiagonal"`: Fitisone-Zhou tridiagonal series. Requires the Gram
    matrix \\V^\top V\\ to be tridiagonal.

  - `"decomposition"`: Fitisone-Zhou signed decomposition into sub-cones
    with positive definite associated matrices.

  - `"mvn"`: integration of a zero-mean multivariate normal with
    covariance \\B^{-1}\\ over the positive orthant via
    [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
    (Genz-Bretz quasi-Monte Carlo).

- max_terms:

  Integer. Maximum number of terms used by the series and decomposition
  backends before truncation. Default `1000`.

- tol:

  Numeric. Convergence tolerance for the iterative backends. Default
  `1e-10`.

- normalize:

  Logical. If `TRUE` (default) the input vectors are rescaled to unit
  Euclidean norm before any geometric computation.

## Value

A single numeric value in \\\[0, 1\]\\ giving the fraction of the unit
sphere \\S^{n-1}\\ covered by the cone. The orthant in dimension \\n\\
returns \\1/2^n\\.

## Details

The dispatcher implements method selection from Fitisone & Zhou (2023)
for series and decomposition routes, the closed-form formulas of Van
Oosterom & Strackee (1983) for low dimensions, and the Genz-Bretz
quasi-Monte Carlo orthant integration of Genz & Bretz (2009) for the
multivariate-normal route.

For \\n = 2\\ the cone reduces to a planar wedge of opening angle
\\\theta\\ and the normalized solid angle is \\\Omega = \theta /
(2\pi)\\.

For \\n = 3\\ the dispatcher uses the Van Oosterom-Strackee formula
\$\$E = 2 \arctan\left(\frac{\|v_1 \cdot (v_2 \times v_3)\|}{1 + v_1
\cdot v_2 + v_2 \cdot v_3 + v_3 \cdot v_1}\right),\$\$ giving \\\Omega =
E / (4\pi)\\.

For \\n \geq 4\\ with positive definite associated matrix the dispatcher
selects Ribando's hypergeometric series (Theorem 1.5 in Ribando 2006);
when the Gram matrix is additionally tridiagonal, the simplified series
of Fitisone & Zhou (2023, Theorem 4.1) is preferred for its lower
per-term cost. When the associated matrix is not positive definite, the
cone is decomposed by Theorem 3.3 and Corollary 3.4 of Fitisone & Zhou
(2023) into signed sub-cones whose individual solid angles are computed
by the series route and aggregated.

The `"mvn"` branch computes \\B = V^\top V\\, recovers the covariance
\\\Sigma = B^{-1}\\, and integrates a centred Gaussian with that
covariance over the positive orthant: \$\$\Omega = \int\_{\[0,
\infty)^n} \phi_n(x;\\ 0,\\ \Sigma)\\ dx,\$\$ where \\\phi_n\\ is the
\\n\\-variate normal density. The integral is evaluated by the
Genz-Bretz lattice rule with random shifts implemented in
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html),
achieving an effective convergence rate close to \\O(N^{-1})\\ for
sufficiently smooth integrands.

Numerical safeguards: when the condition number \\\kappa(B) \> 10^{10}\\
the `"mvn"` branch issues a warning. For ambient dimension \\n \> 20\\ a
high-dimension warning is also emitted, since both the series cost and
the orthant integration accuracy degrade rapidly.

Symmetries of the measure: the normalized solid angle is invariant under
(i) any orthogonal transformation of the ambient space, \\\Omega(QV) =
\Omega(V)\\ for \\Q \in O(n)\\, which includes the antipodal symmetry
\\\Omega(-V) = \Omega(V)\\ as the special case \\Q = -I\\; (ii) positive
scaling of each column, \\\Omega(V\\\mathrm{diag}(d_1, \ldots, d_n)) =
\Omega(V)\\ for any \\d_i \> 0\\, since the cone is determined by its
rays, not by the column lengths; and (iii) permutation of columns, which
only relabels generators. These invariances are exact and hold for every
backend; in the `"mvn"` route they follow at once from the fact that \\V
\mapsto -V\\ and column scaling leave \\V^\top V\\ invariant (after
normalization), so the integrand and the integration region of the
Genz-Bretz orthant probability are unchanged.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane
triangle. *IEEE Transactions on Biomedical Engineering*, 30(2), 125-126.
[doi:10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207)

Genz, A., & Bretz, F. (2009). *Computation of multivariate normal and t
probabilities*. Lecture Notes in Statistics, Vol. 195. Springer-Verlag.
[doi:10.1007/978-3-642-01689-9](https://doi.org/10.1007/978-3-642-01689-9)

## See also

[`compute_solid_angles`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angles.md)
for vectorised computation over a list of cones;
[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
and
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
for the closed-form low-dimensional formulas;
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md),
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md),
[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
for the individual series and decomposition backends;
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the diagnostic S3 wrapper that recommends a method.

## Examples

``` r
# ====================================================================== #
# Example 1: planar wedge in R^2 (45 degree opening) ####
# ====================================================================== #

v1 <- c(1, 0)
v2 <- c(1, 1) / sqrt(2)
V <- cbind(v1, v2)
compute_solid_angle(V)            # 0.125
#> [1] 0.125

# ====================================================================== #
# Example 2: orthogonal octant in R^3 ####
# ====================================================================== #

V <- diag(3)
compute_solid_angle(V)            # 0.125 (one eighth of the sphere)
#> [1] 0.125

# ====================================================================== #
# Example 3: orthant in R^4 via the multivariate-normal branch ####
# ====================================================================== #

if (FALSE) { # \dontrun{
V <- diag(4)
compute_solid_angle(V, method = "mvn")  # 0.0625 = 1/2^4
} # }

# ====================================================================== #
# Example 4: tridiagonal cone via the simplified series ####
# ====================================================================== #

if (FALSE) { # \dontrun{
angles <- c(pi/3, pi/3, pi/3)
V <- create_tridiagonal_cone(angles)
compute_solid_angle(V, method = "tridiagonal")
} # }

# ====================================================================== #
# Example 5: general cone via the decomposition method ####
# ====================================================================== #

if (FALSE) { # \dontrun{
set.seed(123)
V <- matrix(rnorm(9), nrow = 3)
compute_solid_angle(V, method = "auto")
} # }
```
