# Solid angle via Ribando's multivariable hypergeometric series

Computes the normalized solid angle of a simplicial cone using the
multivariable hypergeometric series of Ribando (2006), with the
dimension- specific specialisations of Fitisone & Zhou (2023).
Internally dispatches to a hand-rolled implementation for \\n = 2\\, a
hand-rolled implementation for \\n = 3\\, and the general recursive
enumerator
[`hypergeometric_series_nd`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md)
for \\n \geq 4\\.

## Usage

``` r
hypergeometric_series(V, max_terms = 1000, tol = 1e-10, check_pd = TRUE)
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are the cone
  generators. Need not be unit-normalized; columns are rescaled
  internally.

- max_terms:

  Integer. Hard cap on the number of series terms accumulated before
  truncation. Default `1000`.

- tol:

  Numeric. Per-degree convergence tolerance. The series terminates when
  the contribution of an entire degree level falls below `tol`. Default
  `1e-10`.

- check_pd:

  Logical. When `TRUE` (default), the function verifies that the
  associated matrix \\M_n(C)\\ is positive definite before iterating,
  since this is the necessary and sufficient condition for absolute
  convergence of the series.

## Value

A list with components

- solid_angle:

  Normalized solid angle in \\\[0, 1\]\\.

- n_terms:

  Total number of series terms accumulated.

- converged:

  Logical, `TRUE` when the per-degree contribution fell below `tol`
  before `max_terms` was reached.

- associated_matrix:

  The associated matrix \\M_n(C)\\.

## Details

The Ribando-Aomoto series (equation 5 of Ribando 2006, equation 1 of
Fitisone & Zhou 2023) is \$\$\Omega(C) = \frac{\|\det V\|}{(4\pi)^{n/2}}
\sum\_{a \in \mathbb{N}\_0^N} \frac{(-2)^{\|a\|}}{a!}\\ \prod\_{i=1}^n
\Gamma\\\left(\frac{1 + \deg_i(a)}{2}\right)\\ \alpha^{a},\$\$ where \\N
= \binom{n}{2}\\, \\a = (a\_{ij})\_{i\<j}\\ is a multi-index over pairs
of generator indices, \\\|a\| = \sum\_{i\<j} a\_{ij}\\, \\\deg_i(a) =
\sum\_{k \neq i} a\_{ik}\\, and \\\alpha\_{ij} = v_i \cdot v_j\\. The
series converges absolutely to \\\Omega(C)\\ if and only if the
associated matrix \\M_n(C) = (m\_{ij})\\ with \\m\_{ii} = 1\\ and
\\m\_{ij} = -\|\alpha\_{ij}\|\\ for \\i \neq j\\ is positive definite
(Ribando 2006, Theorem 1.5).

The implementation enumerates terms by total degree \\\|a\|\\ so the
truncation error after the cut-off can be bounded by the contribution of
the next degree level. For \\n = 2, 3\\ the dispatcher delegates to
specialised implementations that exploit the small number of pair
indices.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Aomoto, K. (1977). Analytic structure of Schlafli function. *Nagoya
Mathematical Journal*, 68, 1-16.
<https://projecteuclid.org/euclid.nmj/1118786429>

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`hypergeometric_series_nd`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_nd.md)
for the general recursive enumerator used when \\n \geq 4\\;
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the simplified \\(n-1)\\-variable series available when \\V^\top V\\
is tridiagonal;
[`solid_angle_decomposition`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
for the decomposition route applicable when \\M_n(C)\\ is not positive
definite;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
# 3D orthogonal cone: 1/8 of the sphere
V <- diag(3)
hypergeometric_series(V, max_terms = 100)$solid_angle    # 0.125

# 4D orthogonal cone: 1/16 of the sphere
V <- diag(4)
hypergeometric_series(V, max_terms = 200)$solid_angle    # 0.0625
} # }
```
