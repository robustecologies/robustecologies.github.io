# Solid angle of a simplicial cone via signed decomposition

Computes the normalized solid angle of a simplicial cone by combining
the Fitisone-Zhou line-based decomposition with the hypergeometric or
tridiagonal series on each positive definite sub-cone, then aggregating
the signed contributions. Handles the general case where the associated
Gram matrix is not positive definite.

## Usage

``` r
solid_angle_decomposition(V, max_terms = 500, tol = 1e-10)
```

## Arguments

- V:

  A square \\n \times n\\ numeric matrix whose columns are linearly
  independent cone generators.

- max_terms:

  Integer. Maximum number of series terms used on each positive definite
  sub-cone. Default `500`.

- tol:

  Numeric. Convergence tolerance for the series backends. Default
  `1e-10`.

## Value

A single numeric value in \\\[0, 1\]\\: the normalized solid angle of
the input cone, equal to the sum of signed solid angles of the sub-cones
returned by
[`decompose_cone`](https://robustecologies.github.io/SolidAngleR/reference/decompose_cone.md).

## Details

The function applies the closed-form formula directly when \\n = 2\\ or
when \\n = 3\\ with positive definite associated matrix. Otherwise it
calls
[`decompose_cone`](https://robustecologies.github.io/SolidAngleR/reference/decompose_cone.md)
to obtain a signed family \\\\(C_i, s_i)\\\\ and accumulates
\$\$\Omega(C) = \sum_i s_i \\ \Omega(C_i),\$\$ where each
\\\Omega(C_i)\\ is computed by the closed-form formula in two or three
dimensions, by
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
when the Gram matrix of \\C_i\\ is tridiagonal, or by
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
otherwise. Sub-cones whose rank falls below \\n\\ contribute zero by
Theorem 3.3 case (I).

Computational complexity is at most \\O((n-1)!)\\ sub-cones by Corollary
3.4, but in practice the count is much smaller for well-conditioned
cones. Deep recursion can be triggered by adversarial configurations and
may exhaust the stack; the user-facing dispatcher
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
guards against this by selecting the decomposition route only when the
associated matrix is genuinely not positive definite.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`decompose_cone`](https://robustecologies.github.io/SolidAngleR/reference/decompose_cone.md)
for the underlying signed decomposition;
[`hypergeometric_series`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md),
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the per-sub-cone series backends;
[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md),
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
for the closed-form formulas;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
# 2D cone: closed-form short-circuit
V2 <- cbind(c(1, 0), c(0.5, sqrt(3) / 2))
solid_angle_decomposition(V2)
#> [1] 0.1666667

# 3D orthogonal cone: 1/8 of the sphere
V3 <- diag(3)
solid_angle_decomposition(V3)
#> [1] 0.125

if (FALSE) { # \dontrun{
# 3D cone requiring decomposition (non-PD associated matrix)
set.seed(123)
V3_npd <- normalize_vectors(matrix(rnorm(9), nrow = 3))
solid_angle_decomposition(V3_npd)
} # }
```
