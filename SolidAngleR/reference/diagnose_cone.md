# Diagnose a polyhedral cone and recommend a computation method

Inspects the geometric and spectral structure of a polyhedral cone and
returns a diagnostic S3 object summarising linear independence, positive
definiteness of the associated matrix, tridiagonal structure of the Gram
matrix, eigenvalue spectrum, and the recommended backend for
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md).

## Usage

``` r
diagnose_cone(V)
```

## Arguments

- V:

  A numeric matrix whose columns are the cone generators.

## Value

An object of class `cone_diagnosis`; a list with components

- dimension:

  Integer ambient dimension.

- linearly_independent:

  Logical, `TRUE` when the columns of `V` are linearly independent.

- associated_matrix_PD:

  Logical, positive definiteness of the associated matrix \\M_n(C)\\.

- is_tridiagonal:

  Logical, tridiagonal structure of \\V^\top V\\.

- eigenvalues:

  Numeric vector of eigenvalues of \\M_n(C)\\ in ascending order.

- min_eigenvalue:

  Smallest eigenvalue.

- suggested_method:

  Character string identifying the recommended backend for
  [`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md).

## Details

The diagnostic logic encodes the dispatcher's decision tree explicitly
so that the user can audit method selection without running the
computation. For \\n \leq 3\\ the closed-form formulas are always
recommended. For \\n \geq 4\\ the recommendation depends on the spectrum
of the associated matrix \\M_n(C) = (m\_{ij})\\ with \$\$m\_{ij} =
\langle v_i, v_j \rangle / (\\v_i\\\\ \\v_j\\),\$\$ positive
definiteness of which is required for the convergence of Ribando's
hypergeometric series (Ribando 2006, Theorem 1.5). When \\V^\top V\\ is
additionally tridiagonal the simplified series of Fitisone & Zhou (2023,
Theorem 4.1) reduces the per-term cost from \\O(n^2)\\ to \\O(n)\\.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## See also

[`print.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/print.cone_diagnosis.md),
[`summary.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/summary.cone_diagnosis.md),
[`plot.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
for the S3 method triad;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
as the dispatcher this diagnostic advises;
[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md),
[`is_positive_definite`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md),
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the underlying spectral tests.

## Examples

``` r
if (FALSE) { # \dontrun{
# Full lifecycle: constructor -> print -> summary -> plot
V <- diag(4)
d <- diagnose_cone(V)
print(d)
s <- summary(d)
print(s)
plot(d)
} # }
```
