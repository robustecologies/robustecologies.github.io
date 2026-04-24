# Test whether a symmetric matrix is positive definite

Returns `TRUE` when every eigenvalue of the symmetric matrix `M` exceeds
the tolerance `tol`. Used internally by the dispatcher to decide whether
the Ribando hypergeometric series can be applied directly or whether the
decomposition route is required.

## Usage

``` r
is_positive_definite(M, tol = 1e-10)
```

## Arguments

- M:

  A symmetric numeric matrix.

- tol:

  Numeric. Eigenvalue threshold below which a value is considered
  non-positive. Default `1e-10`.

## Value

A single logical: `TRUE` if every eigenvalue exceeds `tol`, `FALSE`
otherwise.

## Details

Eigenvalues are computed with
[`eigen`](https://rdrr.io/r/base/eigen.html) on the symmetric pathway,
which is numerically stable and faster than the general non-symmetric
routine. The threshold form (\\\lambda \> tol\\) rather than \\\lambda
\> 0\\ guards against eigenvalues that are numerically zero from
finite-precision arithmetic on near-singular matrices.

## References

Horn, R. A., & Johnson, C. R. (2013). *Matrix Analysis*, 2nd edition.
Cambridge University Press. ISBN 978-0521548236. (Sylvester criterion
and eigenvalue characterisation of positive definiteness.)

## See also

[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md)
for the matrix on which this test is most often invoked;
[`is_linearly_independent`](https://robustecologies.github.io/SolidAngleR/reference/is_linearly_independent.md)
for the rank-based degeneracy test;
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the tridiagonal-structure test;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher that uses these predicates.

## Examples

``` r
is_positive_definite(matrix(c(2, -0.5, -0.5, 2), nrow = 2))   # TRUE
#> [1] TRUE
is_positive_definite(matrix(c(1, -2,  -2,   1), nrow = 2))    # FALSE
#> [1] FALSE
```
