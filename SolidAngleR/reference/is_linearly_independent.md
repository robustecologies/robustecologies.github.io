# Test whether the columns of a matrix are linearly independent

Returns `TRUE` when the QR-decomposition rank of `V` equals the number
of columns, signalling that the columns of `V` are linearly independent
and so define a non-degenerate simplicial cone.

## Usage

``` r
is_linearly_independent(V, tol = 1e-10)
```

## Arguments

- V:

  A numeric matrix whose columns are vectors in \\\mathbb{R}^n\\.

- tol:

  Numeric. Reserved for future tolerance control; the current
  implementation delegates the rank test to
  [`qr`](https://rdrr.io/r/base/qr.html), whose default tolerance is
  used.

## Value

A single logical: `TRUE` if `ncol(V)` columns are linearly independent,
`FALSE` otherwise.

## Details

Linear independence of the generators is a precondition for every
solid-angle backend in the package. The dispatcher
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
issues a warning when this test fails because the spanned cone is then
lower-dimensional and has solid-angle measure zero.

## References

Strang, G. (2016). *Introduction to Linear Algebra*, 5th edition.
Wellesley-Cambridge Press. ISBN 978-0980232776.

## See also

[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md),
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for downstream consumers;
[`is_positive_definite`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md),
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for sibling tests.

## Examples

``` r
is_linearly_independent(cbind(c(1, 0), c(0, 1)))   # TRUE
#> [1] TRUE
is_linearly_independent(cbind(c(1, 0), c(2, 0)))   # FALSE
#> [1] FALSE
```
