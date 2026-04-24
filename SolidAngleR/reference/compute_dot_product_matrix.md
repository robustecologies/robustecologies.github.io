# Pairwise dot-product matrix V^T V

Returns the Gram matrix \\V^\top V\\ of pairwise inner products of the
columns of `V`. Used as the discriminant for the tridiagonal- structure
test of
[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md).

## Usage

``` r
compute_dot_product_matrix(V)
```

## Arguments

- V:

  A numeric matrix whose columns are vectors in \\\mathbb{R}^n\\.

## Value

A symmetric numeric matrix of dimension `ncol(V) %*% ncol(V)` with entry
\\(i, j)\\ equal to \\v_i \cdot v_j\\.

## Details

The Gram matrix is positive semidefinite for any `V` and positive
definite when the columns are linearly independent. Tridiagonality of
\\V^\top V\\ is the structural condition that activates the simplified
series of
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
(Fitisone & Zhou 2023, Theorem 4.1).

## References

Strang, G. (2016). *Introduction to Linear Algebra*, 5th edition.
Wellesley-Cambridge Press. ISBN 978-0980232776.

## See also

[`is_tridiagonal`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
for the tridiagonality test on the result;
[`compute_associated_matrix`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md)
for the closely related Ribando matrix \\M_n(C)\\;
[`tridiagonal_series`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
for the tridiagonal-aware backend.

## Examples

``` r
V <- cbind(c(1, 0), c(0, 1), c(1, 1) / sqrt(2))
compute_dot_product_matrix(V)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.0000000 0.7071068
#> [2,] 0.0000000 1.0000000 0.7071068
#> [3,] 0.7071068 0.7071068 1.0000000
```
