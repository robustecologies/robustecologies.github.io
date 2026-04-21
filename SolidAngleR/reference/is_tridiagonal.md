# Check if a matrix is tridiagonal

Checks whether a matrix has only diagonal and first off-diagonal entries
non-zero (up to tolerance).

## Usage

``` r
is_tridiagonal(M, tol = 1e-12)
```

## Arguments

- M:

  A matrix

- tol:

  Tolerance for zero comparison (default: 1e-12)

## Value

Logical: TRUE if tridiagonal, FALSE otherwise

## Examples

``` r
M1 <- matrix(c(1, 0.5, 0, 0.5, 1, 0.3, 0, 0.3, 1), nrow = 3)
is_tridiagonal(M1)  # TRUE
#> [1] TRUE

M2 <- matrix(c(1, 0.5, 0.1, 0.5, 1, 0.3, 0.1, 0.3, 1), nrow = 3)
is_tridiagonal(M2)  # FALSE (M[1,3] != 0)
#> [1] FALSE
```
