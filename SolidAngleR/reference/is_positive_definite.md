# Check if a matrix is positive definite

Determines whether a symmetric matrix is positive definite by checking
if all eigenvalues are positive.

## Usage

``` r
is_positive_definite(M, tol = 1e-10)
```

## Arguments

- M:

  A symmetric matrix

- tol:

  Tolerance for numerical comparison (default: 1e-10)

## Value

Logical value: TRUE if positive definite, FALSE otherwise

## Examples

``` r
# Positive definite matrix
M1 <- matrix(c(2, -0.5, -0.5, 2), nrow = 2)
is_positive_definite(M1)  # TRUE
#> [1] TRUE

# Not positive definite
M2 <- matrix(c(1, -2, -2, 1), nrow = 2)
is_positive_definite(M2)  # FALSE
#> [1] FALSE
```
