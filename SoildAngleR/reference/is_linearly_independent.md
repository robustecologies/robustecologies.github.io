# Check if vectors are linearly independent

Determines if a set of vectors is linearly independent by checking the
rank of the matrix formed by these vectors.

## Usage

``` r
is_linearly_independent(V, tol = 1e-10)
```

## Arguments

- V:

  A matrix where columns are vectors to check

- tol:

  Tolerance for rank determination (default: 1e-10)

## Value

Logical: TRUE if linearly independent, FALSE otherwise

## Examples

``` r
# Linearly independent
V1 <- cbind(c(1, 0), c(0, 1))
is_linearly_independent(V1)  # TRUE
#> [1] TRUE

# Linearly dependent
V2 <- cbind(c(1, 0), c(2, 0))
is_linearly_independent(V2)  # FALSE
#> [1] FALSE
```
