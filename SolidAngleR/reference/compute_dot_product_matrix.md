# Compute dot product matrix

Computes the matrix of all pairwise dot products between column vectors.

## Usage

``` r
compute_dot_product_matrix(V)
```

## Arguments

- V:

  A matrix where columns are vectors

## Value

A symmetric matrix where entry (i,j) is \\v_i \cdot v_j\\

## Examples

``` r
V <- cbind(c(1, 0), c(0, 1), c(1, 1) / sqrt(2))
alpha <- compute_dot_product_matrix(V)
print(alpha)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.0000000 0.7071068
#> [2,] 0.0000000 1.0000000 0.7071068
#> [3,] 0.7071068 0.7071068 1.0000000
```
