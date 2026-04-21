# Create a tridiagonal cone from consecutive angles

Constructs a set of n unit vectors in R^n such that \\v_i \cdot v_j =
0\\ whenever \|i - j\| \> 1, resulting in a tridiagonal V^T V matrix.
This is done by building the desired tridiagonal Gram matrix and
applying a Cholesky factorization. If the Gram matrix is not positive
definite, the function stops with an error.

## Usage

``` r
create_tridiagonal_cone(angles)
```

## Arguments

- angles:

  Vector of n-1 angles (in radians) between consecutive vectors

## Value

An n x n matrix V with columns being the unit vectors

## Examples

``` r
# Create 4D tridiagonal cone with specified consecutive angles
angles <- c(pi/3, pi/3, pi/3)  # 60, 60, 60 degrees
V <- create_tridiagonal_cone(angles)

# Verify tridiagonality
VTV <- t(V) %*% V
is_tridiagonal(VTV)  # TRUE
#> [1] TRUE

# Compute solid angle
result <- tridiagonal_series(V)
#> Warning: ⚠ Reached max_terms without convergence
print(result$solid_angle)
#>          i 
#> 0.01726156 
```
