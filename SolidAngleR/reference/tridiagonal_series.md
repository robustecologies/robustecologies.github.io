# Compute solid angle using tridiagonal hypergeometric series

Implements the simplified hypergeometric series formula (equation 23)
for cones whose associated matrices have tridiagonal structure. This
reduces the number of coordinates from n(n-1)/2 to n-1.

## Usage

``` r
tridiagonal_series(V, max_terms = 1000, tol = 1e-10)
```

## Arguments

- V:

  An n x n matrix where columns are unit vectors

- max_terms:

  Maximum number of terms (default: 1000)

- tol:

  Convergence tolerance (default: 1e-10)

## Value

A list containing:

- `solid_angle`: Normalized solid angle measure

- `n_terms`: Number of terms computed

- `converged`: Logical indicating convergence

- `beta`: Vector of consecutive dot products

- `is_tridiagonal`: Whether V^T V is tridiagonal

## Details

The formula (equation 23) for tridiagonal \\V^T V\\ is:

\$\$T\_{\beta} = \frac{\|det V\|}{(4\pi)^{n/2}} \sum \left\[
\frac{(-2)^{\sum b_i}}{\prod b_i!} \Gamma\left(\frac{1+b_1}{2}\right)
\Gamma\left(\frac{1+b_1+b_2}{2}\right) \cdots
\Gamma\left(\frac{1+b\_{n-2}+b\_{n-1}}{2}\right)
\Gamma\left(\frac{1+b\_{n-1}}{2}\right) \right\] \beta^b\$\$

where \\\beta_i = v_i \cdot v\_{i+1}\\ for \\i = 1, \ldots, n-1\\.

By Theorem 4.1, if V^T V is tridiagonal, then the associated matrix is
automatically positive definite, so the series always converges.

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO), Theorem 4.1 and equation (23).
<https://arxiv.org/abs/2304.11102>

## Examples

``` r
# Create a tridiagonal example
# Vectors with only consecutive dot products non-zero
v1 <- c(1, 0, 0, 0)
v2 <- c(0.8, 0.6, 0, 0)  # Only v1 · v2 != 0
v3 <- c(0, 0.6, 0.8, 0)  # Only v2 · v3 != 0
v4 <- c(0, 0, 0.7, 0.7) / sqrt(0.98)  # Only v3 · v4 != 0

V <- cbind(v1, v2, v3, v4)
result <- tridiagonal_series(V, max_terms = 500)
#> Warning: ⚠ Reached max_terms without convergence
print(result$solid_angle)
#>          i 
#> 0.01231954 
print(result$is_tridiagonal)  # Should be TRUE
#> [1] TRUE
```
