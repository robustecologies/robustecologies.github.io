# Compute solid angle using Ribando's hypergeometric series

Implements the multivariable hypergeometric series formula from Ribando
(2006) and Theorem 1.5 of Fitisone & Zhou (2023) to compute the
normalized solid angle measure of a simplicial cone.

## Usage

``` r
hypergeometric_series(V, max_terms = 1000, tol = 1e-10, check_pd = TRUE)
```

## Arguments

- V:

  An n x n matrix where columns are unit vectors generating the cone

- max_terms:

  Maximum number of terms to compute (default: 1000)

- tol:

  Convergence tolerance (default: 1e-10)

- check_pd:

  Whether to check positive definiteness (default: TRUE)

## Value

A list containing:

- `solid_angle`: The normalized solid angle measure

- `n_terms`: Number of terms computed

- `converged`: Logical indicating convergence

- `associated_matrix`: The associated matrix M_n(C)

## Details

The series formula (equation 5 in the paper) is:

\$\$T\_{\alpha} = \frac{\|det V\|}{(4\pi)^{n/2}} \sum \left\[
\frac{(-2)^{\sum a\_{ij}}}{\prod a\_{ij}!} \prod \Gamma\left(\frac{1 +
\sum\_{m \neq i} a\_{im}}{2}\right) \right\] \alpha^a\$\$

This series converges absolutely to the solid angle measure if and only
if the associated matrix \\M_n(C)\\ is positive definite.

## References

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

Aomoto, K. (1977). Analytic structure of Schlafli function. *Nagoya
Mathematical Journal*, 68, 1-16.
<https://projecteuclid.org/euclid.nmj/1118786429>

## Examples

``` r
if (FALSE) { # \dontrun{
# 3D orthogonal case (should give 1/8)
V <- diag(3)
result <- hypergeometric_series(V, max_terms = 100)
print(result$solid_angle)  # 0.125

# 4D orthogonal case (should give 1/16)
V <- diag(4)
result <- hypergeometric_series(V, max_terms = 200)
print(result$solid_angle)  # 0.0625
} # }
```
