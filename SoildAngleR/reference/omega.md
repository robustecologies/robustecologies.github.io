# Compute normalized solid angle via multivariate normal integration

This function calculates the normalized solid angle of a polyhedral cone
using Ribando's (2006) approach via integration of a multivariate normal
distribution over the positive orthant.

## Usage

``` r
omega(alpha)
```

## Arguments

- alpha:

  A square invertible numeric matrix whose columns define the cone
  generators. The matrix must be square so that \\B = \alpha^T \alpha\\
  is positive definite and invertible.

## Value

A list containing two numeric values:

- `omega`: The raw normalized solid angle measure

- `omega_scaled`: The geometric mean per dimension
  (\\\omega^{\frac{1}{S}}\\)

## Details

**Mathematical foundation (Theorem 2.2 in Ribando 2006):**

The normalized solid angle measure is defined as:

\$\$\tilde{V}\_{\omega} = \frac{\int\_{\omega} f \\
dx}{\int\_{\mathbb{R}^n} f \\ dx}\$\$

where \\f = e^{-r^2}\\ and \\r = \|x\|\\. Through a change of
coordinates and application of the Gaussian integral, this reduces to
integrating a zero-mean multivariate normal distribution with covariance
matrix \\\Sigma = B^{-1}\\ over the positive orthant \\\mathbb{R}\_{\geq
0}^n\\.

**Characteristics of the algorithm:** The algorithm proceeds as follows:

1.  Compute \\B = \alpha^T \alpha\\

2.  Compute \\\Sigma = B^{-1}\\ (requires B to be positive definite)

3.  Integrate \\N(0, \Sigma)\\ over \\\[0, \infty)^n\\ with
    [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html),
    which uses the Genz-Bretz algorithm (Genz & Bretz 2009).

The Genz-Bretz algorithm implements Quasi-Monte Carlo integration with
Halton sequences, uses importance sampling to reduce variance, provides
adaptive subdivision for high-dimensional integrals, and estimates error
via repeated evaluations.

**Convergence:** This algorithm works for any cone where the associated
matrix B is positive definite. It provides an exact (up to numerical
integration tolerance) result without requiring series convergence.
However, if the matrix is ill-conditioned and/or matrix dimension is
large (\>20) results may be unreliable.

**Complexity**

- Matrix inversion: \\O(n^3)\\

- [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html):
  \\O(M \times n^2)\\ where M is number of integration points

- Total: \\O(n^3 + M \times n^2)\\

## References

Genz, A., & Bretz, F. (2009). *Computation of multivariate normal and t
probabilities*. Lecture Notes in Statistics, Vol. 195. Springer-Verlag,
Heidelberg.
[doi:10.1007/978-3-642-01689-9](https://doi.org/10.1007/978-3-642-01689-9)

Ribando, J. M. (2006). Measuring solid angles beyond dimension three.
*Discrete & Computational Geometry*, 36(3), 479-487.
[doi:10.1007/s00454-006-1253-4](https://doi.org/10.1007/s00454-006-1253-4)

## See also

[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for alternative methods including hypergeometric series and
decomposition approaches.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Orthogonal cone in R^3 (octant)
alpha <- diag(3)
result <- omega(alpha)
result$omega  # Should be 0.125 (1/8 of space)
} # }
```
