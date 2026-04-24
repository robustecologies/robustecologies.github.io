# Givens rotation from canonical axis to an arbitrary direction

Rotates a vector \\\hat{x}\\ aligned with the \\n\\-th canonical basis
vector \\\hat{e}\_n\\ so that its image is aligned with an arbitrary
unit direction \\\hat{\mu}\\, in \\O(n)\\ operations using a Givens
rotation in the plane spanned by \\\hat{e}\_n\\ and \\\hat{\mu}\\.

## Usage

``` r
rotate_from_canonical(x, mu_hat)
```

## Arguments

- x:

  Numeric vector. Vector to rotate (should be aligned with n-th
  canonical axis).

- mu_hat:

  Numeric vector. Target direction (unit vector).

## Value

Numeric vector. Rotated vector \\\hat{y}\\.

## Details

This function performs a Givens rotation in the plane containing
\\\hat{e}\_n\\ and \\\hat{\mu}\\. The rotation is: \$\$\hat{y} =
\hat{x} + P(G - I_2)P^T\hat{x}\$\$ where P is an orthonormal basis for
the plane of rotation, G is the 2D Givens rotation matrix, and \\I_2\\
is the 2×2 identity matrix.

The matrices are defined as: \$\$P = \[\hat{e}\_n, \frac{\hat{\mu} -
\mu_n\hat{e}\_n} {\\\hat{\mu} - \mu_n\hat{e}\_n\\}\]\$\$ \$\$G =
\begin{bmatrix} \mu_n & -\sqrt{1-\mu_n^2} \\ \sqrt{1-\mu_n^2} & \mu_n
\end{bmatrix}\$\$ where \\\mu_n = \hat{e}\_n^T\hat{\mu}\\.

This rotation costs only O(n) operations, unlike general n-dimensional
rotations which cost O(n²).

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
for the cone samplers that compose this rotation with a sampled
canonical-axis vector;
[`generate_point_on_sphere`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md)
for the uniform-sphere primitive.

## Examples

``` r
x <- c(0, 0, 1)
mu_hat <- c(1, 1, 1) / sqrt(3)
y <- rotate_from_canonical(x, mu_hat)
sum(y * mu_hat)             # close to 1: y aligns with mu_hat
#> [1] 1
```
