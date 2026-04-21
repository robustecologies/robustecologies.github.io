# Rotate vector from n-th canonical basis to arbitrary orientation

Rotates a vector \\\hat{x}\\ aligned with the n-th canonical basis
vector to align with an arbitrary direction \\\hat{\mu}\\.

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

## Examples

``` r
# Rotate vector from z-axis to arbitrary direction
x <- c(0, 0, 1)  # Aligned with 3rd canonical axis
mu_hat <- c(1, 1, 1) / sqrt(3)  # Target direction
y <- rotate_from_canonical(x, mu_hat)

# Verify alignment
sum(y * mu_hat)  # Should be close to 1
#> [1] 1
```
