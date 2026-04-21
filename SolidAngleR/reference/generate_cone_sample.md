# Generate uniform random point on spherical cap

Generates a random point uniformly distributed on a spherical cap with
axis \\\hat{\mu}\\ and maximum planar angle \\\theta_0\\.

## Usage

``` r
generate_cone_sample(mu_hat, theta0, method = c("inverse", "rejection"))
```

## Arguments

- mu_hat:

  Numeric vector. Central axis direction (will be normalized to unit
  vector).

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians, \\0 \<
  \theta_0 \le \pi\\.

- method:

  Character. Method for generating planar angle: `"inverse"` for inverse
  transform sampling (default), or `"rejection"` for rejection sampling.
  Use `"rejection"` for numerical stability with large n and small
  \\\theta_0\\.

## Value

Numeric vector. Unit vector uniformly distributed on the spherical cap.

## Details

This is the main function implementing the O(n) algorithm. The spherical
cap is defined as the set of all unit vectors \\\hat{x}\\ that satisfy:
\$\$\hat{x} \cdot \hat{\mu} \ge \cos\theta_0\$\$

The algorithm proceeds as follows:

1.  Generate random angle \\\theta \in \[0, \theta_0\]\\ from the
    appropriate distribution

2.  Generate random point on (n-1)-dimensional sphere for the
    perpendicular component

3.  Construct vector aligned with n-th canonical axis: \\\hat{x} =
    \sin\theta \cdot \text{perp} + \cos\theta \cdot \hat{e}\_n\\

4.  Rotate to align with \\\hat{\mu}\\

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A: The Indian
Journal of Statistics*, 87(2), 327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## Examples

``` r
# Generate points on 3D spherical cap (30-degree cone)
mu_hat <- c(0, 0, 1)  # z-axis
points <- replicate(1000, generate_cone_sample(mu_hat, pi/6))

# Verify uniform distribution
angles <- acos(points[3, ])  # Angle from z-axis
hist(angles, main = "Distribution of angles from cone axis")


# High-dimensional example (use rejection method for stability)
mu_hat_100d <- c(rep(0, 99), 1)
point <- generate_cone_sample(mu_hat_100d, pi/12, method = "rejection")
```
