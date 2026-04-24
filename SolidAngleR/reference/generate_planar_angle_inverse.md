# Random planar angle by inverse-transform sampling

Draws a single random planar angle \\\theta \in \[0, \theta_0\]\\ from
the distribution required to produce uniform sampling on a spherical cap
of half-angle \\\theta_0\\ in dimension \\n\\, using the inverse
cumulative distribution map \\\theta = \Theta^{-1}(U)\\ with \\U \sim
\mathrm{Uniform}(0, \Omega_0)\\.

## Usage

``` r
generate_planar_angle_inverse(theta0, n)
```

## Arguments

- theta0:

  Numeric. Maximum planar angle (cone half-angle) in radians.

- n:

  Integer. Dimension of the space.

## Value

Numeric. Random angle \\\theta \in \[0, \theta_0\]\\.

## Details

The cumulative distribution function of \\\theta\\ restricted to \\\[0,
\theta_0\]\\ is \\F\_\theta(\theta) = \Theta(\theta) /
\Theta(\theta_0)\\, where \\\Theta\\ is implemented by
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md).
Inverse-transform sampling generates one random angle per call by
drawing \\U \sim \mathrm{Uniform}(0, \Omega_0)\\ with \\\Omega_0 =
\Theta(\theta_0)\\ and returning \\\Theta^{-1}(U)\\ through
[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md).

This formulation requires both \\\Theta\\ and \\\Theta^{-1}\\, both of
which can underflow for very large \\n\\ and very narrow caps; in that
regime
[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
is preferred.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
for the rejection-based alternative that is numerically stable in high
dimensions;
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md),
[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md)
for the underlying maps;
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the downstream cone sampler that consumes this random angle.

## Examples

``` r
# 30-degree cap in dimension 5: empirical density of planar angles
angles <- replicate(1000, generate_planar_angle_inverse(pi / 6, 5))
summary(angles)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1131  0.3568  0.4281  0.4109  0.4801  0.5235 
```
