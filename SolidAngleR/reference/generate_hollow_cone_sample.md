# Uniform random point on a hollow spherical cap

Returns a single point drawn uniformly from the hollow spherical cap of
axis \\\hat{\mu}\\ and inner / outer planar angles \\\theta_1 \<
\theta_2\\ in dimension \\n\\.

## Usage

``` r
generate_hollow_cone_sample(
  mu_hat,
  theta1,
  theta2,
  method = c("inverse", "rejection")
)
```

## Arguments

- mu_hat:

  Numeric vector. Central axis direction (will be normalized to unit
  vector).

- theta1:

  Numeric. Inner planar angle in radians, \\0 \le \theta_1 \<
  \theta_2\\.

- theta2:

  Numeric. Outer planar angle in radians, \\0 \< \theta_2 \le \pi\\.

- method:

  Character. Method for generating planar angle: `"inverse"` for inverse
  transform sampling (default), or `"rejection"` for rejection sampling.

## Value

Numeric vector. Unit vector uniformly distributed on the hollow cone
surface.

## Details

The hollow cone is defined as the set of all unit vectors \\\hat{x}\\
satisfying: \$\$\cos\theta_2 \le \hat{x} \cdot \hat{\mu} \le
\cos\theta_1\$\$

The algorithm generates a random angle \\\theta\\ from: \$\$\theta =
\Theta^{-1}(U(\Omega_2 - \Omega_1) + \Omega_1)\$\$ where \\U \sim
\text{Uniform}(0,1)\\, \\\Omega_1 = \Theta(\theta_1)\\, and \\\Omega_2 =
\Theta(\theta_2)\\.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md)
for the solid-cap variant;
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md),
[`omega_to_theta`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md)
for the angle/area maps used internally;
[`verify_cone_uniformity`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
for sample-level diagnostics.

## Examples

``` r
mu_hat <- c(1, 0, 0)
samples <- replicate(1000, generate_hollow_cone_sample(mu_hat, pi / 6, pi / 3))
range(acos(samples[1, ]))   # angles fall in [pi/6, pi/3]
#> [1] 0.5261712 1.0466128
```
