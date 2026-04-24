# Uniform random point on a spherical cap

Returns a single point drawn uniformly from the spherical cap of axis
\\\hat{\mu}\\ and maximum planar angle \\\theta_0\\ in dimension \\n\\,
using the O(n) algorithm of Arun & Venkatapathi (2025).

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

The spherical cap is the set of unit vectors \\\hat{x}\\ satisfying
\\\hat{x} \cdot \hat{\mu} \ge \cos\theta_0\\. The algorithm samples a
random planar angle \\\theta \in \[0, \theta_0\]\\ from the appropriate
marginal distribution
([`generate_planar_angle_inverse`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md)
or
[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
depending on `method`), draws a random point on the \\(n-1)\\-sphere for
the perpendicular component
([`generate_point_on_sphere`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md)),
composes the canonical representation \\\hat{x} = \sin\theta \cdot
\mathrm{perp} + \cos\theta \cdot \hat{e}\_n\\, and rotates the result to
align \\\hat{e}\_n\\ with \\\hat{\mu}\\ via
[`rotate_from_canonical`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md).
The total cost is \\O(n)\\ per sample.

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A: The Indian
Journal of Statistics*, 87(2), 327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
for the vectorised C++-accelerated batch sampler;
[`generate_hollow_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md)
for the hollow- cone variant;
[`verify_cone_uniformity`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
for goodness-of-fit tests on the resulting samples;
[`rejection_cost`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
for the expected per-sample cost of the rejection variant.

## Examples

``` r
# 30-degree cap around the z-axis in R^3
mu_hat <- c(0, 0, 1)
samples <- replicate(1000, generate_cone_sample(mu_hat, pi / 6))
summary(acos(samples[3, ]))     # cap-axis angles in [0, pi/6]
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0124  0.2631  0.3689  0.3498  0.4544  0.5233 

# High dimension: use the rejection planar-angle generator
mu_hat_100d <- c(rep(0, 99), 1)
generate_cone_sample(mu_hat_100d, pi / 12, method = "rejection")
#>   [1]  2.165654e-03  2.737157e-02  4.771550e-02 -3.246124e-03  4.470982e-02
#>   [6] -3.565361e-02  9.863948e-03  3.724912e-02 -9.478066e-03 -1.483450e-02
#>  [11] -5.580486e-04 -9.215287e-03  5.047322e-03  9.406004e-03 -2.157529e-02
#>  [16] -2.035962e-02 -5.022154e-02 -4.233905e-02  6.404445e-05  3.998341e-02
#>  [21] -2.114430e-02  1.916251e-02  1.894475e-02  2.534573e-02  7.114902e-02
#>  [26] -4.018155e-02  2.030905e-03 -1.742624e-02  3.125488e-03  2.031583e-02
#>  [31] -3.041791e-02  9.273892e-03  3.400321e-03  8.467942e-03  1.159747e-02
#>  [36] -3.452720e-02  2.365519e-02 -6.411304e-04  1.757033e-02  3.352350e-02
#>  [41] -3.482176e-02 -4.788801e-03  3.124557e-02  2.893887e-03  2.862626e-02
#>  [46]  1.727168e-03  2.067664e-02 -1.586578e-02 -2.265890e-02 -2.597956e-02
#>  [51]  1.273299e-02  4.050435e-02  2.299778e-02 -2.456604e-03 -2.069528e-02
#>  [56] -1.149575e-02 -3.037743e-02 -3.734539e-02 -1.263590e-02 -2.045075e-02
#>  [61] -1.109519e-02  4.818618e-03  2.717746e-02  3.162033e-03 -4.990603e-02
#>  [66] -8.123262e-03  8.557226e-03 -2.943851e-02  6.031834e-03 -1.355745e-03
#>  [71]  3.012027e-02 -2.701522e-02  5.234762e-02 -1.514159e-02  8.960320e-03
#>  [76]  4.781797e-02  3.464497e-02 -1.587921e-02 -8.911850e-03  3.949949e-03
#>  [81]  5.822647e-03  3.142382e-04  3.317067e-02 -2.802916e-02  2.079756e-02
#>  [86] -1.400535e-02  4.734750e-02  2.735475e-02 -3.011533e-02 -1.707576e-02
#>  [91] -4.023621e-02  4.197483e-03 -6.253627e-02 -1.193154e-02 -1.766857e-02
#>  [96]  1.153212e-02  1.279554e-03 -3.494777e-02  1.238288e-03  9.661854e-01
```
