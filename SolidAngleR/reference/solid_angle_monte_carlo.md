# Generic Monte Carlo estimator of a solid angle

Estimates the unnormalized solid angle (in steradians) of an arbitrary
region on the unit sphere by uniformly sampling points and counting
those satisfying the user-supplied membership test. Returns the point
estimate together with a binomial standard error.

## Usage

``` r
solid_angle_monte_carlo(region_test, n_samples = 1e+05, seed = NULL)
```

## Arguments

- region_test:

  Function that takes a 3D point and returns TRUE if inside region

- n_samples:

  Number of Monte Carlo samples (default: 100000)

- seed:

  Random seed for reproducibility (optional)

## Value

List containing:

- estimate:

  Estimated solid angle in steradians

- std_error:

  Standard error of the estimate

- n_samples:

  Number of samples used

- fraction_inside:

  Fraction of samples inside the region

## Details

This function implements Monte Carlo integration on the unit sphere by:

1.  Generating n_samples points uniformly distributed on the unit sphere

2.  Testing each point against the region criterion

3.  Estimating solid angle as: \\\Omega = 4\pi \times (\text{fraction
    inside})\\

4.  Computing standard error using binomial statistics: \\\sigma =
    4\pi\sqrt{p(1-p)/n}\\ where p is the fraction inside

Uniform sampling on the sphere is achieved using the Marsaglia method:
generate 3D Gaussian random vectors and normalize to unit length. This
ensures equal probability density over the sphere surface.

The convergence rate is O(1/√n), which is typical for Monte Carlo
methods. For 1% relative accuracy, approximately n = 10,000 samples are
needed. For 0.1% accuracy, n = 1,000,000 samples are recommended.

This method is particularly useful when the region has complex geometry
not amenable to analytical formulas; when verification of analytical
results is needed; or when approximate solutions with known uncertainty
bounds are acceptable.

## References

Marsaglia, G. (1972). Choosing a point from the surface of a sphere.
*Annals of Mathematical Statistics*, 43(2), 645-646.
[doi:10.1214/aoms/1177692644](https://doi.org/10.1214/aoms/1177692644)

Arvo, J. (2001). Stratified sampling of spherical triangles.
*Proceedings of ACM SIGGRAPH 2001*, 437-438.
[doi:10.1145/383259.383300](https://doi.org/10.1145/383259.383300)

Pharr, M., Jakob, W., & Humphreys, G. (2016). *Physically based
rendering: from theory to implementation*, 3rd edition. Morgan Kaufmann.
ISBN 978-0128006450.

## See also

[`solid_angle_cone`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone.md),
[`solid_angle_cone_segment`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md),
[`solid_angle_intersecting_cones`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
for analytic backends that avoid Monte Carlo error;
[`generate_point_on_sphere`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md),
the single-sample analogue of the uniform-sphere sampling step;
[`verify_cone_uniformity`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
for goodness-of-fit testing on the samples;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
# Estimate solid angle of a cone
cone_test <- function(point) {
  # Test if point is in cone with axis along z and apex angle pi/3
  z_component <- point[3]
  return(z_component >= cos(pi/3))
}
result <- solid_angle_monte_carlo(cone_test, n_samples = 10000)
cat("Estimated solid angle:", result$estimate, "±", result$std_error)
#> Estimated solid angle: 3.124 ± 0.05431203
```
