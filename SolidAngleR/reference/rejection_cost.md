# Expected sample cost of naive sphere-rejection sampling on a cap

Returns the expected number of full-sphere samples needed to obtain one
sample inside the spherical cap of half-angle \\\theta_0\\ in dimension
\\n\\, namely the reciprocal of the cap's normalised solid angle.

## Usage

``` r
rejection_cost(theta0, n)
```

## Arguments

- theta0:

  Numeric. Maximum planar angle of the cone in radians.

- n:

  Integer. Ambient dimension of the sphere (\\n \geq 2\\).

## Value

A single positive numeric value: \\1 / \Theta(\theta_0)\\ where
\\\Theta\\ is the cap-area function implemented by
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md).
Equal to \\1\\ for the full sphere (\\\theta_0 = \pi\\) and growing
rapidly as the cap narrows or the ambient dimension increases.

## Details

Naive rejection sampling on the full sphere is the baseline algorithm
that the O(n) sampler of Arun & Venkatapathi (2025) replaces. The growth
of the rejection cost with dimension is the explicit motivation for the
structured sampler in
[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md).

## References

Arun, I., & Venkatapathi, M. (2025). An O(n) algorithm for generating
uniform random vectors in n-dimensional cones. *Sankhya A*, 87(2),
327-348.
[doi:10.1007/s13171-025-00387-9](https://doi.org/10.1007/s13171-025-00387-9)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
for the structured sampler that avoids this cost;
[`generate_planar_angle_rejection`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
for the one-dimensional rejection step actually used inside the
structured sampler;
[`theta_to_omega`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md)
for the underlying cap-area function.

## Examples

``` r
dims <- c(2, 5, 10, 20, 50, 100)
rates <- sapply(dims, function(n) rejection_cost(pi / 6, n))
rbind(dimension = dims, log10_cost = log10(rates))
#>                 [,1]     [,2]      [,3]      [,4]    [,5]      [,6]
#> dimension  2.0000000 5.000000 10.000000 20.000000 50.0000 100.00000
#> log10_cost 0.7781513 1.890735  3.547474  6.707942 15.9372  31.13893
```
