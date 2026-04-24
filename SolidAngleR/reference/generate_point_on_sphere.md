# Uniform random point on the n-dimensional unit sphere

Returns a single point drawn uniformly from the surface of the
\\n\\-dimensional unit sphere via the Box-Muller / Marsaglia method:
draw \\n\\ independent standard Gaussians and normalize the resulting
vector to unit Euclidean norm.

## Usage

``` r
generate_point_on_sphere(n)
```

## Arguments

- n:

  Integer. Dimension of the sphere.

## Value

Numeric vector of length `n` representing a point on the unit sphere.

## Details

This function uses the Box-Muller transform: if \\Z_1, Z_2, \ldots,
Z_n\\ are standard normal random variables and \\\hat{e}\_1, \hat{e}\_2,
\ldots, \hat{e}\_n\\ are canonical basis vectors, then the uniformly
distributed vector \\\hat{s}\\ on the sphere surface is: \$\$\hat{s} =
\frac{Z_1\hat{e}\_1 + Z_2\hat{e}\_2 + \cdots + Z_n\hat{e}\_n}
{\sqrt{Z_1^2 + Z_2^2 + \cdots + Z_n^2}}\$\$

## References

Box, G. E. P., & Muller, M. E. (1958). A note on the generation of
random normal deviates. *Annals of Mathematical Statistics*, 29(2),
610-611.
[doi:10.1214/aoms/1177706645](https://doi.org/10.1214/aoms/1177706645)

Marsaglia, G. (1972). Choosing a point from the surface of a sphere.
*Annals of Mathematical Statistics*, 43(2), 645-646.
[doi:10.1214/aoms/1177692644](https://doi.org/10.1214/aoms/1177692644)

## See also

[`generate_cone_sample`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
[`generate_cone_samples`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
for cone-restricted sampling that builds on this primitive;
[`rotate_from_canonical`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md)
for the cheap O(n) rotation used to align cone samples with a target
axis;
[`verify_cone_uniformity`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
for goodness-of-fit testing of generated samples.

## Examples

``` r
point <- generate_point_on_sphere(3)
sqrt(sum(point^2))   # equals 1 up to floating-point error
#> [1] 1
```
