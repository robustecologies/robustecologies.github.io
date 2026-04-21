# Generate uniform random point on unit sphere

Generates a random point uniformly distributed on the surface of the
n-dimensional unit sphere using the Box-Muller transform.

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

## Examples

``` r
# Generate point on 3D sphere
point <- generate_point_on_sphere(3)

# Verify it's on unit sphere
sqrt(sum(point^2))  # Should be 1
#> [1] 1
```
