# Compute solid angle for spherical triangle

Alternative representation: the solid angle is the area of the spherical
triangle on the unit sphere formed by the three unit vectors.

## Usage

``` r
spherical_triangle_area(v1, v2, v3)
```

## Arguments

- v1:

  First unit vector

- v2:

  Second unit vector

- v3:

  Third unit vector

## Value

Area of spherical triangle (unnormalized solid angle E)

## Details

This returns E (the actual solid angle in steradians), not the
normalized version. To normalize, divide by \\4\pi\\.

## Examples

``` r
# Octant
v1 <- c(1, 0, 0)
v2 <- c(0, 1, 0)
v3 <- c(0, 0, 1)
E <- spherical_triangle_area(v1, v2, v3)
print(E)  # pi/2 steradians
#> [1] 1.570796
print(E / (4 * pi))  # 0.125 (normalized)
#> [1] 0.125
```
