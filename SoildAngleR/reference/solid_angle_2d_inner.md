# Compute solid angle for a 2D cone using inner product

Alternative method using the inner product to compute the solid angle.

## Usage

``` r
solid_angle_2d_inner(v1, v2)
```

## Arguments

- v1:

  First unit vector

- v2:

  Second unit vector

## Value

Normalized solid angle measure

## Examples

``` r
v1 <- c(1, 0)
v2 <- c(1, 1) / sqrt(2)
solid_angle_2d_inner(v1, v2)  # 0.125
#> [1] 0.125
```
