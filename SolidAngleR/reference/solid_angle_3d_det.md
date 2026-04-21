# Compute solid angle in R^3 using alternative parameterization

Computes solid angle using the determinant form, which is more
numerically stable in some cases.

## Usage

``` r
solid_angle_3d_det(V)
```

## Arguments

- V:

  A 3 x 3 matrix where columns are the unit vectors

## Value

Normalized solid angle measure

## Examples

``` r
V <- cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
solid_angle_3d_det(V)  # 0.125 (orthogonal)
#> [1] 0.125
```
