# Angle between two vectors via the dot product

Computes the angle in \\\[0, \pi\]\\ between two numeric vectors using
the dot-product/`acos` formula, with the dot product clamped to \\\[-1,
1\]\\ for numerical stability near colinear and anti-colinear
configurations.

## Usage

``` r
angle_between(v1, v2)
```

## Arguments

- v1:

  Numeric vector. Unit-normalized inputs are recommended for speed;
  otherwise the caller is responsible for prescaling.

- v2:

  Numeric vector of the same length as `v1`.

## Value

A single numeric value in \\\[0, \pi\]\\: the angle between `v1` and
`v2` in radians.

## Details

For unit vectors the formula is \$\$\theta = \arccos(v_1 \cdot v_2);\$\$
for general vectors, divide the dot product by the product of the norms
before applying `acos`. The implementation assumes pre-normalised input
and only clamps the dot product to the unit interval.

## References

Strang, G. (2016). *Introduction to Linear Algebra*, 5th edition.
Wellesley-Cambridge Press. ISBN 978-0980232776.

## See also

[`lhuilier_angle`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)
for the spherical excess of a triangle whose sides are returned by this
function;
[`spherical_triangle_area`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md)
for the Van Oosterom-Strackee form;
[`cross_product_3d`](https://robustecologies.github.io/SolidAngleR/reference/cross_product_3d.md)
for the associated 3D vector product.

## Examples

``` r
angle_between(c(1, 0, 0), c(0, 1, 0))                # pi / 2
#> [1] 1.570796
angle_between(c(1, 1, 1) / sqrt(3), c(2, 2, 2) / sqrt(12))  # ~ 0
#> [1] 0
angle_between(c(1, 0, 0), c(-1, 0, 0))               # pi
#> [1] 3.141593
```
