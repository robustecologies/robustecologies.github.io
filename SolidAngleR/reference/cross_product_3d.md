# Cross product of two three-dimensional vectors

Returns the cross product \\a \times b\\ of two three-dimensional
numeric vectors. Provided as an explicit helper because base R has no
built-in cross product and inlining it in every consumer obscures the
geometry.

## Usage

``` r
cross_product_3d(a, b)
```

## Arguments

- a:

  Numeric vector of length three.

- b:

  Numeric vector of length three.

## Value

A numeric vector of length three: the cross product \\a \times b\\.

## Details

By construction \\a \times b\\ is orthogonal to both \\a\\ and \\b\\;
its norm equals \\\\a\\\\\\b\\\\\sin\theta\\ where \\\theta\\ is the
angle between the inputs
([`angle_between`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)).
The scalar triple product \\a \cdot (b \times c)\\ that appears in the
Van Oosterom-Strackee formula
([`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md))
is the signed volume of the parallelepiped spanned by \\a, b, c\\.

## References

Strang, G. (2016). *Introduction to Linear Algebra*, 5th edition.
Wellesley-Cambridge Press. ISBN 978-0980232776.

## See also

[`angle_between`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)
for the related angle helper;
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md),
[`spherical_triangle_area`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md)
for the consumers in the 3D solid-angle formulas.

## Examples

``` r
cross_product_3d(c(1, 0, 0), c(0, 1, 0))   # c(0, 0, 1)
#> [1] 0 0 1
```
