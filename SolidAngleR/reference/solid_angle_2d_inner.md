# Normalized solid angle of a planar cone via inner product

Computes the normalized solid angle of a two-dimensional simplicial cone
using the dot-product / `acos` variant of the planar angle formula.
Provided as a complement to
[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
for numerical cross-checks.

## Usage

``` r
solid_angle_2d_inner(v1, v2)
```

## Arguments

- v1:

  Numeric vector of length two; first cone generator. Need not be
  unit-normalized.

- v2:

  Numeric vector of length two; second cone generator.

## Value

A single numeric value in \\\[0, 0.5\]\\: the planar opening angle
divided by \\2\pi\\. Colinear generators return `0`.

## Details

For unit vectors \\v_1, v_2 \in \mathbb{R}^2\\ the opening angle is
\$\$\theta = \arccos(v_1 \cdot v_2),\$\$ giving \\\Omega = \theta /
(2\pi)\\. The dot product is clamped to \\\[-1, 1\]\\ before the `acos`
call to prevent `NaN` from floating-point overshoot. The `acos`
formulation loses precision for vectors very close to colinear or
anti-colinear, so for production use
[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
(`atan2`-based) is preferred; this function is retained as a reference
implementation and for verification tests.

## References

Hanson, A. J. (2006). *Visualizing Quaternions*. Morgan Kaufmann. ISBN
978-0120884001. (Chapter on planar angle conventions.)

## See also

[`solid_angle_2d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_2d.md)
for the numerically preferred `atan2`-based variant;
[`solid_angle_3d`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d.md)
for three dimensions;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
# Right angle (quarter plane)
solid_angle_2d_inner(c(1, 0), c(0, 1))           # 0.25
#> [1] 0.25

# 45-degree opening: same numerical value as solid_angle_2d()
solid_angle_2d_inner(c(1, 0), c(1, 1) / sqrt(2)) # 0.125
#> [1] 0.125
```
