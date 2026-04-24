# Calculate solid angle of a 3D cone from extreme rays

Calculates the solid angle of a 3D cone deterministically using its
extreme rays (V-representation). The method triangulates the spherical
polygon formed by the rays on the unit sphere and sums the areas of the
constituent spherical triangles using L'Huilier's theorem.

## Usage

``` r
solid_angle_3d_from_rays(rays)
```

## Arguments

- rays:

  A numeric matrix with 3 rows, where each column represents a distinct
  extreme ray of the cone. At least 3 rays are required.

## Value

A single numeric value representing the solid angle of the cone in
steradians.

## Details

The solid angle is a measure of the "visual size" of an object from a
point. This function provides a precise value in steradians, avoiding
the sampling error inherent in Monte Carlo methods. The total solid
angle of a full sphere is 4\*pi steradians.

## Note

CRITICAL: This function assumes the input `rays` are ordered
sequentially (e.g., counter-clockwise) around the cone's central axis.
If the rays are unordered, the triangulation will be incorrect and the
result will be invalid. Vertex enumeration algorithms (like in `rcdd`)
do not guarantee such an ordering.

## References

Malekmohammadi, N., & Mostafaee, A. (2017). Obtaining all extreme rays
of a special cone using spanning trees in a complete digraph:
application in DEA. *Journal of the Operational Research Society*,
69(3), 465-472.
[doi:10.1057/s41274-017-0265-9](https://doi.org/10.1057/s41274-017-0265-9)

Todhunter, I. (1886). *Spherical Trigonometry*, 5th edition. Macmillan
and Co., London.

## See also

[`generate_spanning_trees`](https://robustecologies.github.io/SolidAngleR/reference/generate_spanning_trees.md)
for the enumeration of extreme rays that this function consumes;
[`lhuilier_angle`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)
for the per-triangle building block;
[`angle_between`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md)
for the arc-length helper;
[`solid_angle_polyhedral`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md)
for an alternative polyhedral computation;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
# Octant in R^3: solid angle should be pi / 2 steradians
orthant <- matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1), nrow = 3)
solid_angle_3d_from_rays(orthant)        # ~ 1.5708

# Square-base wider cone with four ordered rays
rays <- cbind(c(1,  0.5,  0.5),
              c(1, -0.5,  0.5),
              c(1, -0.5, -0.5),
              c(1,  0.5, -0.5))
solid_angle_3d_from_rays(rays)
} # }
```
