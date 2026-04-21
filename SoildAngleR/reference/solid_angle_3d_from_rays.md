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

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Example 1: A cone covering one orthant of 3D space ---
# The rays are the positive axes. The solid angle should be (4*pi)/8 = pi/2.
orthant_rays <- matrix(c(1, 0, 0,
                         0, 1, 0,
                         0, 0, 1), nrow = 3)
angle <- solid_angle_3d_from_rays(orthant_rays)
cat("Solid angle of one orthant:", angle, "\n")
cat("Expected value (pi/2):", pi/2, "\n")

# --- Example 2: A wider cone with 4 ordered rays ---
# Rays form a square base on the sphere's surface.
ray1 <- c(1, 0.5, 0.5)
ray2 <- c(1, -0.5, 0.5)
ray3 <- c(1, -0.5, -0.5)
ray4 <- c(1, 0.5, -0.5)
wide_cone_rays <- cbind(ray1, ray2, ray3, ray4)

wide_angle <- solid_angle_3d_from_rays(wide_cone_rays)
cat("\nSolid angle of a wider 4-ray cone:", wide_angle, "steradians\n")
} # }
```
