# Normalized solid angle of a right circular cone

Returns the normalized solid angle of a right circular cone with apex
half-angle \\\theta\\ via the closed-form Mazonka (2012, Equation 9)
expression.

## Usage

``` r
solid_angle_cone(theta)
```

## Arguments

- theta:

  Apex **half**-angle of the cone in radians (0 \< theta \< pi)

## Value

Normalized solid angle (0 to 1, where 1 = full sphere)

## Details

A right circular cone with apex at the origin and apex half-angle
\\\theta\\ subtends a solid angle given by:

\$\$\Omega = 2\pi(1 - \cos\theta)\$\$

This implementation returns the normalized value \\\Omega/(4\pi)\\,
which represents the fraction of the full sphere (Equation 9 in Mazonka
2012).

For special cases: when \\\theta = \pi/2\\ (hemisphere), the normalized
solid angle equals 0.5; as \\\theta \to 0\\ (narrow cone), it approaches
0; and as \\\theta \to \pi\\ (full sphere), it approaches 1.

## References

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396 (math.MG). Equation 9.
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## See also

[`solid_angle_cone_segment`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
for a cone cut by a plane;
[`solid_angle_intersecting_cones`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
for the intersection of two cones;
[`solid_angle_polyhedral`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md)
for spherical polygon regions;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher.

## Examples

``` r
# Solid angle of a 60-degree cone
solid_angle_cone(pi/3)
#> [1] 0.25

# Solid angle of a hemisphere (90-degree cone) - should be 0.5
solid_angle_cone(pi/2)
#> [1] 0.5

if (FALSE) { # \dontrun{
# Verify the formula for various cone angles
theta_values <- seq(0.1, pi - 0.1, length.out = 20)
omega_values <- sapply(theta_values, solid_angle_cone)
plot(theta_values, omega_values, type = "l",
     xlab = "Apex angle (radians)", ylab = "Normalized solid angle",
     main = "Solid angle vs. apex angle for circular cones")
abline(h = 0.5, v = pi/2, col = "red", lty = 2)
} # }
```
