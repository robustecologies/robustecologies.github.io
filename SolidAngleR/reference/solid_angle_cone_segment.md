# Compute solid angle of a cone segment

Calculates the normalized solid angle of a cone segment cut by a plane
at angle gamma using Mazonka's analytical formula.

## Usage

``` r
solid_angle_cone_segment(theta, gamma)
```

## Arguments

- theta:

  Apex angle of the cone in radians (0 \< theta \< pi)

- gamma:

  Angle between the cutting plane normal and cone axis in radians (0 \<=
  gamma \<= pi)

## Value

Normalized solid angle (0 to 1, where 1 = full sphere)

## Details

A cone segment is the portion of a circular cone cut by a plane. The
solid angle depends on the cone's apex angle \\\theta\\ and the cutting
plane's orientation \\\gamma\\.

The calculation uses auxiliary angles (Equations 35, 39, 41 in Mazonka
2012):

\$\$\phi = \arccos\left(\frac{\cos\gamma}{\tan\theta}\right)\$\$
\$\$\beta = \arccos\left(\frac{\sin\gamma}{\sin\theta}\right)\$\$
\$\$\Omega = 2(\phi - \cos\theta \cdot \beta)\$\$ If the analytic
expression yields a negative or invalid value due to numerical issues or
boundary geometry, the function falls back to a robust numerical
integration of the equivalent spherical-cap intersection.

Special cases are handled automatically: when \\\gamma = 0\\, the
function returns the full cone solid angle; when \\\gamma \> \pi/2 +
\theta\\, the plane misses the cone entirely and the function returns 0;
and when \\\phi\\ or \\\beta\\ are invalid, geometric constraints
prevent intersection.

## References

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396v2. Section 4.1: Cone's
segment.
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## Examples

``` r
# Full cone (no cut)
solid_angle_cone_segment(pi/3, 0)
#> [1] 0.25

# Cone segment with 45-degree cut
solid_angle_cone_segment(pi/3, pi/4)
#> [1] 0.1340916

if (FALSE) { # \dontrun{
# Visualize how cutting plane angle affects solid angle
theta <- pi/3  # 60-degree cone
gamma_vals <- seq(0, pi/2, length.out = 50)
omega_vals <- sapply(gamma_vals, function(g) solid_angle_cone_segment(theta, g))

plot(gamma_vals * 180/pi, omega_vals,
     type = "l", lwd = 2,
     xlab = "Cutting plane angle (degrees)",
     ylab = "Normalized solid angle",
     main = "Cone segment solid angle vs. cutting plane angle")
abline(h = solid_angle_cone(theta), col = "red", lty = 2)
legend("topright", c("Segment", "Full cone"), lty = c(1, 2), col = c("black", "red"))
} # }
```
