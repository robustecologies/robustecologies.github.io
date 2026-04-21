# Compute solid angle of arbitrary conical surface

Calculates the solid angle of a conical surface defined by a closed
parametric curve on the unit sphere. Currently implemented as a
polyhedral approximation.

## Usage

``` r
solid_angle_conical_surface(curve_points, corner_angles, smooth_segments)
```

## Arguments

- curve_points:

  Matrix where each row represents a point on the unit sphere (n x 3)

- corner_angles:

  Vector of turn angles at corners (in radians). Length should match
  number of corners. Currently not used in the polyhedral approximation.

- smooth_segments:

  Logical vector indicating which segments are smooth (length n-1).
  Currently not used in the polyhedral approximation.

## Value

Normalized solid angle (0 to 1, where 1 = full sphere)

## Details

According to Mazonka (2012), the solid angle of a general conical
surface bounded by a closed curve \\\mathcal{C}\\ on the unit sphere is
given by Equation 3:

\$\$\Omega = \int\_{\mathcal{C}} \kappa \\ ds + \sum\_{i}
\Delta\varphi_i\$\$

where \\\kappa\\ is the geodesic curvature and \\\Delta\varphi_i\\ are
the exterior angles at corners.

**Current implementation**: This function currently uses a polyhedral
approximation via
[`solid_angle_polyhedral`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md),
treating the curve points as vertices of a polygon. This provides
accurate results when the curve is already piecewise linear, but is an
approximation for smooth curves.

Future versions may implement the full geodesic curvature integral for
smooth conical surfaces.

## References

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396v2. Section 2: Conical
surfaces.
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## Examples

``` r
# Simple circular cone approximated by polygon
n_points <- 8
theta <- pi/3  # 60 degree cone
angles <- seq(0, 2*pi, length.out = n_points + 1)[1:n_points]
curve_points <- cbind(
  sin(theta) * cos(angles),
  sin(theta) * sin(angles),
  rep(cos(theta), n_points)
)
smooth_segments <- rep(TRUE, n_points)
corner_angles <- rep(2*pi/n_points, n_points)
solid_angle_conical_surface(curve_points, corner_angles, smooth_segments)
#> [1] 0.2399796

if (FALSE) { # \dontrun{
# Compare polygonal approximation convergence
theta <- pi/4
n_vals <- c(4, 8, 16, 32, 64)
omega_approx <- numeric(length(n_vals))

for (i in seq_along(n_vals)) {
  n <- n_vals[i]
  angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
  curve_points <- cbind(
    sin(theta) * cos(angles),
    sin(theta) * sin(angles),
    rep(cos(theta), n)
  )
  omega_approx[i] <- solid_angle_conical_surface(
    curve_points,
    rep(2*pi/n, n),
    rep(TRUE, n)
  )
}

plot(n_vals, omega_approx,
     type = "b", log = "x",
     xlab = "Number of polygon sides",
     ylab = "Solid angle",
     main = "Convergence of polygonal approximation")
abline(h = solid_angle_cone(theta), col = "red", lty = 2)
legend("bottomright", c("Approximation", "Exact"),
       pch = c(1, NA), lty = c(1, 2), col = c("black", "red"))
} # }
```
