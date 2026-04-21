# Compute solid angle of intersection of two cones

Calculates the normalized solid angle of intersection between two cones
with the same apex by integrating the overlap of spherical caps on the
unit sphere.

## Usage

``` r
solid_angle_intersecting_cones(
  theta1,
  theta2,
  alpha,
  method = c("integral", "analytic", "mc"),
  fallback = c("integral", "mc", "none"),
  mc_samples = 2e+05,
  mc_seed = NULL
)
```

## Arguments

- theta1:

  Apex angle of first cone in radians (0 \< theta1 \< pi)

- theta2:

  Apex angle of second cone in radians (0 \< theta2 \< pi)

- alpha:

  Angle between cone axes in radians (0 \<= alpha \<= pi)

- method:

  Computation method: `"integral"` (default), `"analytic"`, or `"mc"`
  for Monte Carlo.

- fallback:

  Fallback method if analytic computation fails: `"integral"`, `"mc"`,
  or `"none"`.

- mc_samples:

  Number of Monte Carlo samples if `method = "mc"` or `fallback = "mc"`
  (default: 200000).

- mc_seed:

  Optional random seed for Monte Carlo.

## Value

Normalized solid angle of the intersection region (0 to 1, where 1 =
full sphere)

## Details

The intersection of two circular cones corresponds to the intersection
of two spherical caps on the unit sphere. In a coordinate system where
the first axis is the north pole and the second axis is separated by
angle \\\alpha\\, a point on the sphere is \$\$x(\theta,\varphi) =
(\sin\theta\cos\varphi,\\ \sin\theta\sin\varphi,\\ \cos\theta)\$\$ and
belongs to the second cone if \$\$\sin\alpha \sin\theta \cos\varphi +
\cos\alpha \cos\theta \ge \cos\theta_2.\$\$ For each \\\theta \in
\[0,\theta_1\]\\, the admissible azimuth range is \$\$\varphi \in
\[-\arccos C(\theta),\\ \arccos C(\theta)\]\$\$ where \$\$C(\theta) =
\frac{\cos\theta_2 - \cos\alpha \cos\theta}{\sin\alpha \sin\theta}.\$\$

The solid angle (steradians) of the intersection is then \$\$\Omega =
\int_0^{\theta_1} L(\theta)\\\sin\theta\\d\theta\$\$ with \\L(\theta)\\
the admissible azimuth length. The result is normalized by \\4\pi\\.

All methods handle several special cases: when there is no overlap (if
\\\alpha \ge \theta_1 + \theta_2\\), the function returns 0; when one
cone contains the other (if \\\alpha \le \|\theta_1 - \theta_2\|\\), it
returns the smaller cone; for co-directed axes (\\\alpha \approx 0\\),
it returns the smaller cone; for opposite axes (\\\alpha \approx \pi\\),
it uses a closed-form band area; when two hemispheres intersect, it
returns \\(\pi - \alpha)/(2\pi)\\; and when one cone is a hemisphere, it
reduces to
[`solid_angle_cone_segment`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
with \\\gamma = \alpha\\.

## References

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396v2. Section 4:
Intersecting cones.
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## Examples

``` r
# Two 60-degree cones with 30-degree separation
solid_angle_intersecting_cones(pi/3, pi/3, pi/6)
#> [1] 0.1781148

# Two hemispheres at 90 degrees
solid_angle_intersecting_cones(pi/2, pi/2, pi/2)
#> [1] 0.25

# Two identical cones with varying separation
if (FALSE) { # \dontrun{
theta <- pi/4  # 45-degree cones
alpha_vals <- seq(0, pi, length.out = 50)
omega_vals <- sapply(alpha_vals, function(a) {
  solid_angle_intersecting_cones(theta, theta, a)
})

plot(alpha_vals * 180/pi, omega_vals,
     type = "l", lwd = 2,
     xlab = "Separation angle (degrees)",
     ylab = "Normalized intersection solid angle",
     main = "Intersection of two 45-degree cones")
abline(h = solid_angle_cone(theta), col = "red", lty = 2)
legend("topright", c("Intersection", "Single cone"),
       lty = c(1, 2), col = c("black", "red"))
} # }
```
