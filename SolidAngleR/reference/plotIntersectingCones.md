# Interactive 3D visualisation of two intersecting cones

Builds an interactive plotly widget showing two cones with common apex
on the unit sphere, their intersection circle and intersection points,
the spherical excess region, and the auxiliary angles (\\\alpha, \beta,
\gamma\\) that appear in the Mazonka (2012) analytic formula for the
solid angle of the intersection.

## Usage

``` r
plotIntersectingCones(
  theta1 = pi/6,
  theta2 = pi/4,
  phi = pi/3,
  radius = 1,
  color_palette = "viridis",
  show_rays = FALSE,
  cone_opacity = 0.2
)
```

## Arguments

- theta1:

  Apex angle of first cone in radians (default: pi/6)

- theta2:

  Apex angle of second cone in radians (default: pi/4)

- phi:

  Separation angle between cone axes in radians (default: pi/3)

- radius:

  Sphere radius (default: 1)

- color_palette:

  Color palette name: "viridis", "plasma", "inferno", "magma",
  "cividis", "rainbow" (default: "viridis")

- show_rays:

  Logical; if TRUE, show cone edge rays (default: FALSE)

- cone_opacity:

  Opacity of cone surfaces, range \$\[0, 1\]\$ (default: 0.2)

## Value

A plotly object with interactive 3D visualization

## Details

The visualization shows two cones with apex angles \\\theta_1\\ and
\\\theta_2\\, separated by angle \\\phi\\ between cone axes;
intersection points on the unit sphere are marked, along with great
circle arcs and planar circle arcs connecting them. The display includes
angles \\\alpha\\, \\\beta\\, and \\\gamma\\ relevant to solid angle
computation, providing a complete geometric view of the intersecting
cone configuration.

## References

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396 (math.MG).
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

Sievert, C. (2020). *Interactive Web-Based Data Visualization with R,
plotly, and Shiny*. Chapman and Hall/CRC. ISBN 978-1138331457.

## See also

[`solid_angle_intersecting_cones`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
for the closed-form measure of the region shown by the visualisation;
[`plot_cone_3d`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md)
for the simplicial-cone visualisation;
[`shinySolidAngleR`](https://robustecologies.github.io/SolidAngleR/reference/shinySolidAngleR.md)
for the dashboard that embeds these widgets.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default configuration with semi-transparent cones
plotIntersectingCones(color_palette = "viridis")

# Custom angles and palette
plotIntersectingCones(theta1 = pi/4, theta2 = pi/3, phi = pi/5,
                      color_palette = "magma", cone_opacity = 0.25)
} # }
```
