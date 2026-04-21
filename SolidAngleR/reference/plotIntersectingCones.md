# Interactive 3D visualization of intersecting cones

Creates an interactive plotly visualization of two intersecting cones on
a unit sphere, showing their intersection points, circles, arcs, and
angles. This visualization helps understand the geometry of intersecting
cones as described in Mazonka (2012).

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
<https://arxiv.org/abs/1205.1396>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Default visualization with semi-transparent cones
fig1 <- plotIntersectingCones(color_palette = "viridis")
fig1

# Example 2: Show cone rays with very transparent surfaces
fig2 <- plotIntersectingCones(show_rays = TRUE, cone_opacity = 0.1,
                              color_palette = "plasma")
fig2

# Example 3: More opaque cones without rays
fig3 <- plotIntersectingCones(cone_opacity = 0.5, color_palette = "inferno")
fig3

# Example 4: Custom angles and color palette
fig4 <- plotIntersectingCones(theta1 = pi/4, theta2 = pi/3, phi = pi/5,
                              color_palette = "magma", cone_opacity = 0.25)
fig4
} # }
```
