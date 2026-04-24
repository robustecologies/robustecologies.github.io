# Interactive 3D plotly visualisation of a simplicial cone

Builds an interactive plotly widget showing a simplicial cone in
\\\mathbb{R}^3\\: the three generator rays, the planar faces, and the
spherical triangle traced by the cone on the unit sphere. The widget is
suitable for embedding in pkgdown articles and Shiny apps.

## Usage

``` r
plot_cone_3d(
  V,
  ray_colors = NULL,
  surface_color = "orange",
  opacity = 0.5,
  show_sphere = TRUE,
  resolution = 60,
  show_legend = TRUE
)
```

## Arguments

- V:

  A 3x3 matrix where columns are the generating vectors of the cone.

- ray_colors:

  Vector of 3 colors for the generator rays. If NULL, defaults to "red",
  "green", "blue".

- surface_color:

  Color for the cone surface (default: "orange").

- opacity:

  Opacity of the cone surface (default: 0.3).

- show_sphere:

  Logical; whether to show a reference unit sphere (default: TRUE).

- resolution:

  Integer; resolution for the arcs and surface mesh (default: 60).

- show_legend:

  Logical; whether to show the legend (default: TRUE).

## Value

A plotly object containing the visualization.

## Details

The visualization includes generator rays (unit vectors extending from
the origin, each clickable to trigger JavaScript console logs and Shiny
events), spherical arcs (geodesic arcs connecting the tips of the
vectors, perfectly aligned with the surface mesh), and the spherical
triangle (the region on the unit sphere enclosed by the cone,
constructed using barycentric subdivision and spherical projection).
Each ray displays tooltips with vector coordinates.

The plot title includes the calculated normalized solid angle measure.

## References

Sievert, C. (2020). *Interactive Web-Based Data Visualization with R,
plotly, and Shiny*. Chapman and Hall/CRC. ISBN 978-1138331457.

Mazonka, O. (2012). Solid angle of conical surfaces, polyhedral cones,
and intersecting spherical caps. arXiv:1205.1396 (math.MG).
[doi:10.48550/arXiv.1205.1396](https://doi.org/10.48550/arXiv.1205.1396)

## See also

[`plotIntersectingCones`](https://robustecologies.github.io/SolidAngleR/reference/plotIntersectingCones.md)
for the visualisation of two intersecting circular cones;
[`plot.cone_diagnosis`](https://robustecologies.github.io/SolidAngleR/reference/plot.cone_diagnosis.md)
for the diagnostic eigenvalue plot;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the scalar measure shown in the title;
[`shinySolidAngleR`](https://robustecologies.github.io/SolidAngleR/reference/shinySolidAngleR.md)
for the dashboard that embeds these plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Octant in R^3
plot_cone_3d(diag(3), ray_colors = c("orange", "purple", "blue"))

# Narrow cone with custom rays
V2 <- matrix(c(1, 0,   0,
               1, 0.2, 0,
               1, 0,   0.2), nrow = 3, byrow = FALSE)
plot_cone_3d(V2, ray_colors = c("green", "red", "blue"), opacity = 0.5)
} # }
```
