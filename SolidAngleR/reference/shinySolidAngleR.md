# Launch the SolidAngleR Shiny dashboard

Starts an interactive shinydashboard application that exposes the main
functionality of the package: low-dimensional formulas, series methods,
decomposition, Monte Carlo, dispatcher diagnostics, and 3D
visualisation. Intended as an exploratory tool and a teaching resource.

## Usage

``` r
shinySolidAngleR(launch.browser = TRUE)
```

## Arguments

- launch.browser:

  Logical. If `TRUE` (default), opens the running app in the system
  default browser; if `FALSE`, the server runs in the current R session
  and the URL is printed.

## Value

`NULL`, invisibly. The function is called for the side effect of
starting the Shiny server.

## Details

The dashboard is structured into thematic tabs covering 2D and 3D
visualisation modules, circular and intersecting cone explorers,
high-dimensional methods (hypergeometric and tridiagonal series), the
decomposition route for general cones, Monte Carlo estimation with
uncertainty quantification, ecological feasibility domains via the
multivariate-normal branch of
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md),
uniform cone-sampling diagnostics, and a theory tab with the
mathematical background.

Required Suggested packages: `shiny`, `shinydashboard`, `shinyWidgets`,
`ggplot2`, `DT`, `plotly`. The function checks for them at start-up and
stops with an installation hint if any are missing.

## References

Chang, W., Cheng, J., Allaire, J. J., Sievert, C., Schloerke, B., Xie,
Y., Allen, J., McPherson, J., Dipert, A., & Borges, B. (2024). *shiny:
Web Application Framework for R*. R package version 1.8 or later.
<https://shiny.posit.co/>.

## See also

[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md),
[`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)
for the analytical engines exposed in the dashboard;
[`plot_cone_3d`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md),
[`plotIntersectingCones`](https://robustecologies.github.io/SolidAngleR/reference/plotIntersectingCones.md)
for the embedded 3D visualisations.

## Examples

``` r
if (FALSE) { # \dontrun{
shinySolidAngleR()
} # }
```
