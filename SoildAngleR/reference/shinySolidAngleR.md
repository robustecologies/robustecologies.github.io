# Launch interactive SolidAngleR Shiny application

Launches a comprehensive Shiny dashboard application demonstrating all
major functions of the SolidAngleR package for computing solid angles of
polyhedral cones in arbitrary dimensions.

## Usage

``` r
shinySolidAngleR(launch.browser = TRUE)
```

## Arguments

- launch.browser:

  Logical. If TRUE (default), opens the app in a browser.

## Value

Launches a Shiny application (invisibly returns NULL).

## Details

The application provides interactive exploration of solid angle
computation across all dimensions, including 2D and 3D visualization
modules, circular and intersecting cone explorers, high-dimensional
methods such as hypergeometric series and tridiagonal optimization, the
decomposition method for general cones, Monte Carlo estimation with
uncertainty quantification, ecological feasibility domain analysis via
the omega function, uniform cone sampling visualization, and
comprehensive theory vignettes with mathematical background.

## Examples

``` r
if (FALSE) { # \dontrun{
shinySolidAngleR()
} # }
```
