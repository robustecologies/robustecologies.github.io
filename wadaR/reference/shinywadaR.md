# Launch interactive wadaR Shiny application

Launches a comprehensive Shiny dashboard application demonstrating all
major functions of the wadaR package for Wada basin detection and
analysis.

## Usage

``` r
shinywadaR(launch.browser = TRUE)
```

## Arguments

- launch.browser:

  Logical. If TRUE (default), opens the app in a browser.

## Value

Launches a Shiny application (invisibly returns NULL).

## Details

The application includes:

- Basin computation for multiple dynamical systems (pendulum,
  Henon-Heiles, Newton fractals, multispecies competition)

- Wada detection methods (grid, merging, saddle-straddle)

- Basin entropy calculation and visualization

- Merged basin analysis

- Theory vignettes with mathematical background

## Examples

``` r
if (FALSE) { # \dontrun{
shinywadaR()
} # }
```
