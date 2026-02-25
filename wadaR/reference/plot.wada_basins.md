# Plot method for wada_basins (S3)

Creates publication-quality visualizations of basins of attraction using
ggplot2. Supports custom color palettes and boundary overlay.

## Usage

``` r
# S3 method for class 'wada_basins'
plot(
  x,
  colors = NULL,
  show_boundary = FALSE,
  boundary_color = "black",
  unclassified_color = "gray80",
  title = "Basins of attraction",
  ...
)
```

## Arguments

- x:

  wada_basins object from
  [`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md).

- colors:

  Character vector of colors for each basin. If NULL, uses the HCL
  "Set2" palette. Length should match number of attractors.

- show_boundary:

  Logical. If TRUE, overlays boundary points as dots.

- boundary_color:

  Color for boundary point overlay.

- unclassified_color:

  Color for unclassified points (basin = 0).

- title:

  Plot title.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object that can be printed or further modified.

## Details

The plot uses
[`geom_raster()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
for efficient rendering of the basin grid. Basin colors are assigned
using a discrete scale, with unclassified points (those that did not
converge to any attractor) shown in gray by default.

For Wada basins, the fractal boundary structure can be highlighted by
setting `show_boundary = TRUE`, which overlays boundary points detected
by
[`get_boundary`](https://robustecologies.github.io/wadaR/reference/get_boundary.md).

The returned ggplot2 object can be further customized using standard
ggplot2 functions (e.g., adding annotations, changing themes).

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins,
[`get_boundary`](https://robustecologies.github.io/wadaR/reference/get_boundary.md)
for boundary extraction,
[`plot.wada_grid_result`](https://robustecologies.github.io/wadaR/reference/plot.wada_grid_result.md)
for grid method results.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Basic basin plot                                                      #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)
plot(basins)

# ===================================================================== #
# Custom colors with boundary overlay                                   #
# ===================================================================== #
plot(basins,
     colors = c("#E41A1C", "#377EB8", "#4DAF4A"),
     show_boundary = TRUE,
     boundary_color = "white",
     title = "Wada basins with fractal boundary")

# ===================================================================== #
# Add custom annotations to the ggplot object                           #
# ===================================================================== #
library(ggplot2)
p <- plot(basins, title = "Forced damped pendulum")
p + annotate("text", x = 0, y = 2.5, label = "F = 1.66",
             color = "white", size = 5) +
    theme(legend.position = "bottom")

# ===================================================================== #
# Compare different parameter values                                    #
# ===================================================================== #
library(patchwork)  # For combining plots

p1 <- plot(compute_basins(forced_damped_pendulum(forcing = 1.66),
           c(-pi, pi), c(-3, 3), resolution = 200),
           title = "F = 1.66 (Wada)")
p2 <- plot(compute_basins(forced_damped_pendulum(forcing = 1.71),
           c(-pi, pi), c(-3, 3), resolution = 200),
           title = "F = 1.71 (Partial Wada)")
p1 + p2
} # }
```
