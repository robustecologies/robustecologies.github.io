# Plot grid method results

S3 method to visualize the boundary classification from the grid method
using ggplot2. Colors indicate the number of distinct basins bordering
each grid cell.

## Usage

``` r
# S3 method for class 'wada_grid_result'
plot(
  x,
  basins = NULL,
  show_wada_only = FALSE,
  plotly = FALSE,
  caption = TRUE,
  ...
)
```

## Arguments

- x:

  A `wada_grid_result` object from
  [`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md).

- basins:

  Original basins object (`wada_basins`) for coordinate mapping. If
  NULL, uses pixel indices.

- show_wada_only:

  Logical. If TRUE, only displays Wada boundary points (cells bordering
  all \\N_A\\ basins) in black. Default is FALSE.

- plotly:

  Logical. If TRUE, returns an interactive plotly plot instead of a
  static ggplot2 plot. Default is FALSE.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name and primary citation. Set to FALSE to suppress the
  caption when the figure is composed with other panels.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot2` object (if `plotly = FALSE`) or a `plotly` object (if
`plotly = TRUE`) that can be further customized.

## Details

The plot shows the boundary classification matrix where each cell is
colored according to the number of distinct basins in its neighborhood:

- **White**: Interior points (only one basin nearby)

- **Light colors**: Border between 2 basins

- **Dark colors**: Border between 3+ basins

- **Darkest**: Wada boundary (borders all \\N_A\\ basins)

M-values from 2 to \\N_A\\ count the number of distinct basins reachable
from each refined boundary box under Sweet-Ott-Yorke refinement.

When `show_wada_only = TRUE`, only the Wada boundary points are shown in
black against a white background, highlighting the fractal structure of
the Wada boundary.

## See also

[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
for the main analysis function,
[`plot.wada_basins`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md)
for plotting the basins themselves.

## Examples

``` r
if (FALSE) { # \dontrun{
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)
result <- wada_grid_method(basins)

# Standard boundary classification plot
plot(result, basins = basins)

# Only Wada boundary points
plot(result, basins = basins, show_wada_only = TRUE)

# Customize with ggplot2
library(ggplot2)
plot(result, basins = basins) +
    theme_dark() +
    labs(title = "Wada boundary structure")
} # }
```
