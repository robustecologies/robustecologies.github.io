# Get basin boundary points

Identifies points on the boundary between basins using a parallel Rcpp
implementation, returning coordinates and neighbor-basin counts for each
boundary cell.

## Usage

``` r
get_boundary(basins, x_grid = NULL, y_grid = NULL)
```

## Arguments

- basins:

  Matrix of basin assignments or wada_basins object.

- x_grid:

  Optional vector of x coordinates.

- y_grid:

  Optional vector of y coordinates.

## Value

Data frame with columns:

- x:

  X coordinate of boundary point

- y:

  Y coordinate of boundary point

- basin:

  Basin assignment of the point

- n_neighbors:

  Number of distinct basins in the neighborhood

## Details

The boundary detection uses 4-connectivity (von Neumann neighborhood): a
point is on the boundary if any of its four cardinal neighbors belongs
to a different basin. The number of distinct neighboring basins is also
computed, which is useful for identifying Wada points (points adjacent
to three or more basins).

For Wada basins, all boundary points should eventually touch all basins
at sufficiently fine resolution, since the boundary is shared by all
basins simultaneously.

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
[`merge_basins`](https://robustecologies.github.io/wadaR/reference/merge_basins.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 200)

# Extract boundary points
boundary <- get_boundary(basins)
head(boundary)

# Count points touching 2, 3, or more basins
table(boundary$n_neighbors)

# Visualize boundary colored by number of adjacent basins
library(ggplot2)
ggplot(boundary, aes(x = x, y = y, color = factor(n_neighbors))) +
    geom_point(size = 0.3) +
    scale_color_viridis_d(name = "Adjacent\nbasins") +
    coord_fixed() +
    theme_minimal()
} # }
```
