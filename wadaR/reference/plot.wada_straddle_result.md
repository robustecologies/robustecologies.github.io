# Plot saddle-straddle method results

S3 method to visualize the computed chaotic saddles from the
saddle-straddle method using ggplot2. Can display all saddles overlaid
or in separate facets.

## Usage

``` r
# S3 method for class 'wada_straddle_result'
plot(x, basins = NULL, overlay = TRUE, ...)
```

## Arguments

- x:

  A `wada_straddle_result` object from
  [`wada_straddle_method`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md).

- basins:

  Optional `wada_basins` object to display as background. When provided,
  saddles are overlaid on a semi-transparent basin plot.

- overlay:

  Logical. If TRUE (default), displays all saddles overlaid in different
  colors on a single plot. If FALSE, shows each saddle in a separate
  faceted panel.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot2` object that can be further customized.

## Details

**Overlay mode (`overlay = TRUE`):**

Displays all \\N_A\\ computed saddle sets in different colors on a
single plot. For Wada basins, all saddles should visually overlap (same
chaotic set). Deviations indicate multiple distinct saddles (partial or
non-Wada).

When a `basins` object is provided, the saddles are overlaid on a
semi-transparent basin plot, showing the relationship between the
chaotic saddle and basin structure.

**Faceted mode (`overlay = FALSE`):**

Displays each saddle in a separate panel using ggplot2 facets. Useful
for comparing individual saddle geometries and identifying differences.

## See also

[`wada_straddle_method`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
for the main analysis function,
[`plot.wada_basins`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md)
for plotting basins.

## Examples

``` r
if (FALSE) { # \dontrun{
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)
result <- wada_straddle_method(pendulum,
                               x_range = c(-pi, pi),
                               y_range = c(-3, 3),
                               n_points = 5000)

# Overlay all saddles (default)
plot(result)

# Saddles overlaid on basins
plot(result, basins = basins)

# Individual saddles in faceted plot
plot(result, overlay = FALSE)

# Customize with ggplot2
library(ggplot2)
plot(result, basins = basins) +
    theme_dark() +
    labs(title = "Chaotic saddle structure")
} # }
```
