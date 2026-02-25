# Plot merging method results

S3 method to visualize the slim boundaries from the merging method using
ggplot2. Can display all boundaries overlaid or compare specific pairs.

## Usage

``` r
# S3 method for class 'wada_merging_result'
plot(x, basins = NULL, show_all = TRUE, boundary1 = 1, boundary2 = 2, ...)
```

## Arguments

- x:

  A `wada_merging_result` object from
  [`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md).

- basins:

  Optional original basins object for background (not currently used).

- show_all:

  Logical. If TRUE (default), displays all slim boundaries overlaid in
  different colors. If FALSE, compares two specific boundaries.

- boundary1:

  Integer. Index of first boundary for comparison when
  `show_all = FALSE`. Default is 1.

- boundary2:

  Integer. Index of second boundary for comparison when
  `show_all = FALSE`. Default is 2.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot2` object that can be further customized.

## Details

**Overlay mode (`show_all = TRUE`):**

Displays all \\N_A\\ slim boundaries in different colors. For Wada
basins, these boundaries should visually coincide (appear as a single
overlapping curve). Significant deviations indicate partial or non-Wada
basins.

**Comparison mode (`show_all = FALSE`):**

Compares two specific slim boundaries, showing:

- Both boundary point sets in different colors

- A dashed circle indicating the Hausdorff distance

- The point pair achieving the maximum distance

## See also

[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for the main analysis function,
[`hausdorff_distance`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md)
for the distance metric.

## Examples

``` r
if (FALSE) { # \dontrun{
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 300)
result <- wada_merging_method(basins)

# Show all slim boundaries overlaid (default)
plot(result)

# Compare boundaries 1 and 2
plot(result, show_all = FALSE, boundary1 = 1, boundary2 = 2)

# Compare boundaries 1 and 3
plot(result, show_all = FALSE, boundary1 = 1, boundary2 = 3)

# Customize with ggplot2
library(ggplot2)
plot(result) +
    theme_dark() +
    labs(title = "Slim boundaries (merging method)")
} # }
```
