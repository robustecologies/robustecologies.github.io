# Plot method for merged_basins

Creates a visualization of merged basins showing the kept basin versus
the merged basin, with optional slim boundary overlay. The slim boundary
is the boundary between the two merged regions, which is identical for
all merging configurations in true Wada basins.

## Usage

``` r
# S3 method for class 'merged_basins'
plot(
  x,
  colors = NULL,
  show_boundary = FALSE,
  boundary_color = "black",
  boundary_size = 0.3,
  title = NULL,
  plotly = FALSE,
  caption = TRUE,
  ...
)
```

## Arguments

- x:

  A `merged_basins` object from
  [`merge_basins`](https://robustecologies.github.io/wadaR/reference/merge_basins.md).

- colors:

  Character vector of length 2. Colors for the kept basin and merged
  basin respectively. Default uses complementary colors.

- show_boundary:

  Logical. If TRUE, overlays the slim boundary in black. Default is
  FALSE.

- boundary_color:

  Color for the slim boundary overlay. Default is "black".

- boundary_size:

  Numeric. Size of boundary points. Default is 0.3.

- title:

  Character. Plot title. If NULL, generates an informative title.

- plotly:

  Logical. If TRUE, returns an interactive plotly plot instead of a
  static ggplot2 plot. Default is FALSE.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name and primary citation. Set to FALSE to suppress the
  caption when the figure is composed with other panels.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot2` object (if `plotly = FALSE`) or a `plotly` object (if
`plotly = TRUE`) that can be further customized.

## Details

The plot shows the two-color merged basin map where:

- The **kept basin** \\B_k\\ is shown in the first color

- The **merged basin** \\M_k = \bigcup\_{j \neq k} B_j\\ is shown in the
  second color

The **slim boundary** between these two regions is the key element for
the merging method. For Wada basins, this boundary is identical
regardless of which basin is kept separate. This can be verified by
comparing plots with `show_boundary = TRUE` for different values of
`keep_basin`.

## References

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining when
a basin is Wada: The merging method. *Scientific Reports*, 8, 9954.
[doi:10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

## See also

[`merge_basins`](https://robustecologies.github.io/wadaR/reference/merge_basins.md)
for creating merged basin objects,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for the full Wada detection algorithm.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Basic merged basin visualization                            #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 200)

merged <- merge_basins(basins, keep_basin = 1)
plot(merged)

# ===================================================================== #
# Example 2: Show slim boundary                                          #
# ===================================================================== #
plot(merged, show_boundary = TRUE, boundary_color = "white")

# ===================================================================== #
# Example 3: Compare different merging configurations                    #
# ===================================================================== #
# For Wada basins, slim boundaries should be identical
library(patchwork)

p1 <- plot(merge_basins(basins, 1), show_boundary = TRUE,
           title = "B1 vs (B2+B3)")
p2 <- plot(merge_basins(basins, 2), show_boundary = TRUE,
           title = "B2 vs (B1+B3)")
p3 <- plot(merge_basins(basins, 3), show_boundary = TRUE,
           title = "B3 vs (B1+B2)")

p1 + p2 + p3

# ===================================================================== #
# Example 4: Custom colors                                               #
# ===================================================================== #
plot(merged, colors = c("#264653", "#E9C46A"),
     title = "Merged basins with custom palette")
} # }
```
