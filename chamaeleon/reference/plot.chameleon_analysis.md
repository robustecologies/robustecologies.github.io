# Plot chameleon analysis summary

Create a multi-panel summary plot for chameleon analysis results. This
S3 method provides a comprehensive visualization of the complete
analysis, including the time series, attractor projections, and
scale-dependent metrics.

## Usage

``` r
# S3 method for class 'chameleon_analysis'
plot(x, interactive = FALSE, ...)
```

## Arguments

- x:

  Object of class `chameleon_analysis`.

- interactive:

  Logical. Use plotly for interactive 3D attractor visualization.
  Default is FALSE.

- ...:

  Additional arguments passed to plotting functions.

## Value

If `interactive=TRUE`, returns a plotly htmlwidget for the 3D attractor.
Otherwise returns a patchwork ggplot2 composition invisibly.

## Details

The default base graphics output shows five panels: the original (and
filtered) time series, two 2D projections of the attractor (XY and XZ),
and the scale-dependent dimension and persistence.

For interactive 3D visualization of the attractor, use
`interactive=TRUE` which returns a plotly widget with the attractor
colored by dimension.

The ggplot2 summary provides a quick overview of the analysis across
five panels:

- **Panel 1**: Time series showing original and filtered data

- **Panels 2-3**: 2D projections of the reconstructed attractor

- **Panels 4-5**: Scale-dependent dimension and persistence

The interactive mode creates a 3D scatter plot where points are colored
by their instantaneous dimension, allowing exploration of which regions
of the attractor show different dynamical properties.

## See also

[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for running the analysis,
[`plot.takens_embedding`](https://robustecologies.github.io/chamaeleon/reference/plot.takens_embedding.md)
for standalone attractor plots,
[`plot.scale_metrics`](https://robustecologies.github.io/chamaeleon/reference/plot.scale_metrics.md)
for detailed scale metrics visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- chameleon_analysis(x, verbose = TRUE)
plot(result)
plot(result, interactive = TRUE)
} # }
```
