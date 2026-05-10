# Plot 3D attractor from Takens embedding

Visualize the reconstructed attractor in 3D phase space using
interactive plotly visualization, optionally colored by instantaneous
dimension or inverse persistence. This is the S3 plot method for objects
of class `takens_embedding` returned by
[`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md).

## Usage

``` r
# S3 method for class 'takens_embedding'
plot(
  x,
  color_by = NULL,
  color_label = "",
  palette = "viridis",
  point_size = 2,
  alpha = 0.7,
  axes = 1:3,
  main = NULL,
  ...
)
```

## Arguments

- x:

  Object of class `takens_embedding` (a matrix with at least 3 columns
  representing the embedded trajectory).

- color_by:

  Optional numeric vector. Values to use for coloring points. Length
  must match `nrow(x)`. Common choices are `d` or `theta` from
  [`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md).

- color_label:

  Character. Label for the colorbar (e.g., "D(t)" or "theta(t)").

- palette:

  Character. Color palette name: "viridis" (default), "plasma",
  "inferno", "magma", "cividis", "turbo", or a vector of colors.

- point_size:

  Numeric. Point size for the scatter plot (default 2).

- alpha:

  Numeric. Point transparency, between 0 (fully transparent) and 1
  (fully opaque). Default is 0.7.

- axes:

  Integer vector of length 3. Which columns of the embedding matrix to
  use as x, y, z coordinates. Default is `1:3`.

- main:

  Character. Plot title. If NULL, a default title is generated.

- ...:

  Additional arguments passed to plotly.

## Value

A plotly htmlwidget object (interactive 3D plot). If plotly is not
available, returns NULL invisibly after creating a base R plot.

## Details

The function creates an interactive 3D scatter plot using plotly, where
each point represents a state in the reconstructed phase space. Points
can be colored by any metric vector of matching length (typically the
instantaneous dimension d or inverse persistence theta from EVT
analysis).

If plotly is not available, the function falls back to a 2D projection
using base R graphics.

The Takens embedding reconstructs the attractor from a scalar time
series. Visualizing the 3D projection reveals the geometric structure of
the underlying dynamical system. Coloring by EVT metrics (dimension d or
persistence theta) highlights regions with different dynamical
properties, which is essential for detecting chameleon behavior.

The axis labels use Greek theta notation to indicate delay coordinates.

## See also

[`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
for creating embeddings,
[`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
for computing coloring values,
[`plot.chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_analysis.md)
for complete analysis visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Simple quasi-periodic attractor (torus-like)
set.seed(42)
t <- seq(0, 100, by = 0.02)
x <- sin(t) + 0.5*sin(sqrt(2)*t) + 0.3*sin(sqrt(3)*t) + 0.05*rnorm(length(t))
embedded <- takens_embed(x, m = 3, tau = 15)

# Inspect and visualize
print(embedded)
summary(embedded)
plot(embedded)

# Example 2: Bistable system with intermittent switching (two-lobe attractor)
set.seed(123)
n <- 10000
x <- numeric(n)
x[1] <- 0.1
for (i in 2:n) {
  x[i] <- x[i-1] + 0.1*(x[i-1] - x[i-1]^3) + 0.15*rnorm(1)
}
embedded <- takens_embed(x, m = 3, tau = 5)
plot(embedded, main = "Two-lobe attractor")

# Example 3: With EVT metrics coloring
metrics <- evt_metrics(embedded)
plot(embedded, color_by = metrics$d, color_label = "D(t)",
     palette = "Jet", point_size = 3)
plot(embedded, color_by = metrics$theta, color_label = "\u03B8(t)",
     palette = "Hot")
} # }
```
