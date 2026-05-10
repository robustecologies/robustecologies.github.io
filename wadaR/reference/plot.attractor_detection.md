# Plot method for attractor_detection objects

Plot method for attractor_detection objects

## Usage

``` r
# S3 method for class 'attractor_detection'
plot(
  x,
  show_endpoints = TRUE,
  show_circles = TRUE,
  point_alpha = 0.5,
  caption = TRUE,
  ...
)
```

## Arguments

- x:

  An attractor_detection object.

- show_endpoints:

  Logical. Show trajectory endpoints.

- show_circles:

  Logical. Show attractor convergence radii.

- point_alpha:

  Numeric. Transparency for endpoints.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name. Set to FALSE to suppress the caption when the figure is
  composed with other panels.

- ...:

  Additional arguments passed to ggplot.

## Value

A ggplot object.
