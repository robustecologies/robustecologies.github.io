# Plot method for compiled system objects

Displays a schematic of the dynamical system showing attractor locations
in a 2D projection of state space.

## Usage

``` r
# S3 method for class 'compiled_system'
plot(x, caption = TRUE, ...)
```

## Arguments

- x:

  A compiled_system object.

- caption:

  Logical. If TRUE (default), render a one-line caption with the
  function name. Set to FALSE to suppress the caption when the figure is
  composed with other panels.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object, invisibly.
