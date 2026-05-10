# Plot method for mlmc_estimate objects

Visualizes the MLMC estimation diagnostics. The default "convergence"
type shows log-variance and log-cost as functions of level, which should
be approximately linear. The "allocation" type shows the optimal vs.
actual number of samples per level.

## Usage

``` r
# S3 method for class 'mlmc_estimate'
plot(x, type = c("convergence", "allocation"), ...)
```

## Arguments

- x:

  A `mlmc_estimate` object

- type:

  Type of plot: "convergence" (default) or "allocation"

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
