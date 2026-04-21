# Plot method for rare_event_estimate objects

Plots the distribution of log-weights or the convergence of the
probability estimate.

## Usage

``` r
# S3 method for class 'rare_event_estimate'
plot(x, type = c("weights", "convergence"), ...)
```

## Arguments

- x:

  A `rare_event_estimate` object

- type:

  Type of plot: "weights" (default) or "convergence"

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
