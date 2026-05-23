# Plot method for qsd_estimate objects

Plots the marginal quasi-stationary distribution as histograms.

## Usage

``` r
# S3 method for class 'qsd_estimate'
plot(x, type = c("marginal", "survival"), ...)
```

## Arguments

- x:

  A `qsd_estimate` object

- type:

  Type of plot: "marginal" (default) or "survival"

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
