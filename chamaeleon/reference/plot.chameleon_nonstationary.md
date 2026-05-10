# Plot method for chameleon_nonstationary objects

Visualize time-varying chameleon behavior.

## Usage

``` r
# S3 method for class 'chameleon_nonstationary'
plot(x, type = c("time_series", "heatmap", "classification", "all"), ...)
```

## Arguments

- x:

  Object of class "chameleon_nonstationary".

- type:

  Character. Type of plot:

  "time_series"

  :   (Default) Chameleon index over time with significance threshold
      and change points.

  "heatmap"

  :   Heat map of scale-dependent metrics over time.

  "classification"

  :   Binary classification of windows as chameleon/non-chameleon over
      time.

  "all"

  :   Multi-panel summary.

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot2 or patchwork object.

## See also

[`chameleon_test_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test_nonstationary.md)
for the constructor;
[`print.chameleon_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_nonstationary.md),
[`summary.chameleon_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_nonstationary.md).
