# Plot method for chameleon_test objects

Visualize the results of the chameleon attractor statistical test.

## Usage

``` r
# S3 method for class 'chameleon_test'
plot(x, type = c("null_dist", "metrics", "diagnostics", "all"), ...)
```

## Arguments

- x:

  Object of class "chameleon_test".

- type:

  Character. Type of plot:

  "null_dist"

  :   (Default) Histograms of null distributions with observed values
      marked as vertical lines.

  "metrics"

  :   Scale-dependent D(f) and theta(f) with confidence bands from
      surrogate distribution.

  "diagnostics"

  :   Diagnostic plots showing robustness check results.

  "all"

  :   Multi-panel summary combining all plot types.

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot2 or patchwork object.

## See also

[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
for the constructor;
[`print.chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_test.md),
[`summary.chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_test.md).
