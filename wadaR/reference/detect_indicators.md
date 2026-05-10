# Detect binary indicator parameters in a posterior

Scans the columns of a posterior matrix and identifies those whose
samples consist exclusively of values in {0, 1} (within floating point
tolerance). These columns are candidate spike-and-slab inclusion
indicators.

## Usage

``` r
detect_indicators(posterior, tol = 1e-08)
```

## Arguments

- posterior:

  A numeric matrix of posterior samples (rows = iterations, columns =
  parameters).

- tol:

  Numeric tolerance for testing whether values are 0 or 1. Default is
  `1e-8`.

## Value

Integer vector of column indices corresponding to binary indicator
parameters. Returns `integer(0)` if none found.

## Details

Implementation of `detect_indicators`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`plot.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/plot.ssvs_summary.md),
[`print.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/print.ssvs_summary.md),
[`summary.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/summary.ssvs_summary.md).

## Examples

``` r
if (FALSE) { # \dontrun{
post <- fit$Posterior2
binary_cols <- detect_indicators(post)
} # }
```
