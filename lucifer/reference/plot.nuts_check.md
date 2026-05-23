# Plot NUTS diagnostics

Produces a publication-quality 4-panel diagnostic plot for NUTS
sampling. When multiple chains are present, each chain is rendered in a
distinct color so convergence and mixing can be assessed visually.

## Usage

``` r
# S3 method for class 'nuts_check'
plot(x, ...)
```

## Arguments

- x:

  A `nuts_check` object from
  [`check_nuts`](https://robustecologies.github.io/lucifer/reference/check_nuts.md).

- ...:

  Ignored.

## Value

Invisible `NULL`.

## Details

Implementation of `plot.nuts_check`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`print.nuts_check`](https://robustecologies.github.io/lucifer/reference/print.nuts_check.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.nuts_check
} # }
```
