# Plots of posterior predictive checks for SMC

Produces posterior predictive check plots for an object of class
`smc.ppc`, delegating to the shared internal PPC engine.

## Usage

``` r
# S3 method for class 'smc.ppc'
plot(x, Style = NULL, Data = NULL, Rows = NULL, PDF = FALSE, ...)
```

## Arguments

- x:

  An object of class `smc.ppc`.

- Style:

  Character string selecting the plot style. Defaults to `"Density"`.
  See
  [`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)
  for the full list of available styles.

- Data:

  Optional data list used when fitting the model. Required for certain
  plot styles.

- Rows:

  Optional integer vector of row numbers to plot. Defaults to all rows.

- PDF:

  Logical; save plots to a PDF file (default `FALSE`).

- ...:

  Additional arguments passed to the internal plotting engine, including
  `Group`, `stat_fun`, `loo`, and `forecast_start`.

## Value

Invisibly returns the plot object(s).

## Details

Implementation of `plot.smc.ppc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md),
[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## See the SMC function for an example.
} # }
```
