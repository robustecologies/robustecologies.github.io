# Print method for coupled MCMC objects

Displays a concise summary of a coupled MCMC run including meeting
status, meeting time, chain dimensions, and wall-clock time.

## Usage

``` r
# S3 method for class 'coupled_mcmc'
print(x, ...)
```

## Arguments

- x:

  An object of class `coupled_mcmc`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `x`.

## Details

Implementation of `print.coupled_mcmc`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md),
[`plot.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/plot.coupled_mcmc.md),
[`summary.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/summary.coupled_mcmc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.coupled_mcmc
} # }
```
