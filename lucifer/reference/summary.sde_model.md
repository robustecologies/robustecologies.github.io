# Summary method for SDE model objects

Prints an extended summary of the SDE model specification, including
data statistics, method details, and initial parameter values.

## Usage

``` r
# S3 method for class 'sde_model'
summary(object, ...)
```

## Arguments

- object:

  An object of class `sde_model`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the summary information as a list.

## Details

Implementation of `summary.sde_model`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
[`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md),
[`plot.sde_model`](https://robustecologies.github.io/lucifer/reference/plot.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`simulate.sde_model`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving summary.sde_model
} # }
```
