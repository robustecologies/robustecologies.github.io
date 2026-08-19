# Summarize a parameter sweep

Marginal success fractions per swept parameter value, the incumbent
cell, and the stable region extent.

## Usage

``` r
# S3 method for class 'sym_sweep'
summary(object, ...)

# S3 method for class 'summary.sym_sweep'
print(x, ...)
```

## Arguments

- object:

  A `sym_sweep` object.

- ...:

  Ignored.

- x:

  A `summary.sym_sweep` object.

## Value

An object of class `summary.sym_sweep`.

## See also

[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`print.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/print.sym_sweep.md),
[`plot.sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_sweep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sw <- sym_sweep(sym_objective(sym_benchmark("quadratic", d = 4)), "slc_expo",
                grid = list(C = 10^(-2:2), h = 10^(-2:1)), x0 = rep(2, 4))
summary(sw)
} # }
```
