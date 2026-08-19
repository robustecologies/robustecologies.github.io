# Summarize a multi-start ensemble

Per-start table statistics and the incumbent fit summary.

## Usage

``` r
# S3 method for class 'sym_ensemble'
summary(object, ...)

# S3 method for class 'summary.sym_ensemble'
print(x, ...)
```

## Arguments

- object:

  A `sym_ensemble` object.

- ...:

  Ignored.

- x:

  A `summary.sym_ensemble` object.

## Value

An object of class `summary.sym_ensemble`.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`print.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/print.sym_ensemble.md),
[`plot.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_ensemble.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ens <- sym_optim(sym_objective(sym_benchmark("rosenbrock", d = 4)),
                 x0 = rep(0, 4), method = "rgd", n_starts = 8, seed = 1)
summary(ens)
} # }
```
