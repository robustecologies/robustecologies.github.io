# Summarize a symplectic optimization fit

Richer diagnostics than
[`print.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/print.sym_fit.md):
parameter estimates, trace tail behaviour, restart and trust-region
statistics where the method provides them.

## Usage

``` r
# S3 method for class 'sym_fit'
summary(object, ...)

# S3 method for class 'summary.sym_fit'
print(x, ...)
```

## Arguments

- object:

  A `sym_fit` object.

- ...:

  Ignored.

- x:

  A `summary.sym_fit` object.

## Value

An object of class `summary.sym_fit`.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`print.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/print.sym_fit.md),
[`plot.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- sym_optim(sym_objective(sym_benchmark("quadratic", d = 5)),
                 x0 = rep(0, 5), method = "rgd")
summary(fit)
} # }
```
