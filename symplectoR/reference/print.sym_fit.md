# Print a symplectic optimization fit

One-screen report of a
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
run: method, incumbent value, convergence state and evaluation budget.

## Usage

``` r
# S3 method for class 'sym_fit'
print(x, ...)
```

## Arguments

- x:

  A `sym_fit` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`summary.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_fit.md),
[`plot.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- sym_optim(sym_objective(sym_benchmark("quadratic", d = 5)),
                 x0 = rep(0, 5), method = "rgd")
print(fit)
} # }
```
