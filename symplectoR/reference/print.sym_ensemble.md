# Print a multi-start ensemble

One-screen report of a multi-start
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
run: incumbent value, start statistics and thread usage.

## Usage

``` r
# S3 method for class 'sym_ensemble'
print(x, ...)
```

## Arguments

- x:

  A `sym_ensemble` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`summary.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_ensemble.md),
[`plot.sym_ensemble()`](https://robustecologies.github.io/symplectoR/reference/plot.sym_ensemble.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- sym_compile("return arma::dot(x, x);", "return 2 * x;")
obj <- sym_objective(co, dim = 5, lower = -3, upper = 3)
ens <- sym_optim(obj, method = "rgd", n_starts = 16, n_threads = 4, seed = 1)
print(ens)
} # }
```
