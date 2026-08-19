# Adapt a symplectoR method to an optim-style local solver

Returns a function with the signature `function(par, fn, lower, upper)`
returning `list(par, value, convergence)`, the drop-in shape expected by
multi-start global optimizers such as the local-descent slot of
[`RElabverse::stiff_globopt()`](https://rdrr.io/pkg/RElabverse/man/stiff_globopt.html).
This lets any symplectoR method serve as the incumbent generator inside
an existing box-constrained search loop.

## Usage

``` r
as_incumbent_solver(method = "rgd", control = sym_control(), ...)
```

## Arguments

- method:

  A trajectory method accepted by
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
  (`"qhd"` is excluded: it is itself a global method and needs the full
  objective specification rather than a bare function).

- control:

  A `sym_control` for the chosen method.

- ...:

  Further arguments forwarded to
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
  (for example `n_starts`, `seed`).

## Value

A function `(par, fn, lower, upper)` returning a list with elements
`par` (the incumbent), `value` (its objective value) and `convergence`
(`0L` on convergence, `1L` otherwise), matching the
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) convention.

## Details

When both bounds are finite the adapter evaluates the objective exactly
and folds iterates into the box by reflection; otherwise the problem is
treated as unconstrained. The wrapped function is called through the
R-closure route, so the adapter is serial by construction; parallelism
belongs to the outer search loop that calls it.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`sym_control()`](https://robustecologies.github.io/symplectoR/reference/sym_control.md),
[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md)

## Examples

``` r
if (FALSE) { # \dontrun{
solver <- as_incumbent_solver("rgd", sym_control("rgd", eps = 0.05, max_iter = 2000))
res <- solver(par = c(0, 0), fn = function(x) sum((x - 1)^2),
              lower = c(-5, -5), upper = c(5, 5))
str(res)
} # }
```
