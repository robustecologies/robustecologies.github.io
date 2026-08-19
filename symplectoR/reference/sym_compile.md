# Compile a C++ objective for thread-safe evaluation

Compiles a user-supplied C++ objective body (and optionally a gradient
body) into function pointers with the fixed symplectoR signatures,
enabling the OpenMP-parallel evaluation routes of
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
multi-start ensembles and
[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md).
Compilation happens once per session per source; repeated calls with
identical code return the cached pointers.

## Usage

``` r
sym_compile(cpp_obj, cpp_grad = NULL, params = NULL)
```

## Arguments

- cpp_obj:

  Character string with the C++ body of the objective. The body sees
  `const arma::vec& x` (the point) and `const arma::vec& par` (the fixed
  parameter vector) and must return a `double`.

- cpp_grad:

  Optional character string with the C++ body of the gradient. The body
  sees the same arguments and must return an `arma::vec` of the same
  length as `x`.

- params:

  Optional numeric vector bound to `par` at evaluation time.

## Value

An object of class `sym_compiled_objective` carrying the external
pointers, the parameter vector and the source code. Pass it as the `f`
argument of
[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
which inherits the gradient and parameters.

## Details

The generated translation unit wraps the bodies as

    double sym_user_f(const arma::vec& x, const arma::vec& par) { <cpp_obj> }

and exports factory functions returning `Rcpp::XPtr` handles to them.
The compiled functions must not touch the R API; this is what makes them
safe to call inside OpenMP parallel regions, where R closures are
structurally excluded.

## References

Eddelbuettel, D., & Sanderson, C. (2014). RcppArmadillo: Accelerating R
with high-performance C++ linear algebra. *Computational Statistics &
Data Analysis*, 71, 1054-1063.
[doi:10.1016/j.csda.2013.02.005](https://doi.org/10.1016/j.csda.2013.02.005)

## See also

[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md),
[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
co <- sym_compile(
  cpp_obj = "return 0.5 * arma::dot(x, x);",
  cpp_grad = "return x;"
)
obj <- sym_objective(co, dim = 10, name = "half sum of squares")
fit <- sym_optim(obj, x0 = rep(1, 10), method = "rgd", n_starts = 8, n_threads = 4)
print(fit)
plot(fit, type = "dashboard")
} # }
```
