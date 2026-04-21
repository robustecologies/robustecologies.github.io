# Discrete Lyapunov function for a map

Builds a local quadratic Lyapunov function for a discrete-time dynamical
system specified via
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
with the `map` argument. Given a fixed point \\x^\*\\ and the
linearisation \\A = DF(x^\*)\\, solves the discrete Lyapunov equation
\\A^\top P A - P = -Q\\ and returns \\V(x) = (x - x^\*)^\top P (x -
x^\*)\\, which certifies local asymptotic stability when the spectral
radius \\\rho(A) \< 1\\.

## Usage

``` r
lyapunov_discrete(
  model,
  x_star = NULL,
  Q = NULL,
  tol = sqrt(.Machine$double.eps),
  verbose = TRUE,
  ...
)
```

## Arguments

- model:

  A `model_spec` object with `is_map = TRUE`.

- x_star:

  Numeric vector, the fixed point. If `NULL`, tries the model's default
  initial condition as Newton starting point.

- Q:

  Symmetric positive-definite weight matrix. Defaults to the identity.

- tol:

  Tolerance for the Hurwitz (discrete) test; the spectral radius must be
  below `1 - tol`. Defaults to `sqrt(.Machine$double.eps)`.

- verbose:

  Logical. When `TRUE`, narrates the construction using coloured
  symbols. Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function` with the same structure as the
continuous-time constructors; the `type` field is `"discrete"`. The
`certificate$method_details` field includes the spectral radius and the
residual of the discrete Lyapunov equation.

## Details

The discrete Lyapunov equation \\A^\top P A - P = -Q\\ has a unique
symmetric positive-definite solution \\P\\ if and only if the
eigenvalues of \\A\\ lie strictly inside the open unit disk (Khalil,
2002, Theorem 4.11). The solver uses the Bartels-Stewart variant
`bartels_stewart_discrete_rcpp` already in the package.

For a nonlinear map the resulting \\V\\ is a local certificate
(Lyapunov's indirect method for discrete systems). The certificate label
is therefore `"local_algebraic"` when `is_map` and the map is nonlinear,
and `"algebraic"` when the map is linear.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

Antsaklis, P. J., & Michel, A. N. (2006). *Linear Systems*. Birkhauser.
[doi:10.1007/0-8176-4435-0](https://doi.org/10.1007/0-8176-4435-0)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
logistic <- model_spec(
    map = list(x ~ r * x * (1 - x)),
    state_names = "x",
    parms = list(r = 2.8),
    init = c(x = 0.5)
)
lf <- lyapunov_discrete(logistic, x_star = 1 - 1/2.8)
print(lf); summary(lf)
} # }
```
