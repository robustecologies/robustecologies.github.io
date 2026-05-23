# Lyapunov-Krasovskii functional for a DDE

Constructs a Lyapunov-Krasovskii functional for a delay differential
equation specified via
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
with the `delays` argument. The DDE is linearised around a selected
steady state to \\dx/dt = A_0 x(t) + \sum_k A_k x(t - \tau_k)\\, reduced
to the single-delay case when \\A_k = 0\\ for \\k \ge 2\\, and a
sufficient LMI condition is solved to certify global exponential
stability.

## Usage

``` r
lyapunov_krasovskii(model, x_star = NULL, epsilon = 1e-06, verbose = TRUE, ...)
```

## Arguments

- model:

  A `system_spec` object with `is_dde = TRUE`.

- x_star:

  Numeric vector, the linearisation point. If `NULL`, the model's
  default initial condition is used (frozen-system Jacobian). Newton
  refinement on the frozen system could be added but is not required for
  a linear LMI test.

- epsilon:

  Strict inequality slack used in the LMI. Defaults to `1e-6`.

- verbose:

  Logical. When `TRUE`, narrates the construction. Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function`. The `type` field is
`"krasovskii"`. The `params` field stores the frozen-system matrices
`A0` and `A1`, the delay `tau`, and the LMI solution `P` when feasible.

## Details

The LMI sufficient condition used here is the classical simple
Lyapunov-Krasovskii test \\\left\[\begin{array}{cc} A_0^\top P + P A_0 +
S & P A_1 \\ A_1^\top P & -S \end{array}\right\] \prec 0\\, which
certifies delay-independent exponential stability when feasible (Gu,
Kharitonov and Chen, 2003, Theorem 5.6). Feasibility is checked with
CVXR, which is declared as a Suggests dependency. If CVXR is unavailable
or the LMI is infeasible, the constructor returns `feasible = FALSE`
with a clear rationale.

Current scope: single delay. For multi-delay systems, the code collapses
additional delays by summing their coefficient matrices, which is a
conservative but valid upper bound for stability testing.

## References

Gu, K., Kharitonov, V. L., & Chen, J. (2003). *Stability of Time-Delay
Systems*. Birkhauser.
[doi:10.1007/978-1-4612-0039-0](https://doi.org/10.1007/978-1-4612-0039-0)

Fridman, E. (2014). *Introduction to Time-Delay Systems*. Birkhauser.
[doi:10.1007/978-3-319-09393-2](https://doi.org/10.1007/978-3-319-09393-2)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dde <- system_spec(
    rhs = list(x ~ -a * x + b * xtau),
    delays = list(xtau = list(state = "x", tau = 1)),
    state_names = "x",
    parms = list(a = 1, b = 0.2),
    init = c(x = 0)
)
lf <- lyapunov_krasovskii(dde)
print(lf); summary(lf)
} # }
```
