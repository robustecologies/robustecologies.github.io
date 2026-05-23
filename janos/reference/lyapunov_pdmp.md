# Common-quadratic Lyapunov function for a PDMP

Constructs a common quadratic Lyapunov function \\V(x) = (x - x^\*)
^\top P (x - x^\*)\\ for a piecewise deterministic Markov process
specified via
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
with the `regimes` and `transitions` arguments. Each regime is
linearised around \\x^\*\\ to obtain \\A_i = Df_i(x^\*)\\, and the LMI
\\A_i^\top P + P A_i \prec 0\\ is solved simultaneously for all modes.
When the LMI is feasible, \\V\\ certifies mean-square exponential
stability under arbitrary switching (a fortiori for the jump Markov
generator).

## Usage

``` r
lyapunov_pdmp(model, x_star = NULL, epsilon = 1e-06, verbose = TRUE, ...)
```

## Arguments

- model:

  A `system_spec` object with `is_pdmp = TRUE`.

- x_star:

  Numeric vector, the linearisation point. If `NULL`, the model's
  default initial condition is used.

- epsilon:

  Strict inequality slack for each LMI block. Defaults to `1e-6`.

- verbose:

  Logical. When `TRUE`, narrates the construction. Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function` with `type = "pdmp"`. The
`params` field stores the per-regime Jacobians `A_list` and the common
LMI solution `P`.

## Details

The constructor uses CVXR (Suggests) to pose and solve the stack of
LMIs. When CVXR is unavailable the result is returned as infeasible with
a clear diagnostic. The common-quadratic test is conservative;
mode-dependent Lyapunov functions would admit broader stability regions
but require full knowledge of the transition matrix as an LMI variable
and are out of scope for this constructor.

## References

Costa, O. L. V., Fragoso, M. D., & Marques, R. P. (2005). *Discrete-Time
Markov Jump Linear Systems*. Springer.
[doi:10.1007/b138575](https://doi.org/10.1007/b138575)

Davis, M. H. A. (1984). Piecewise-deterministic Markov processes.
*Journal of the Royal Statistical Society. Series B*, 46(3), 353-388.
[doi:10.1111/j.2517-6161.1984.tb01308.x](https://doi.org/10.1111/j.2517-6161.1984.tb01308.x)

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
pdmp <- system_spec(
    regimes = list(
        mode1 = list(x ~ -a1 * x, y ~ -b1 * y),
        mode2 = list(x ~ -a2 * x, y ~ -b2 * y)
    ),
    transitions = list(
        list(from = "mode1", to = "mode2", rate = ~ 0.5),
        list(from = "mode2", to = "mode1", rate = ~ 0.5)
    ),
    state_names = c("x", "y"),
    parms = list(a1 = 1, b1 = 1, a2 = 0.5, b2 = 2),
    init = c(x = 0, y = 0)
)
lf <- lyapunov_pdmp(pdmp)
print(lf); summary(lf)
} # }
```
