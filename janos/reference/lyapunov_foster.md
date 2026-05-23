# Foster-Lyapunov drift function for a CTMC

Attempts to build a Foster-Lyapunov function for a continuous-time
Markov chain reaction network specified via
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
with `stoichiometry` and `propensities`. The construction routes through
the fluid-limit ODE \\y' = S\\ v(y)\\, where \\S\\ is the stoichiometry
matrix and \\v\\ the propensity vector evaluated at continuous state,
and recursively invokes
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
on that ODE model. If the fluid limit admits a quadratic or
gLV-logarithmic Lyapunov function, the same \\V\\ is lifted to the CTMC
generator, and the Foster inequality is verified numerically on a grid
of reachable states.

## Usage

``` r
lyapunov_foster(
  model,
  x_star = NULL,
  grid_radius = 2,
  n_grid = 64L,
  verbose = TRUE,
  ...
)
```

## Arguments

- model:

  A `system_spec` object with `!is.null(stoich_matrix)` (CTMC / SSA).

- x_star:

  Numeric vector, the fluid-limit equilibrium. If `NULL`, uses the
  fluid-limit model's
  [`equilibrium()`](https://robustecologies.github.io/janos/reference/equilibrium.md).

- grid_radius:

  Numeric. Scaling of the neighbourhood around \\x^\*\\ where the Foster
  inequality is checked. Defaults to 2.

- n_grid:

  Integer, number of probe points. Defaults to 64.

- verbose:

  Logical. Narrate the construction. Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function` with `type = "foster"`. When
successful, the `params` field stores the inner ODE `lyapunov_function`
object, the empirical drift exponent `epsilon_emp`, and the constant
`b_emp`.

## Details

The CTMC generator acting on a smooth \\V\\ is \\QV(x) = \sum_j a_j(x)\\
\[V(x + \nu_j) - V(x)\]\\, where \\\nu_j\\ is the \\j\\th column of the
stoichiometry matrix and \\a_j\\ the associated propensity. For large
populations the functional-central-limit approximation gives \\QV(x)
\approx \nabla V(x) \cdot f(x) + \tfrac{1}{2} \sum_j a_j \nu_j^\top
\nabla^2 V \nu_j\\ where \\f = S\\ v\\, which we recognise as the Ito
generator applied to the fluid limit plus a Laplacian-like correction.
We therefore validate the lifted \\V\\ by evaluating \\QV(x)\\ exactly
on a grid around \\x^\*\\ and asking for \\QV(x) \< 0\\ outside a small
compact set.

## References

Meyn, S. P., & Tweedie, R. L. (2009). *Markov Chains and Stochastic
Stability* (2nd ed.). Cambridge University Press.
[doi:10.1017/CBO9780511626630](https://doi.org/10.1017/CBO9780511626630)

Anderson, D. F., & Kurtz, T. G. (2015). *Stochastic Analysis of
Biochemical Systems*. Springer.
[doi:10.1007/978-3-319-16895-1](https://doi.org/10.1007/978-3-319-16895-1)

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
bd <- system_spec(
    stoichiometry = list(birth = c(X = 1L), death = c(X = -1L)),
    propensities = list(birth ~ lambda, death ~ mu * X),
    state_names = "X",
    parms = list(lambda = 2, mu = 1),
    init = c(X = 5L)
)
lf <- lyapunov_foster(bd)
print(lf); summary(lf)
} # }
```
