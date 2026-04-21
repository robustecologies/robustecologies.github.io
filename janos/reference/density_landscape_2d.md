# Compute an empirical 2D probability density landscape from ensemble data

Estimates the joint stationary density of two state variables from the
terminal states of an ensemble simulation using kernel density
estimation. This provides a Monte Carlo alternative to solving the 2D
Fokker-Planck equation directly, and is applicable to any stochastic
model type (SDE, SSA, tau-leap, jump-diffusion, PDMP).

## Usage

``` r
density_landscape_2d(
  ensemble,
  states = NULL,
  xlim = NULL,
  ylim = NULL,
  n_grid = 100L,
  log_scale = FALSE
)
```

## Arguments

- ensemble:

  An `ensemble_sim` object (from
  [`ensemble_sim`](https://robustecologies.github.io/janos/reference/ensemble_sim.md)).

- states:

  Character vector of length 2 specifying which state variables to use.
  Defaults to the first two.

- xlim:

  Numeric vector of length 2 for the x-axis range. If `NULL`, determined
  from the data with 10 percent padding.

- ylim:

  Numeric vector of length 2 for the y-axis range. If `NULL`, determined
  from the data with 10 percent padding.

- n_grid:

  Integer; resolution of the KDE grid in each dimension. Default: 100.

- log_scale:

  Logical; if `TRUE`, return \\\log\_{10} p(x, y)\\ instead of \\p(x,
  y)\\. Default: FALSE.

## Value

An S3 object of class `density_landscape_2d` with elements `x_grid`,
`y_grid` (grid vectors), `density` (matrix), `state_names`, `log_scale`,
`n_replicates`, and `terminal_states`.

## Details

The function extracts the terminal state of each replicate from an
[`ensemble_sim`](https://robustecologies.github.io/janos/reference/ensemble_sim.md)
object and applies two-dimensional Gaussian kernel density estimation
via [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html). The
bandwidth is selected automatically unless the user provides explicit
axis limits that constrain the domain. For log-scale densities (useful
when the probability mass spans many orders of magnitude), the returned
density matrix contains \\\log\_{10} p(x, y)\\ instead of \\p(x, y)\\.

The resulting object has S3 methods for `print`, `summary`, and `plot`.
The default plot is a filled-contour heatmap with white contour lines;
`type = "surface"` produces a 3D plotly surface.

## References

Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with
S* (4th ed.). Springer. ISBN 978-0-387-21706-2.

## Examples

``` r
if (FALSE) { # \dontrun{
# SDE model
m <- model_spec(
    rhs = list(x ~ -x, y ~ -y),
    diffusion = list(x ~ 0.5, y ~ 0.5),
    state_names = c("x", "y"),
    parms = list(),
    init = c(x = 1, y = 1)
)

# Ensemble
ens <- ensemble_sim(m, n_replicates = 500, t_max = 20,
                    solver = solver_euler_maruyama(dt = 0.01))

# Density landscape
dl <- density_landscape_2d(ens)
print(dl)
summary(dl)
plot(dl)
plot(dl, type = "surface")
} # }
```
