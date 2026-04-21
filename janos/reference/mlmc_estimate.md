# Estimate expectations of stochastic models via multi-level Monte Carlo

Estimates \\E\[Q(X(T))\]\\ for a quantity of interest \\Q\\ evaluated at
the terminal state of a CTMC model using the Anderson-Higham multi-level
Monte Carlo (MLMC) algorithm. MLMC achieves the same accuracy as
brute-force ensemble simulation with substantially fewer samples by
exploiting a telescoping identity across tau-leap resolution levels.

## Usage

``` r
mlmc_estimate(
  model,
  t_max,
  quantity,
  tol = 1,
  max_level = 5L,
  R = 2L,
  tau_base = NULL,
  n_warmup = 100L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with stoichiometry and propensities.

- t_max:

  Terminal time at which the quantity of interest is evaluated.

- quantity:

  A function taking a named numeric vector (terminal state) and
  returning a scalar. Examples: `function(x) x["S"]` for population of
  species S, `function(x) as.numeric(x["I"] == 0)` for extinction
  indicator.

- tol:

  Target root-mean-square error (default: 1.0).

- max_level:

  Maximum number of MLMC levels (default: 5).

- R:

  Refinement factor between levels (default: 2).

- tau_base:

  Base tau for level 0. If NULL (default), uses 1 percent of t_max as
  the coarsest step.

- n_warmup:

  Number of warm-up samples per level to estimate variance and cost
  (default: 100).

- seed:

  Random seed (default: 42).

- verbose:

  Logical; print progress (default: TRUE).

## Value

An S3 object of class `mlmc_estimate` containing:

- mean:

  Point estimate of the expected quantity at terminal time

- variance:

  Estimated mean squared error of the estimator

- ci:

  95 percent confidence interval (mean +/- 1.96 \* sqrt(variance))

- level_stats:

  Data frame with per-level statistics: level, n_samples, mean_diff,
  variance, cost, n_optimal

- model:

  The model_spec used

- metadata:

  List of algorithm parameters (t_max, tol, R, etc.)

## Details

The MLMC estimator decomposes the expectation as \\E\[Q_L\] = E\[Q_0\] +
\sum\_{l=1}^{L} E\[Q_l - Q\_{l-1}\]\\, where \\Q_l\\ denotes the
quantity of interest computed from a tau-leap trajectory at level \\l\\
with step size \\\tau_l = \tau_0 / R^l\\. The variance of the difference
\\Q_l - Q\_{l-1}\\ decreases with \\l\\ because the coupled coarse and
fine trajectories share Poisson increments, so fewer samples are needed
at finer levels.

The coupling follows Algorithm 3.1 of Anderson and Higham (2012): at
each coarse step of size \\\tau\_{l-1} = R \tau_l\\, the fine trajectory
advances \\R\\ sub-steps while the coarse trajectory updates using the
aggregated Poisson counts from those sub-steps.

The optimal sample allocation minimizes total cost for a target MSE of
\\\epsilon^2\\ using the formula \\N_l = \lceil 2\epsilon^{-2} \sqrt{V_l
/ C_l} \sum\_{k} \sqrt{V_k C_k} \rceil\\, where \\V_l\\ and \\C_l\\ are
the variance and cost at level \\l\\.

## References

Anderson, D. F., & Higham, D. J. (2012). Multilevel Monte Carlo for
continuous time Markov chains, with applications in biochemical
kinetics. *Multiscale Modeling and Simulation*, 10(1), 146-179.
[doi:10.1137/110840546](https://doi.org/10.1137/110840546)

Giles, M. B. (2015). Multilevel Monte Carlo methods. *Acta Numerica*,
24, 259-328.
[doi:10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X)

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`estimate_qsd`](https://robustecologies.github.io/janos/reference/estimate_qsd.md),
[`estimate_extinction`](https://robustecologies.github.io/janos/reference/estimate_extinction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Pure death process: E[N(T)] = N0 * exp(-mu * T)
death <- model_spec(
    stoichiometry = list(death = c(N = -1L)),
    propensities = list(death ~ mu * N),
    state_names = "N",
    parms = list(mu = 0.1),
    init = c(N = 100)
)

# MLMC estimate of E[N(10)]
est <- mlmc_estimate(death, t_max = 10,
                      quantity = function(x) x["N"],
                      tol = 1.0, seed = 42)
print(est)
summary(est)
plot(est)
plot(est, type = "allocation")
} # }
```
