# Estimate rare event probability via importance sampling

Estimates the probability that a CTMC trajectory reaches a target state
(e.g., extinction) using weighted stochastic simulation (importance
sampling). Each trajectory runs a modified SSA with biased propensities,
accumulating a likelihood ratio weight. The unbiased probability
estimate is the weighted average of the indicator function.

## Usage

``` r
estimate_extinction(
  model,
  target,
  bias = NULL,
  n_runs = 1000L,
  t_max = 100,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with stoichiometry and propensities.

- target:

  A one-sided formula specifying the target condition (e.g.,
  `~ N == 0`).

- bias:

  Named numeric vector of propensity bias factors, one per reaction.
  Factors \> 1 increase the rate; \< 1 decrease it. Must be positive. If
  NULL (default), all biases are 1.0 (no importance sampling).

- n_runs:

  Number of independent trajectories (default: 1000).

- t_max:

  Maximum simulation time per trajectory (default: 100).

- seed:

  Random seed (default: 42).

- verbose:

  Logical; print progress (default: TRUE).

## Value

An S3 object of class `rare_event_estimate` containing:

- probability:

  Point estimate of the rare event probability

- se:

  Standard error of the estimate

- ci_95:

  95 percent confidence interval

- log_weights:

  Vector of log-weights for trajectories that reached the target

- reached:

  Logical vector: did each trajectory reach the target?

- ess:

  Effective sample size

- model:

  The model_spec used

- parameters:

  List of algorithm parameters

## Details

For each of `n_runs` trajectories, the algorithm runs a standard SSA but
applies a bias factor `bias` to the propensities of reactions that move
the system toward the target. The likelihood ratio weight for each step
is:

\$\$w\_{\text{step}} = \frac{a_j}{a_j^{\text{biased}}} \cdot \exp\left\[
-(a_0 - a_0^{\text{biased}}) \tau \right\]\$\$

and the total weight for a trajectory is the product of step weights. If
the target is reached, the trajectory contributes its weight to the
probability estimate; otherwise it contributes zero.

The bias vector controls which reactions are over- or under-weighted.
For extinction events, increasing the bias on death reactions and
decreasing it on birth reactions accelerates sampling of rare paths.
Setting all biases to 1.0 recovers standard SSA.

## References

Kuwahara, H., & Mura, I. (2008). An efficient and exact stochastic
simulation method to analyze rare events in biochemical systems.
*Journal of Chemical Physics*, 129, 165101.
[doi:10.1063/1.2987701](https://doi.org/10.1063/1.2987701)

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`estimate_qsd`](https://robustecologies.github.io/janos/reference/estimate_qsd.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Birth-death extinction probability
bd <- model_spec(
    stoichiometry = list(
        birth = c(N = 1L),
        death = c(N = -1L)
    ),
    propensities = list(
        birth ~ lambda * N,
        death ~ mu * N
    ),
    state_names = "N",
    parms = list(lambda = 2.0, mu = 1.0),
    init = c(N = 5)
)

# Use bias to accelerate extinction sampling
prob <- estimate_extinction(bd, target = ~ N == 0,
                             bias = c(birth = 0.5, death = 2.0),
                             n_runs = 1000)
print(prob)
summary(prob)
plot(prob)
} # }
```
