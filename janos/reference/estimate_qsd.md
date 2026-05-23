# Quasi-stationary distribution via Fleming-Viot particles

Estimates the quasi-stationary distribution (QSD) of a CTMC model using
the Fleming-Viot particle approximation. The algorithm maintains a
population of N particles evolving via exact SSA; whenever a particle
reaches the absorbing state, it is replaced by a copy of a randomly
selected surviving particle. After a burn-in period, particle states are
sampled to construct the empirical QSD.

## Usage

``` r
estimate_qsd(
  model,
  absorbing,
  n_particles = 100L,
  t_max = 1000,
  burn_in = 0.5,
  sample_interval = 1,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec` with stoichiometry and propensities.

- absorbing:

  A one-sided formula specifying the absorbing condition (e.g.,
  `~ I == 0`).

- n_particles:

  Number of Fleming-Viot particles (default: 100).

- t_max:

  Maximum simulation time per particle synchronization interval
  (default: 1000).

- burn_in:

  Fraction of t_max to discard as burn-in (default: 0.5).

- sample_interval:

  Time interval between QSD samples (default: 1.0).

- seed:

  Random seed (default: 42).

- verbose:

  Logical; print progress (default: TRUE).

## Value

An S3 object of class `qsd_estimate` containing:

- samples:

  Matrix of sampled particle states (rows = samples, cols = state
  variables)

- marginals:

  List of marginal empirical distributions per state

- n_replacements:

  Total number of particle replacements (absorptions)

- absorption_rate:

  Estimated absorption rate per unit time

- model:

  The system_spec used

- parameters:

  List of algorithm parameters

## Details

The QSD is the limiting distribution of the process conditioned on
non-absorption. For a birth-death chain with birth rate \\\lambda\\ and
death rate \\\mu \> \lambda\\, the QSD of the population size is
geometric: \\P(N = k) = (1 - \rho) \rho^k\\ where \\\rho = \lambda /
\mu\\.

The Fleming-Viot particle system approximation converges to the true QSD
as `n_particles` increases. The absorbing condition is specified as a
formula of the form `~ condition`, where the condition is an expression
involving state variable names and comparison operators. Examples:
`~ I == 0`, `~ N <= 0`, `~ S + I == 0`.

## References

Ferrari, P. A., & Maric, N. (2007). Quasi stationary distributions and
Fleming-Viot processes. *Electronic Journal of Probability*, 12,
684-702. [doi:10.1214/EJP.v12-415](https://doi.org/10.1214/EJP.v12-415)

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`estimate_extinction`](https://robustecologies.github.io/janos/reference/estimate_extinction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Birth-death process
bd <- system_spec(
    stoichiometry = list(
        birth = c(N = 1L),
        death = c(N = -1L)
    ),
    propensities = list(
        birth ~ lambda * N,
        death ~ mu * N
    ),
    state_names = "N",
    parms = list(lambda = 0.5, mu = 1.0),
    init = c(N = 10)
)

# Estimate QSD
qsd <- estimate_qsd(bd, absorbing = ~ N == 0,
                     n_particles = 100, t_max = 1000)
print(qsd)
summary(qsd)
plot(qsd)
} # }
```
