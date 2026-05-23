# Add observation noise to a simulation

Applies an observation error model to the attractor of a completed
simulation, producing synthetic observed data. The true trajectory is
preserved alongside the noisy observations, enabling comparison and
signal-to-noise diagnostics.

## Usage

``` r
observe(x, obs_model, variables = NULL, seed = 42)
```

## Arguments

- x:

  A `dyn_sim` object.

- obs_model:

  A list specifying the noise model with at minimum a `type` element
  (one of `"gaussian"`, `"lognormal"`, `"poisson"`, `"neg_binomial"`)
  and the relevant parameters.

- variables:

  Character vector of variable names to observe. Defaults to `NULL`,
  which observes all state variables in the model.

- seed:

  Integer seed for reproducible noise generation (default: 42).

## Value

An S3 object of class `c("observed_dyn_sim", "dyn_sim")` inheriting all
fields from the input `dyn_sim`, with additional:

- observed:

  Data frame with the same structure as `attractor`, containing observed
  (noisy) values for the specified variables

- obs_model:

  The observation model specification used

- obs_variables:

  Character vector of observed variable names

- obs_seed:

  The seed used for noise generation

## Details

The observation layer is a post-processing step applied after
integration, not part of the dynamical model itself. This separation
keeps the model specification clean and allows the same simulation to be
observed under different noise regimes without re-integrating.

Supported noise models and their parameterizations:

**Gaussian**: \\y\_{obs} = y + \varepsilon\\, where \\\varepsilon \sim
N(0, \sigma^2)\\. Parameterized by `sigma`.

**Lognormal**: \\y\_{obs} = y \cdot \exp(\varepsilon - \sigma^2/2)\\,
where \\\varepsilon \sim N(0, \sigma^2)\\. The mean correction ensures
\\E\[y\_{obs}\] = y\\. Parameterized by `sigma`. Requires \\y \> 0\\.

**Poisson**: \\y\_{obs} = \text{Pois}(\lambda \cdot y) / \lambda\\. The
scaling by `lambda` controls the effective sample size. Parameterized by
`lambda`. Requires \\y \ge 0\\.

**Negative binomial**: \\y\_{obs} = \text{NB}(\mu = \lambda \cdot y,
\text{size}) / \lambda\\. Parameterized by `lambda` and `size` (the
dispersion parameter; smaller `size` means more overdispersion).
Requires \\y \ge 0\\.

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate a Lotka-Volterra system
model <- system_spec(
    rhs = list(
        prey ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init = c(prey = 40, predator = 9)
)
sim <- dyn_sim(model, t_max = 50, discard_transient = 30)

# Add Gaussian observation noise
obs <- observe(sim, obs_model = list(type = "gaussian", sigma = 0.01))
print(obs)
summary(obs)
plot(obs)

# Lognormal noise (multiplicative)
obs_ln <- observe(sim, obs_model = list(type = "lognormal", sigma = 0.2))
plot(obs_ln)

# Poisson counting noise
obs_pois <- observe(sim, obs_model = list(type = "poisson", lambda = 100))
plot(obs_pois)
} # }
```
