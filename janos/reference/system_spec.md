# Specify a dynamical system model

Creates a model specification object that describes a dynamical system.
The model may be a deterministic ODE system (specified via `rhs`), a
continuous-time Markov chain (CTMC) reaction network (specified via
`stoichiometry` and `propensities`), a discrete map (specified via
`map`), a delay differential equation (DDE, specified via `rhs` plus
`delays`), or a stochastic differential equation (SDE, specified via
`rhs` plus `diffusion`), or a jump-diffusion process (specified via
`rhs` plus `diffusion` plus `jumps`). The specification provides a
single portable description that can be passed to different solvers via
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md).

## Usage

``` r
system_spec(
  rhs = NULL,
  state_names,
  parms = list(),
  init = NULL,
  positive_states = FALSE,
  constraint = NULL,
  env_noise = NULL,
  demo_noise = NULL,
  compiled_id = NULL,
  meta = list(),
  stoichiometry = NULL,
  propensities = NULL,
  map = NULL,
  delays = NULL,
  diffusion = NULL,
  jumps = NULL,
  noise = NULL,
  regimes = NULL,
  transitions = NULL,
  initial_regime = NULL,
  pde = NULL,
  spatial = NULL,
  rhs_explicit = NULL,
  rhs_implicit = NULL,
  jac_implicit = NULL
)
```

## Arguments

- rhs:

  Either a function with signature `function(t, y, parms)` returning a
  numeric vector of derivatives, or a list of two-sided formulas
  specifying the ODE symbolically. In formula mode, each formula has the
  form `state_var ~ expression`, where the left-hand side is a state
  variable name and the right-hand side is the derivative expression in
  terms of state variables, parameters, and time (`t`). Formula mode
  triggers automatic C++ compilation for high performance. For
  pre-compiled models (C++ backend), this may be `NULL` when the model
  is identified by `compiled_id`. For CTMC models using
  `stoichiometry`/`propensities`, set to `NULL`.

- state_names:

  Character vector of state variable names, in the order matching the
  state vector `y` and the derivative output of `rhs`.

- parms:

  Named list of default parameter values. These serve as defaults that
  can be overridden at simulation time.

- init:

  Default initial conditions. For ODE/SDE/DDE/map models, a named
  numeric vector of length `n_states`. For 1D PDE models, one of (a) a
  function `function(x)` returning a scalar (single-state) or a numeric
  vector of length `n_states`, (b) a named list of per-state functions
  `function(x)` each returning a scalar, (c) a named list of per-state
  numeric vectors of length `n_grid` (or length one, broadcast across
  the grid), or (d) a numeric vector of length `n_states` (broadcast) or
  `n_states * n_grid` (fully expanded).

- positive_states:

  Specifies which state variables are constrained to be non-negative.
  Accepts `FALSE` (default, no clamping; appropriate for sign-free
  systems such as Lorenz, Rossler, Thomas, Halvorsen, Aizawa,
  Rabinovich-Fabrikant, Van der Pol, Duffing, Newton-Leipnik, Chua and
  any ODE with oscillations across zero), `TRUE` (clamp every state to
  \\\ge 0\\ after each accepted step; appropriate for population
  densities in Lotka-Volterra, May-Leonard, Hastings-Powell and similar
  ecological models), or a character vector naming the subset of state
  variables to clamp (for mixed systems, e.g. populations coupled to
  temperature). The clamp, when active, is emitted inside the compiled
  C++ integrator and participates in the compilation hash.

- constraint:

  Optional function `function(y)` that projects state `y` onto the valid
  region after each integration step. Must return a numeric vector of
  the same length as `y`. Default: `NULL` (no constraint).

- env_noise:

  Optional named list describing environmental noise structure. Default:
  `NULL`.

- demo_noise:

  Optional named list describing demographic noise structure. Default:
  `NULL`.

- compiled_id:

  Optional string identifying a compiled C++ model backend. When set,
  solvers dispatch to the compiled implementation instead of calling
  `rhs` from R. Default: `NULL`.

- meta:

  Optional named list of metadata (description, references, source,
  etc.). Default: empty list.

- stoichiometry:

  Optional list of named integer vectors or an integer matrix specifying
  the stoichiometry of a CTMC reaction network. When provided as a list,
  each element is a named integer vector giving the state change for one
  reaction. When provided as a matrix, rows correspond to state
  variables and columns to reactions. Default: `NULL`.

- propensities:

  Optional list of two-sided formulas specifying reaction propensity
  functions. Each formula has the form `reaction_name ~ expression`
  where the expression may reference state variables and parameters. The
  number of propensities must match the number of reactions in
  `stoichiometry`. Default: `NULL`.

- map:

  Optional list of two-sided formulas specifying a discrete map y(n+1) =
  F(y(n)). Each formula has the form `state_var ~ expression`. Mutually
  exclusive with `rhs` and `stoichiometry`. Default: `NULL`.

- delays:

  Optional named list of delay specifications for DDEs. Each element is
  a list with `state` (character, the delayed state variable) and `tau`
  (numeric, the delay duration). Requires `rhs`. Default: `NULL`.

- diffusion:

  Optional list of two-sided formulas specifying the diffusion
  coefficients for an SDE. Each formula has the form
  `state_var ~ expression` giving g_i(x). Requires `rhs`. Default:
  `NULL`.

- jumps:

  Optional list of two-sided formulas specifying Poisson jump channels
  for a jump-diffusion process. Each formula has the form
  `state ~ list(intensity = ~ expr, size = ~ expr)` for deterministic
  jump sizes, or
  `state ~ list(intensity = ~ expr, size_distribution = "normal", size_mean = ~ expr, size_sd = ~ expr)`
  for stochastic jump sizes. Supported distributions: "normal",
  "exponential", "uniform". Requires `rhs`. Default: `NULL`.

- noise:

  Optional noise specification for SDE models. A `noise_spec` object
  created by
  [`correlated_noise`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
  [`levy_noise`](https://robustecologies.github.io/janos/reference/levy_noise.md),
  [`fbm_noise`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
  or
  [`colored_noise`](https://robustecologies.github.io/janos/reference/colored_noise.md),
  or a list of composable noise_specs (correlated + Levy only). Requires
  `diffusion`. Default: `NULL` (standard independent Gaussian noise).

- regimes:

  Optional named list of regime specifications for piecewise
  deterministic Markov processes (PDMP). Each element is a named list of
  formulas defining the ODE right-hand side within that regime. Mutually
  exclusive with `rhs`, `map`, `stoichiometry`, and `pde`. Requires
  `transitions`. Default: `NULL`.

- transitions:

  Optional list of transition specifications for PDMP models. Each
  element is a list with `from` (regime name), `to` (regime name), and
  `rate` (one-sided formula for transition intensity). Requires
  `regimes`. Default: `NULL`.

- initial_regime:

  Optional character string naming the initial regime for PDMP models.
  Must be one of the regime names. Default: `NULL` (uses the first
  regime).

- pde:

  Optional list of two-sided formulas specifying a 1D PDE system. Each
  formula has the form `state_var ~ expression` where the expression may
  contain spatial derivative operators `d1x(state)` and `d2x(state)` for
  first and second spatial derivatives. Mutually exclusive with `rhs`,
  `map`, and `stoichiometry`. Requires `spatial`. Default: `NULL`.

- spatial:

  Optional list specifying the spatial domain for PDE models, with
  elements `domain` (numeric vector of length 2, left and right
  boundaries), `n_grid` (integer, number of grid points, at least 3),
  and `bc` (named list of boundary conditions, one per state). Each
  boundary condition is a list with `type` ("dirichlet", "neumann", or
  "periodic") and, for Dirichlet/Neumann, `left` and `right` values.
  Default: `NULL`.

- rhs_explicit:

  Optional function `(t, y, parms) -> numeric` giving the explicit
  (non-stiff) component of an additive RHS split, consumed by
  [`solver_imex_ark`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md).
  Must be supplied together with `rhs_implicit`; the total derivative is
  the sum of the two halves. Default: `NULL`.

- rhs_implicit:

  Optional function `(t, y, parms) -> numeric` giving the implicit
  (stiff) component of the additive RHS split. Default: `NULL`.

- jac_implicit:

  Optional function `(t, y, parms) -> matrix` or a constant numeric
  matrix supplying the Jacobian of the implicit half. If absent, a
  forward-difference Jacobian of `rhs_implicit` is used. Default:
  `NULL`.

## Value

An S3 object of class `system_spec` containing:

- rhs:

  The right-hand side function (NULL for CTMC models)

- state_names:

  Character vector of state variable names

- n_states:

  Integer number of state variables

- parms:

  Named list of default parameters

- init:

  Named numeric vector of default initial conditions

- constraint:

  Constraint projection function or NULL

- env_noise:

  Environmental noise specification or NULL

- demo_noise:

  Demographic noise specification or NULL

- compiled_id:

  Compiled backend identifier or NULL

- meta:

  Metadata list

- stoich_matrix:

  Integer stoichiometry matrix or NULL (CTMC only)

- reaction_names:

  Character vector of reaction names or NULL

- n_reactions:

  Integer number of reactions or NULL

- r_propensities:

  R fallback propensity function or NULL

## Details

For ODE models, the model is defined by its deterministic drift (the ODE
right-hand side) together with optional stochastic components. The `rhs`
function must accept arguments `(t, y, parms)` where `t` is the current
time, `y` is a named numeric vector of state variables, and `parms` is a
named list of parameters. It must return a numeric vector of derivatives
with length equal to `length(state_names)`.

For CTMC models, the system is defined by a stoichiometry matrix and
propensity functions. The `stoichiometry` argument specifies how each
reaction changes the state, while `propensities` provides formulas for
the reaction rates as functions of the current state and parameters.
Propensity formulas are compiled to C++ for high-performance stochastic
simulation via the SSA family of solvers
([`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_ssa_nrm`](https://robustecologies.github.io/janos/reference/solver_ssa_nrm.md),
[`solver_ssa_mnrm`](https://robustecologies.github.io/janos/reference/solver_ssa_mnrm.md)).

For discrete maps, the `map` argument provides a list of formulas
defining the iterated map y(n+1) = F(y(n), parms). Map formulas are
compiled to C++ and iterated with
[`solver_map`](https://robustecologies.github.io/janos/reference/solver_map.md).

For DDEs, the `delays` argument specifies delayed state lookups. Each
element is a named list with `state` (the state variable to delay) and
`tau` (the delay duration). The delay variable names become available in
the `rhs` formulas. DDEs are integrated via
[`solver_dde`](https://robustecologies.github.io/janos/reference/solver_dde.md).

For SDEs, the `diffusion` argument provides formulas for the diffusion
coefficient g(x) in the Ito SDE dx = f(x)dt + g(x)dW. Use
[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md)
or
[`solver_milstein`](https://robustecologies.github.io/janos/reference/solver_milstein.md).

For jump-diffusion processes, the `jumps` argument adds Poisson-driven
jump channels to an SDE. Each jump channel specifies a target state, an
intensity (arrival rate) formula, and either a deterministic or
stochastic jump size. The diffusion component is optional; if omitted, a
zero diffusion is created automatically (pure jump-drift process). Use
[`solver_jump_diffusion`](https://robustecologies.github.io/janos/reference/solver_jump_diffusion.md).

The `constraint` function, when supplied, is called after each
integration step to project the state back onto a valid region. It
receives the state vector and must return a corrected state vector of
the same length.

Noise structure is specified through `env_noise` and `demo_noise`, which
describe extrinsic (environmental) and intrinsic (demographic)
stochasticity respectively. These are lists containing the noise type
and its parameters. The interpretation of these fields is
solver-dependent.

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`solver_rk4`](https://robustecologies.github.io/janos/reference/solver_rk4.md),
[`solver_rk45`](https://robustecologies.github.io/janos/reference/solver_rk45.md),
[`solver_mol`](https://robustecologies.github.io/janos/reference/solver_mol.md),
[`solver_map`](https://robustecologies.github.io/janos/reference/solver_map.md),
[`solver_dde`](https://robustecologies.github.io/janos/reference/solver_dde.md),
[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
[`solver_milstein`](https://robustecologies.github.io/janos/reference/solver_milstein.md),
[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_ssa_nrm`](https://robustecologies.github.io/janos/reference/solver_ssa_nrm.md),
[`solver_ssa_mnrm`](https://robustecologies.github.io/janos/reference/solver_ssa_mnrm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Lotka-Volterra with formulas (auto-compiled to C++)
lv <- system_spec(
    rhs = list(
        prey ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init = c(prey = 40, predator = 9)
)
print(lv)

# Simulate (first call triggers C++ compilation)
result <- dyn_sim(lv, t_max = 10, discard_transient = 2,
                  solver = solver_rk4())
print(result)
summary(result)
plot(result)

# SIR epidemic as a CTMC (stochastic simulation)
sir <- system_spec(
    stoichiometry = list(
        infection = c(S = -1L, I = 1L, R = 0L),
        recovery  = c(S = 0L,  I = -1L, R = 1L)
    ),
    propensities = list(
        infection ~ beta * S * I,
        recovery  ~ gamma * I
    ),
    state_names = c("S", "I", "R"),
    parms = list(beta = 0.001, gamma = 0.1),
    init = c(S = 999, I = 1, R = 0)
)
print(sir)

# Simulate with Gillespie direct method
result <- dyn_sim(sir, t_max = 200, solver = solver_ssa_direct(),
                  discard_transient = 0)
print(result)
summary(result)
plot(result)

# 1D heat equation (PDE via method of lines)
heat <- system_spec(
    pde = list(u ~ D * d2x(u)),
    state_names = "u",
    parms = list(D = 0.01),
    spatial = list(
        domain = c(0, 1), n_grid = 101,
        bc = list(u = list(type = "dirichlet", left = 0, right = 0))
    ),
    init = function(x) sin(pi * x)
)
print(heat)

# Simulate and visualize
result <- dyn_sim(heat, t_max = 2, solver = solver_mol(dt = 0.001),
                  discard_transient = 0)
print(result)
summary(result)
plot(result)
plot(result, type = "snapshot", times = c(0, 0.5, 1, 2))
} # }
```
