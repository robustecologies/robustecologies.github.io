# Portrait of a piecewise deterministic Markov process

Computes and assembles the qualitative geometric structure of a PDMP
(piecewise deterministic Markov process) defined by a `model_spec` with
`regimes` and `transitions`. The resulting object captures per-regime
vector fields, nullclines, equilibria (classified by stability type),
and stochastic sample trajectories colored by the active regime.

## Usage

``` r
pdmp_portrait(
  model,
  xlim = NULL,
  ylim = NULL,
  states = NULL,
  parms = NULL,
  n_grid = 25L,
  n_nullcline = 200L,
  per_regime_field = TRUE,
  per_regime_nullclines = TRUE,
  per_regime_equilibria = TRUE,
  transition_rates = TRUE,
  sample_trajectories = TRUE,
  n_trajectories = 3L,
  n_eq_grid = 10L,
  eq_tol = 1e-10,
  eq_merge_tol = 1e-06,
  traj_length = 50,
  traj_dt = 0.01,
  discard_transient = 0,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with `regimes` and `transitions`.

- xlim:

  Numeric vector of length 2 giving the domain range for the first
  selected state. If `NULL`, inferred from the model's initial
  conditions.

- ylim:

  Numeric vector of length 2 giving the domain range for the second
  selected state. If `NULL`, inferred from the model's initial
  conditions.

- states:

  Character vector of 2 state names to use. If `NULL`, uses the first 2
  states.

- parms:

  Named list of parameter values. If `NULL`, uses the model defaults.

- n_grid:

  Integer grid density for each regime's vector field (default 25).

- n_nullcline:

  Integer grid density for nullcline computation (default 200).

- per_regime_field:

  Logical; compute the vector field for each regime (default TRUE).

- per_regime_nullclines:

  Logical; compute nullclines for each regime (default TRUE).

- per_regime_equilibria:

  Logical; find and classify equilibria for each regime (default TRUE).

- transition_rates:

  Logical; reserved for future implementation of transition rate field
  visualization (default TRUE, currently ignored).

- sample_trajectories:

  Logical; simulate switching trajectories via
  [`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
  and
  [`solver_pdmp`](https://robustecologies.github.io/janos/reference/solver_pdmp.md)
  (default TRUE).

- n_trajectories:

  Integer number of sample switching trajectories to simulate (default
  3). Each uses a distinct random seed.

- n_eq_grid:

  Integer grid points per dimension for equilibrium search initial
  guesses (default 10).

- eq_tol:

  Numeric convergence tolerance for Newton iteration (default 1e-10).

- eq_merge_tol:

  Numeric distance below which two equilibria are merged (default 1e-6).

- traj_length:

  Numeric maximum simulation time for sample trajectories (default 50).

- traj_dt:

  Numeric output time step for sample trajectories (default 0.01).

- discard_transient:

  Numeric fraction between 0 (inclusive) and 1 (exclusive) of each
  trajectory to discard as transient (default 0).

- parallel:

  Logical; dispatch the independent work units across workers. On Unix
  this uses
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (fork); on Windows a PSOCK cluster is created transparently. Default:
  TRUE.

- n_cores:

  Integer or NULL; number of workers. When NULL, defaults to
  `max(1, parallel::detectCores() - 1)`. Silently clamped to 1 inside
  `R CMD check` (`_R_CHECK_LIMIT_CORES_`). The work is chunked so that
  pressing Esc returns a partial result rather than discarding the run;
  the returned S3 object carries `$interrupted` and
  `$metadata$n_completed` when this happens.

- verbose:

  Logical; print progress messages (default TRUE).

## Value

An S3 object of class `pdmp_portrait` containing:

- regime_fields:

  Named list of per-regime vector field data frames, each with columns
  x, y, dx, dy, magnitude.

- regime_nullclines:

  Named list of per-regime nullcline data, each a two-element list (one
  per state) containing contour data frames.

- regime_equilibria:

  Named list of per-regime equilibrium data frames, each with columns x,
  y, classification, stable.

- trajectories:

  Data frame with columns x, y, time, id, regime recording the switching
  trajectories, or `NULL`.

- switch_points:

  Data frame with columns x, y, time, id, regime marking the locations
  of regime transitions, or `NULL`.

- model:

  The input `model_spec`.

- parms:

  Parameter values used.

- states:

  Character vector of selected state names.

- regime_names:

  Character vector of regime names.

- xlim, ylim:

  Domain limits used.

- n_dim:

  Integer, always 2.

## Details

A PDMP evolves by following deterministic ODE dynamics within a regime
(also called a mode) and switching between regimes at random times
governed by state-dependent transition rate functions. This function
analyses each regime's ODE skeleton independently and then computes full
switching trajectories via
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
with
[`solver_pdmp`](https://robustecologies.github.io/janos/reference/solver_pdmp.md).

The computation proceeds as follows.

**Per-regime vector fields.** For each regime \\r\\, the deterministic
drift \\\mathbf{f}\_r(\mathbf{x})\\ is evaluated on a regular grid
spanning the domain. Arrow length in the plot is proportional to the
local flow speed \\\lVert \mathbf{f}\_r(\mathbf{x}) \rVert\\. Comparing
fields across regimes reveals how the qualitative dynamics change upon
switching.

**Per-regime nullclines.** Within each regime the nullclines are the
zero-level sets of each component of \\\mathbf{f}\_r\\. For a
two-dimensional system \\\dot{x} = f\_{r,1}(x,y),\\ \dot{y} =
f\_{r,2}(x,y)\\, the \\x\\-nullcline of regime \\r\\ is the curve where
\\f\_{r,1} = 0\\ and the \\y\\-nullcline where \\f\_{r,2} = 0\\. Their
intersections locate regime-specific equilibrium candidates. Nullclines
are extracted by evaluating each component on a fine grid and computing
zero-contours via
[`grDevices::contourLines()`](https://rdrr.io/r/grDevices/contourLines.html).

**Per-regime equilibria.** For each regime, equilibria
\\\mathbf{f}\_r(\mathbf{x}^\*) = \mathbf{0}\\ are found by a grid search
followed by Newton iteration. The Jacobian is estimated via central
finite differences; each equilibrium is classified by the eigenvalues
\\\lambda_i\\ of the Jacobian at \\\mathbf{x}^\*\\:

- **Stable node**: all \\\lambda_i\\ real with \\\mathrm{Re}(\lambda_i)
  \< 0\\

- **Unstable node**: all \\\lambda_i\\ real with
  \\\mathrm{Re}(\lambda_i) \> 0\\

- **Stable focus**: complex \\\lambda_i\\ with \\\mathrm{Re}(\lambda_i)
  \< 0\\

- **Unstable focus**: complex \\\lambda_i\\ with
  \\\mathrm{Re}(\lambda_i) \> 0\\

- **Center**: purely imaginary eigenvalues

- **Saddle**: real eigenvalues with mixed signs

- **Saddle-focus**: complex eigenvalues with mixed real-part signs

These equilibria belong to the individual deterministic skeletons and
need not correspond to any stationary behaviour of the full PDMP; they
serve as landmarks for interpreting the regime-specific flows.

**Sample trajectories.** Full switching trajectories are obtained by
calling
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
with
[`solver_pdmp`](https://robustecologies.github.io/janos/reference/solver_pdmp.md).
Each trajectory records both the continuous state evolution and the
active regime at every output time. Regime switches are identified from
consecutive changes in the regime label, and switch points are stored
for visual annotation. Multiple trajectories can be drawn from different
seeds to illustrate the stochastic variability of the switching
behaviour.

This function is restricted to systems with exactly two continuous state
variables.

## References

Davis, M. H. A. (1984). Piecewise-deterministic Markov processes: a
general class of non-diffusion stochastic models. *J. R. Stat. Soc. Ser.
B*, 46(3), 353-388.

Benaim, M., Le Borgne, S., Malrieu, F. & Zitt, P. A. (2012). Qualitative
properties of certain piecewise deterministic Markov processes. *Ann.
Inst. Henri Poincare Probab. Stat.*, 51(3), 1024-1057.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`solver_pdmp`](https://robustecologies.github.io/janos/reference/solver_pdmp.md),
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Gene expression with promoter switching (telegraph model)
telegraph <- model_spec(
    regimes = list(
        on  = list(x ~ k1 - d1 * x,  y ~ k2 * x - d2 * y),
        off = list(x ~ -d1 * x,       y ~ k2 * x - d2 * y)
    ),
    transitions = list(
        list(from = "on",  to = "off", rate = ~ lambda_off),
        list(from = "off", to = "on",  rate = ~ lambda_on)
    ),
    state_names = c("x", "y"),
    parms = list(k1 = 10, k2 = 5, d1 = 1, d2 = 0.5,
                 lambda_on = 0.5, lambda_off = 0.3),
    init = c(x = 5, y = 20),
    initial_regime = "on"
)

pp <- pdmp_portrait(telegraph)
print(pp)
summary(pp)
plot(pp)
plot(pp, type = "trajectories")
plot(pp, type = "regime_on")
} # }
```
