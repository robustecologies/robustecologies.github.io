# Phase portrait of a delay differential equation system

Computes and assembles the qualitative geometric structure of a delay
differential equation (DDE) defined by a `system_spec` with delay terms.
The resulting object captures the frozen vector field, nullclines,
equilibrium points (classified under the frozen-system approximation),
an optional delay-induced characteristic spectrum, stable and unstable
manifolds, and actual DDE trajectories.

## Usage

``` r
dde_portrait(
  model,
  xlim = NULL,
  ylim = NULL,
  states = NULL,
  parms = NULL,
  n_grid = 25L,
  n_nullcline = 200L,
  frozen_field = TRUE,
  nullclines = TRUE,
  equilibria = TRUE,
  manifolds = TRUE,
  trajectories = TRUE,
  streamlines = FALSE,
  delay_spectrum = FALSE,
  n_spectrum_roots = 20L,
  n_streamlines = 12L,
  streamline_length = 50,
  streamline_dt = 0.01,
  manifold_eps = 0.01,
  manifold_length = 50,
  manifold_dt = 0.01,
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

  A `system_spec` with DDE formulas (created with `delays` argument).

- xlim:

  Numeric vector of length 2 giving the domain range for the first
  selected state. If `NULL`, inferred from trajectory extent.

- ylim:

  Numeric vector of length 2 giving the domain range for the second
  selected state. If `NULL`, inferred from trajectory extent.

- states:

  Character vector of exactly 2 state names to project onto. If `NULL`,
  uses the first 2 states.

- parms:

  Named list of parameter values. If `NULL`, uses the model defaults.

- n_grid:

  Integer grid density for the frozen vector field (default 25).

- n_nullcline:

  Integer grid density for nullcline computation (default 200).

- frozen_field:

  Logical; compute the frozen vector field (default TRUE).

- nullclines:

  Logical; compute nullclines of the frozen system (default TRUE).

- equilibria:

  Logical; find and classify equilibria of the frozen system (default
  TRUE).

- manifolds:

  Logical; compute stable/unstable manifolds for frozen saddle
  equilibria (default TRUE).

- trajectories:

  Logical; compute actual DDE trajectories via `dyn_sim` with
  `solver_dde` (default TRUE).

- streamlines:

  Logical; compute streamlines of the frozen system from boundary seeds
  (default FALSE).

- delay_spectrum:

  Logical; compute the characteristic spectrum at each equilibrium via
  Chebyshev collocation (default FALSE).

- n_spectrum_roots:

  Integer number of Chebyshev collocation points for the delay spectrum
  computation (default 20).

- n_streamlines:

  Integer number of seed points for streamlines (default 12).

- streamline_length:

  Numeric maximum integration time for streamlines (default 50).

- streamline_dt:

  Numeric time step for streamline integration (default 0.01).

- manifold_eps:

  Numeric perturbation magnitude along eigenvectors for manifold
  computation (default 0.01).

- manifold_length:

  Numeric maximum integration time for manifolds (default 50).

- manifold_dt:

  Numeric time step for manifold integration (default 0.01).

- n_eq_grid:

  Integer grid points per dimension for equilibrium search initial
  guesses (default 10).

- eq_tol:

  Numeric convergence tolerance for Newton iteration (default 1e-10).

- eq_merge_tol:

  Numeric distance below which two equilibria are merged (default 1e-6).

- traj_length:

  Numeric maximum simulation time for DDE trajectories (default 50).

- traj_dt:

  Numeric time step for DDE trajectory integration (default 0.01).

- discard_transient:

  Numeric fraction between 0 (inclusive) and 1 (exclusive) of each
  trajectory to discard as transient (default 0). A value of 0.2 removes
  the first 20 percent of each trajectory.

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

An S3 object of class `dde_portrait` containing:

- frozen_field:

  Data frame with columns for coordinates, derivatives, and magnitude
  from the frozen system, or `NULL` if not computed.

- nullclines:

  List of nullcline data per state, each with contour data frames, or
  `NULL`.

- equilibria:

  Data frame of equilibrium coordinates, frozen classification, and
  stability, or `NULL`.

- delay_spectrum:

  List of spectral data per equilibrium, each containing the approximate
  characteristic roots, delay Jacobians, and delay-corrected stability
  assessment, or `NULL`.

- manifolds:

  List of manifold branch data from the frozen system (equilibrium
  index, type, coordinate data frame), or `NULL`.

- streamlines:

  Data frame of streamline coordinates and IDs from the frozen system,
  or `NULL`.

- trajectories:

  Data frame of actual DDE trajectory coordinates, or `NULL`.

- model:

  The input `system_spec`.

- parms:

  Parameter values used.

- states:

  Character vector of selected state names.

- delays:

  Named list of delay specifications from the model.

- xlim, ylim:

  Domain limits used.

- n_dim:

  Integer, always 2 for DDE portraits.

## Details

Delay differential equations of the form \$\$\dot{\mathbf{x}}(t) =
\mathbf{f}\bigl(\mathbf{x}(t),\\ \mathbf{x}(t - \tau_1),\\ \ldots,\\
\mathbf{x}(t - \tau_k)\bigr)\$\$ possess infinite-dimensional state
spaces and correspondingly richer dynamics than their ODE counterparts.
A full phase-plane analysis is therefore not possible in the same sense
as for ODEs, but useful approximations can be built by combining the
"frozen" (zero-delay) skeleton with delay-specific spectral information.
This function produces such a combined portrait.

**Frozen vector field.** The frozen system replaces every delayed
variable with the current state, reducing the DDE to the ODE
\\\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}, \mathbf{x}, \ldots,
\mathbf{x})\\. The vector field of this ODE is evaluated on a regular
grid spanning the domain, producing direction and magnitude at each
point. This approximation captures the instantaneous flow structure and
is exact at equilibria (where delayed and current states coincide).

**Nullclines** are the zero-level sets of each component of the frozen
right-hand side. Their intersections locate candidate equilibrium
points. Because the frozen system shares its equilibria with the full
DDE (any fixed point satisfies \\\mathbf{x}^\* = \mathbf{x}^\*(t -
\tau)\\), these nullclines are geometrically exact for equilibrium
location, even though the flow field between equilibria is approximate.

**Equilibrium points** are found by Newton iteration on the frozen
system, using the same grid-search strategy as
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md).
Each equilibrium is classified by the eigenvalues of the frozen Jacobian
\\\mathbf{J}\_{\mathrm{frozen}}\\. All classifications carry the label
"(frozen)" in summary output to emphasize that the true stability of a
DDE equilibrium depends on the full characteristic spectrum, not merely
the frozen eigenvalues. A frozen-stable equilibrium can be destabilized
by delay, and conversely a frozen-unstable equilibrium can (rarely) be
stabilized.

**Delay characteristic spectrum** (optional, enabled by
`delay_spectrum = TRUE`). For each equilibrium, the function
approximates the leading roots of the characteristic quasi-polynomial
\$\$\det\bigl(s\mathbf{I} - \mathbf{A}\_0 - \sum\_{k} \mathbf{A}\_k \\
e^{-s\tau_k}\bigr) = 0\$\$ where \\\mathbf{A}\_0\\ is the instantaneous
Jacobian and \\\mathbf{A}\_k\\ is the Jacobian contribution from the
\\k\\-th delay term. The roots are computed via pseudospectral
discretization of the infinitesimal generator on Chebyshev nodes,
following Breda, Maset, and Vermiglio (2005). The delay interval
\\\[-\tau, 0\]\\ is discretized with `n_spectrum_roots` collocation
points, yielding a companion matrix whose eigenvalues approximate the
characteristic roots. This method is accurate for the rightmost
(stability-determining) roots but becomes less reliable for roots deep
in the left half-plane. Results are labeled "(approximate)" throughout.

**Manifolds** are computed for saddle-type equilibria of the frozen
system, using the same eigenvector perturbation and RK4 integration
approach as
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md).
Because these manifolds belong to the frozen ODE rather than the full
DDE, they provide qualitative guidance on the separatrix geometry but
should not be interpreted as exact invariant manifolds of the DDE.

**Trajectories** are computed by running the actual DDE solver
(`dyn_sim` with `solver_dde`) from the model's initial conditions.
Unlike the frozen-system components, these trajectories incorporate the
full delay dynamics and accurately represent the system's time
evolution. They are the primary visual element of the portrait and will
show phenomena absent from the frozen system, such as delay-induced
oscillations and limit cycles.

**Limitations.** This analysis is restricted to two-dimensional state
projections. The frozen-system analysis provides correct equilibrium
locations but only approximate stability information; enable
`delay_spectrum` for a more reliable stability assessment. The Chebyshev
collocation approach handles single-delay and multi-delay systems but
accuracy degrades for very large delay values or highly stiff systems.

## References

Erneux, T. (2009). *Applied delay differential equations*. Springer.
ISBN: 978-0-387-74371-4.

Breda, D., Maset, S. and Vermiglio, R. (2005). Pseudospectral
differencing methods for characteristic roots of delay differential
equations. *SIAM Journal on Scientific Computing*, 27(2), 482-495.
[doi:10.1137/030601600](https://doi.org/10.1137/030601600)

Smith, H. (2011). *An introduction to delay differential equations with
applications to the life sciences*. Springer. ISBN: 978-1-4419-7645-1.

## See also

[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`solver_dde`](https://robustecologies.github.io/janos/reference/solver_dde.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Mackey-Glass equation: delay-induced oscillations
mg <- system_spec(
    rhs = list(x ~ a * x_delay / (1 + x_delay^c) - b * x),
    delays = list(x_delay = list(state = "x", tau = 17.0)),
    state_names = "x",
    parms = list(a = 0.2, b = 0.1, c = 10.0),
    init = c(x = 0.9)
)
# Note: dde_portrait requires 2D systems; for 1D DDEs, add a
# dummy state or use dyn_sim directly.

# Two-species delayed predator-prey
pp_dde <- system_spec(
    rhs = list(
        N ~ r * N * (1 - N_delay / K) - a * N * P,
        P ~ b * N * P - d * P
    ),
    delays = list(N_delay = list(state = "N", tau = 2.0)),
    state_names = c("N", "P"),
    parms = list(r = 1.0, K = 10.0, a = 0.1, b = 0.05, d = 0.5),
    init = c(N = 5.0, P = 2.0)
)
dp <- dde_portrait(pp_dde, traj_length = 100)
print(dp)
summary(dp)
plot(dp)

# With delay spectrum analysis
dp2 <- dde_portrait(pp_dde, delay_spectrum = TRUE,
                    traj_length = 200, discard_transient = 0.3)
summary(dp2)
plot(dp2, type = "equilibria")
} # }
```
