# Compute phase portrait of an ODE system

Computes and assembles the qualitative geometric structure of a
continuous dynamical system defined by a `model_spec`. The resulting
object captures the vector field, nullclines, equilibrium points
(classified by stability type), streamlines, stable and unstable
manifolds, and a solution trajectory for any two- or three-dimensional
ODE system.

## Usage

``` r
phase_portrait(
  model,
  xlim = NULL,
  ylim = NULL,
  zlim = NULL,
  states = NULL,
  parms = NULL,
  n_grid = 25L,
  n_nullcline = 200L,
  vector_field = TRUE,
  nullclines = TRUE,
  equilibria = TRUE,
  streamlines = FALSE,
  manifolds = TRUE,
  trajectories = TRUE,
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

  A `model_spec` with formula-based or function-based ODE right-hand
  side.

- xlim:

  Numeric vector of length 2 giving the domain range for the first
  selected state. If `NULL`, inferred from the model.

- ylim:

  Numeric vector of length 2 giving the domain range for the second
  selected state. If `NULL`, inferred from the model.

- zlim:

  Numeric vector of length 2 for the third state (3D only). If `NULL`,
  inferred from the model.

- states:

  Character vector of 2 or 3 state names to use. If `NULL`, uses the
  first 2 (or 3 if `zlim` is provided) states.

- parms:

  Named list of parameter values. If `NULL`, uses the model defaults.

- n_grid:

  Integer grid density for the vector field (default 25).

- n_nullcline:

  Integer grid density for nullcline computation (default 200).

- vector_field:

  Logical; compute the flow field (default TRUE).

- nullclines:

  Logical; compute nullclines (default TRUE, 2D only).

- equilibria:

  Logical; find and classify equilibria (default TRUE).

- streamlines:

  Logical; compute streamlines from boundary seeds (default FALSE).

- manifolds:

  Logical; compute stable/unstable manifolds for saddle equilibria
  (default TRUE).

- trajectories:

  Controls trajectory computation. `TRUE` (default) integrates from the
  model's initial conditions (`model$init`). A list of named numeric
  vectors (or a matrix with one row per initial condition) specifies
  custom starting points. `FALSE` or `NULL` disables trajectory
  computation.

- n_streamlines:

  Integer number of seed points for streamlines (default 12).

- streamline_length:

  Numeric maximum integration time for streamlines (default 50).

- streamline_dt:

  Numeric time step for RK4 integration of streamlines (default 0.01).

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

  Numeric maximum integration time for trajectories (default 50).

- traj_dt:

  Numeric time step for trajectory integration (default 0.01).

- discard_transient:

  Numeric fraction between 0 (inclusive) and 1 (exclusive) of each
  streamline and trajectory to discard as transient (default 0). A value
  of 0.2 removes the first 20 percent of each integrated path.

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

An S3 object of class `phase_portrait` containing:

- vector_field:

  Data frame with columns for coordinates, derivatives, and magnitude,
  or `NULL` if not computed.

- nullclines:

  List of nullcline data per state, each with contour data frames, or
  `NULL`.

- equilibria:

  Data frame of equilibrium coordinates, classification, eigenvalues,
  and stability, or `NULL`.

- streamlines:

  Data frame with streamline coordinates and IDs, or `NULL`.

- manifolds:

  List of manifold branch data (equilibrium index, type, coordinate data
  frame), or `NULL`.

- trajectories:

  Data frame of trajectory coordinates, or `NULL`.

- model:

  The input `model_spec`.

- parms:

  Parameter values used.

- states:

  Character vector of selected state names.

- xlim, ylim, zlim:

  Domain limits used.

- n_dim:

  Integer, 2 or 3.

## Details

The computation proceeds in stages, each controlled by a boolean flag or
parameter.

**Trajectory.** By default, a solution trajectory is integrated from the
model's initial conditions (`model$init`) using a fourth-order
Runge-Kutta scheme. This is the primary visual element of the phase
portrait: it shows the orbit structure in 2D (approach to limit cycles,
spirals into foci) and traces out the attractor geometry in 3D (e.g. the
Rossler or Lorenz butterfly). The `trajectories` argument controls this
behavior: `TRUE` (default) uses `model$init`, a list of named numeric
vectors provides custom initial conditions, and `FALSE` disables
trajectory computation. Integration length and resolution are controlled
by `traj_length` and `traj_dt`. When `discard_transient > 0`, the
initial fraction of each trajectory is removed so that only the
asymptotic behavior is displayed; this is essential for clean attractor
visualization in chaotic systems.

The **vector field** evaluates the right-hand side
\\\mathbf{f}(\mathbf{x})\\ on a regular grid spanning the domain,
producing direction and magnitude at each point. Arrow length in the
plot is proportional to the local flow speed \\\lVert
\mathbf{f}(\mathbf{x}) \rVert\\.

**Nullclines** are the zero-level sets of each component of
\\\mathbf{f}\\. For a two-dimensional system \\\dot{x} = f_1(x,y),\\
\dot{y} = f_2(x,y)\\, the \\x\\-nullcline is the curve where \\f_1 = 0\\
and the \\y\\-nullcline where \\f_2 = 0\\. Their intersections locate
equilibrium points. Nullclines are extracted by evaluating each
component on a fine grid and computing zero-contours via
[`grDevices::contourLines()`](https://rdrr.io/r/grDevices/contourLines.html).

**Equilibrium points** are found by a grid search followed by Newton's
method. The algorithm distributes `n_eq_grid`\\^{n}\\ initial guesses
uniformly across the domain (where \\n\\ is the number of selected
states), then applies Newton iteration \$\$\mathbf{y}\_{k+1} =
\mathbf{y}\_k -
\mathbf{J}(\mathbf{y}\_k)^{-1}\\\mathbf{f}(\mathbf{y}\_k)\$\$ using the
symbolically compiled Jacobian \\\mathbf{J}\\. Duplicates closer than
`eq_merge_tol` are merged. Each equilibrium is classified by the
eigenvalues \\\lambda_i\\ of \\\mathbf{J}\\:

- **Stable node**: all \\\lambda_i\\ real with \\\mathrm{Re}(\lambda_i)
  \< 0\\

- **Unstable node**: all \\\lambda_i\\ real with
  \\\mathrm{Re}(\lambda_i) \> 0\\

- **Stable focus**: complex \\\lambda_i\\ with \\\mathrm{Re}(\lambda_i)
  \< 0\\

- **Unstable focus**: complex \\\lambda_i\\ with
  \\\mathrm{Re}(\lambda_i) \> 0\\

- **Center**: purely imaginary eigenvalues (\\\mathrm{Re}(\lambda_i) =
  0\\)

- **Saddle**: real eigenvalues with mixed signs

- **Saddle-focus**: complex eigenvalues with mixed real-part signs (3D
  only)

**Streamlines** are auxiliary solution curves seeded from evenly spaced
points on the domain boundary and integrated forward. They illustrate
the global flow structure but are not shown by default. Integration
terminates when the trajectory exits the domain, the velocity \\\lVert
\mathbf{f} \rVert\\ drops below a threshold, or the maximum integration
time is reached.

**Manifolds** are computed for saddle-type equilibria. At each saddle,
the eigenvectors of \\\mathbf{J}\\ identify the stable and unstable
subspaces. The equilibrium is perturbed by a small displacement
\\\varepsilon \mathbf{v}\_i\\ along each eigenvector (both positive and
negative directions), then integrated forward (for unstable manifolds)
or backward in time (for stable manifolds, by negating \\\mathbf{f}\\).
For complex eigenvectors the real and imaginary parts are used as
independent perturbation directions.

For models with more than two or three state variables, the `states`
argument selects which dimensions to project onto. Unselected states are
held fixed at their initial condition values. Equilibria are found for
the full system and projected into the selected subspace.

## References

Strogatz, S. H. (2015). *Nonlinear dynamics and chaos: with applications
to physics, biology, chemistry, and engineering*. Westview Press. ISBN:
978-0-8133-4910-7.

Kuznetsov, Y. A. (2004). *Elements of applied bifurcation theory*.
Springer. ISBN: 978-0-387-21906-6.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`equilibrium`](https://robustecologies.github.io/janos/reference/equilibrium.md),
[`continuation`](https://robustecologies.github.io/janos/reference/continuation.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Lotka-Volterra predator-prey: auto-domain, trajectory, manifolds
lv <- model_spec(
    rhs = list(x ~ a * x - b * x * y, y ~ d * x * y - c * y),
    state_names = c("x", "y"),
    parms = list(a = 1.0, b = 0.1, d = 0.075, c = 1.5),
    init = c(x = 10, y = 5)
)
pp <- phase_portrait(lv)
print(pp)
summary(pp)
plot(pp)

# Lotka-Volterra competition with bistability: two stable states
# separated by the stable manifold of an interior saddle
comp <- model_spec(
    rhs = list(
        N1 ~ r1 * N1 * (1 - N1 / K1 - a12 * N2 / K1),
        N2 ~ r2 * N2 * (1 - N2 / K2 - a21 * N1 / K2)
    ),
    state_names = c("N1", "N2"),
    parms = list(r1 = 1, r2 = 1, K1 = 100, K2 = 100,
                 a12 = 1.5, a21 = 1.5),
    init = c(N1 = 10, N2 = 80)
)
pp2 <- phase_portrait(comp,
    trajectories = list(c(N1 = 10, N2 = 80),
                        c(N1 = 80, N2 = 10)),
    manifold_length = 200)
print(pp2)
summary(pp2)
plot(pp2)

# 3D Rossler attractor with transient removal
rossler <- model_spec(
    rhs = list(x ~ -(y + z),
               y ~ x + a * y,
               z ~ b + z * (x - cc)),
    state_names = c("x", "y", "z"),
    parms = list(a = 0.2, b = 0.2, cc = 5.7),
    init = c(x = 0.1, y = 0.1, z = 0.1)
)
pp3d <- phase_portrait(rossler, n_grid = 6,
                       traj_length = 300,
                       discard_transient = 0.3)
print(pp3d)
summary(pp3d)
plot(pp3d)
} # }
```
