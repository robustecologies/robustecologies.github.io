# Stochastic phase portrait of an SDE system

Computes and assembles the qualitative geometric structure of a
stochastic differential equation system defined by a `model_spec`. The
resulting object captures the deterministic skeleton (drift field,
nullclines, equilibria, manifolds), the diffusion intensity landscape,
stationary covariance ellipses at stable equilibria, and sample paths
from the Euler-Maruyama scheme. The portrait reveals how noise reshapes
the phase-space organization of the underlying deterministic dynamics.

## Usage

``` r
sde_portrait(
  model,
  xlim = NULL,
  ylim = NULL,
  zlim = NULL,
  states = NULL,
  parms = NULL,
  n_grid = 25L,
  n_nullcline = 200L,
  drift_field = TRUE,
  diffusion_field = TRUE,
  nullclines = TRUE,
  equilibria = TRUE,
  manifolds = TRUE,
  streamlines = FALSE,
  sample_paths = TRUE,
  confidence_ellipses = TRUE,
  n_streamlines = 12L,
  streamline_length = 50,
  streamline_dt = 0.01,
  manifold_eps = 0.01,
  manifold_length = 50,
  manifold_dt = 0.01,
  n_eq_grid = 10L,
  eq_tol = 1e-10,
  eq_merge_tol = 1e-06,
  n_paths = 5L,
  path_length = 50,
  path_dt = 0.01,
  discard_transient = 0,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with SDE formulation (requires `diffusion`).

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

  Integer grid density for the drift and diffusion fields (default 25).

- n_nullcline:

  Integer grid density for nullcline computation (default 200).

- drift_field:

  Logical; compute the drift vector field (default TRUE).

- diffusion_field:

  Logical; compute the diffusion intensity field (default TRUE).

- nullclines:

  Logical; compute nullclines (default TRUE, 2D only).

- equilibria:

  Logical; find and classify deterministic equilibria (default TRUE).

- manifolds:

  Logical; compute stable/unstable manifolds for saddle equilibria
  (default TRUE).

- streamlines:

  Logical; compute streamlines of the drift from boundary seeds (default
  FALSE).

- sample_paths:

  Logical; simulate Euler-Maruyama sample paths (default TRUE).

- confidence_ellipses:

  Logical; compute stationary covariance ellipses at stable equilibria
  (default TRUE, 2D only).

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

- n_paths:

  Integer number of independent Euler-Maruyama sample paths (default 5).

- path_length:

  Numeric maximum integration time for sample paths (default 50).

- path_dt:

  Numeric time step for Euler-Maruyama integration (default 0.01).

- discard_transient:

  Numeric fraction between 0 (inclusive) and 1 (exclusive) of each
  streamline and sample path to discard as transient (default 0). A
  value of 0.2 removes the first 20 percent of each path.

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

An S3 object of class `sde_portrait` containing:

- drift_field:

  Data frame with columns for coordinates, derivatives, and magnitude,
  or `NULL` if not computed.

- diffusion_field:

  Data frame with columns for coordinates and diffusion intensity, or
  `NULL` if not computed.

- nullclines:

  List of nullcline data per state, each with contour data frames, or
  `NULL`.

- equilibria:

  Data frame of equilibrium coordinates, classification, eigenvalues,
  and stability, or `NULL`.

- confidence_ellipses:

  List of data frames, one per stable equilibrium, each containing 68\\
  `NULL`.

- manifolds:

  List of manifold branch data (equilibrium index, type, coordinate data
  frame), or `NULL`.

- streamlines:

  Data frame with streamline coordinates and IDs, or `NULL`.

- sample_paths:

  Data frame of sample path coordinates with columns x, y (and z for
  3D), time, and id, or `NULL`.

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

- noise_type:

  Character string describing the noise type.

## Details

The computation proceeds in stages, each controlled by a boolean flag.

**Drift field.** The deterministic skeleton \\f(\mathbf{x})\\ is
evaluated on a regular grid spanning the domain using the drift
evaluator and `compute_vector_field_pp()`, producing direction and
magnitude at each point. Arrow length in the plot is proportional to the
local drift speed \\\lVert f(\mathbf{x}) \rVert\\.

**Diffusion intensity field.** For each grid point, the diffusion
coefficient vector \\g(\mathbf{x})\\ is evaluated and its Euclidean norm
\\\lVert g(\mathbf{x}) \rVert\\ computed. When the noise specification
includes a Cholesky correlation factor \\L\\ (from
[`correlated_noise`](https://robustecologies.github.io/janos/reference/correlated_noise.md)),
the effective diffusion vector is \\L\\g(\mathbf{x})\\ and the intensity
is \\\lVert L\\g(\mathbf{x}) \rVert\\. The intensity heatmap highlights
state-dependent noise amplification or quenching.

**Nullclines** are the zero-level sets of each drift component \\f_i =
0\\, extracted via
[`grDevices::contourLines()`](https://rdrr.io/r/grDevices/contourLines.html)
on a fine grid. Their intersections locate the equilibria of the
deterministic skeleton.

**Equilibria** of the deterministic skeleton are found by grid search
followed by Newton iteration on \\f(\mathbf{x}) = 0\\, using the
symbolically compiled Jacobian. Each equilibrium is classified by the
eigenvalues of the Jacobian as stable node, unstable node, stable focus,
unstable focus, center, saddle, or saddle-focus.

**Confidence ellipses.** At each asymptotically stable equilibrium
\\\mathbf{x}^\*\\, the stationary covariance of the linearized SDE is
obtained by solving the continuous Lyapunov equation \$\$A\\\Sigma +
\Sigma\\A^\top + Q = 0\$\$ where \\A = J(\mathbf{x}^\*)\\ is the
Jacobian of the drift and \\Q = B\\B^\top\\ with \\B =
\mathrm{diag}(g(\mathbf{x}^\*))\\ (or \\B =
L\\\mathrm{diag}(g(\mathbf{x}^\*))\\ for correlated noise). The solution
is computed via Kronecker vectorization: \$\$\mathrm{vec}(\Sigma) = -(A
\otimes I + I \otimes A)^{-1}\\ \mathrm{vec}(Q).\$\$ The 68\\
\\(\mathbf{x} - \mathbf{x}^\*)^\top \Sigma^{-1} (\mathbf{x} -
\mathbf{x}^\*) = \chi^2_2(p)\\ for \\p \in \\0.68, 0.95\\\\. These are
valid only for 2D projections and only when \\A\\ is Hurwitz (all
eigenvalues with negative real part).

**Manifolds** are computed for saddle-type equilibria of the
deterministic skeleton. At each saddle, eigenvectors of the Jacobian
identify the stable and unstable subspaces. The equilibrium is perturbed
along each eigenvector and integrated forward (unstable) or backward
(stable) in time using RK4.

**Sample paths** are generated by Euler-Maruyama integration of the full
SDE: \$\$\mathbf{y}\_{n+1} = \mathbf{y}\_n + f(\mathbf{y}\_n)\\\Delta
t + g(\mathbf{y}\_n)\\\sqrt{\Delta t}\\\mathbf{Z}\_n\$\$ where
\\\mathbf{Z}\_n \sim N(0, I)\\. For correlated noise, the increment is
rotated: \\L\\\mathrm{diag}(g(\mathbf{y}\_n))\\ \sqrt{\Delta
t}\\\mathbf{Z}\_n\\. Each of `n_paths` independent realizations starts
from the model's initial condition. Integration is wrapped in error
handling to prevent NaN or overflow from terminating the computation.

## References

Gardiner, C. W. (2009). *Stochastic methods: a handbook for the natural
and social sciences* (4th ed.). Springer. ISBN: 978-3-540-70712-7.

Kloeden, P. E. & Platen, E. (1992). *Numerical solution of stochastic
differential equations*. Springer. ISBN: 978-3-540-54062-5. DOI:
[doi:10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)
.

Risken, H. (1996). *The Fokker-Planck equation: methods of solution and
applications* (2nd ed.). Springer. ISBN: 978-3-540-61530-9.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Ornstein-Uhlenbeck on 2D with multiplicative noise
ou <- model_spec(
    rhs = list(x ~ -a * x + c * y,
               y ~ c * x - b * y),
    diffusion = list(x ~ sigma * (1 + x^2)^0.5,
                     y ~ sigma * (1 + y^2)^0.5),
    state_names = c("x", "y"),
    parms = list(a = 1.0, b = 1.5, c = 0.3, sigma = 0.3),
    init = c(x = 0.5, y = 0.5)
)

set.seed(42)
sp <- sde_portrait(ou, xlim = c(-3, 3), ylim = c(-3, 3))

# Inspect results
print(sp)
summary(sp)

# Visualize
plot(sp)
plot(sp, type = "drift")
plot(sp, type = "diffusion")
plot(sp, type = "skeleton")
plot(sp, type = "stochastic")

# Bistable SDE with additive noise
bistable <- model_spec(
    rhs = list(x ~ x - x^3 - d * y,
               y ~ -g * y + d * x),
    diffusion = list(x ~ sigma, y ~ sigma),
    state_names = c("x", "y"),
    parms = list(d = 0.5, g = 1.0, sigma = 0.2),
    init = c(x = 0.5, y = 0.0)
)

set.seed(123)
sp2 <- sde_portrait(bistable, n_paths = 10, path_length = 80)
print(sp2)
summary(sp2)
plot(sp2)
} # }
```
