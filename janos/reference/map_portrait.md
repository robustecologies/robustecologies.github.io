# Compute portrait of a discrete-time map

Computes and assembles the qualitative geometric structure of a discrete
dynamical system defined by a `model_spec` with `map` formulas. The
resulting object captures the displacement field, isoclines, fixed
points (classified by eigenvalue modulus), stable and unstable
manifolds, orbits, and (for one-dimensional maps) a cobweb diagram.

## Usage

``` r
map_portrait(
  model,
  xlim = NULL,
  ylim = NULL,
  states = NULL,
  parms = NULL,
  n_grid = 25L,
  n_isocline = 200L,
  displacement_field = TRUE,
  isoclines = TRUE,
  fixed_points = TRUE,
  manifolds = TRUE,
  orbits = TRUE,
  orbit_scatter = FALSE,
  n_fp_grid = 10L,
  fp_tol = 1e-10,
  fp_merge_tol = 1e-06,
  n_iter = 500L,
  discard_transient = 0,
  manifold_eps = 0.01,
  manifold_n_iter = 500L,
  n_preimage_iter = 50L,
  cobweb = NULL,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with `map` formulas.

- xlim:

  Numeric vector of length 2 giving the domain range for the first
  selected state. If `NULL`, inferred from the model.

- ylim:

  Numeric vector of length 2 giving the domain range for the second
  selected state (2D maps only). If `NULL`, inferred from the model.

- states:

  Character vector of 1 or 2 state names to use. If `NULL`, uses all
  states (1D or 2D depending on the model).

- parms:

  Named list of parameter values. If `NULL`, uses the model defaults.

- n_grid:

  Integer grid density for the displacement field (default 25).

- n_isocline:

  Integer grid density for isocline computation (default 200).

- displacement_field:

  Logical; compute the displacement field (default TRUE). Ignored for 1D
  maps.

- isoclines:

  Logical; compute isoclines (default TRUE, 2D only).

- fixed_points:

  Logical; find and classify fixed points (default TRUE).

- manifolds:

  Logical; compute stable/unstable manifolds for saddle fixed points
  (default TRUE, 2D only).

- orbits:

  Logical; iterate orbits from initial conditions (default TRUE).

- orbit_scatter:

  Logical; if TRUE, plots long-run orbit points as a scatter cloud to
  reveal attractor structure (default FALSE).

- n_fp_grid:

  Integer grid points per dimension for fixed point search initial
  guesses (default 10).

- fp_tol:

  Numeric convergence tolerance for Newton iteration (default 1e-10).

- fp_merge_tol:

  Numeric distance below which two fixed points are merged (default
  1e-6).

- n_iter:

  Integer number of map iterations for orbits (default 500).

- discard_transient:

  Numeric fraction between 0 and 1 of each orbit to discard as transient
  (default 0).

- manifold_eps:

  Numeric perturbation magnitude along eigenvectors for manifold
  computation (default 0.01).

- manifold_n_iter:

  Integer number of iterations for manifold tracing (default 500).

- n_preimage_iter:

  Integer maximum Newton iterations for finding preimages in stable
  manifold computation (default 50).

- cobweb:

  Logical; compute cobweb diagram data. If `NULL` (default),
  automatically enabled for 1D maps.

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

An S3 object of class `map_portrait` containing:

- displacement_field:

  Data frame with columns for coordinates, displacements, and magnitude,
  or `NULL`.

- isoclines:

  List of isocline data per state, each with contour data frames, or
  `NULL`.

- fixed_points:

  Data frame of fixed point coordinates, classification, eigenvalue
  moduli, and stability, or `NULL`.

- manifolds:

  List of manifold branch data, or `NULL`.

- orbits:

  Data frame of orbit coordinates and iteration IDs, or `NULL`.

- orbit_scatter_data:

  Data frame of long-run scatter points, or `NULL`.

- cobweb:

  List with `curve`, `diagonal`, and `staircase` data frames, or `NULL`.

- model:

  The input `model_spec`.

- parms:

  Parameter values used.

- states:

  Character vector of selected state names.

- xlim, ylim:

  Domain limits used.

- n_dim:

  Integer, 1 or 2.

- is_cobweb:

  Logical, TRUE for 1D maps.

## Details

This function is the discrete-time analogue of
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
adapting each concept from continuous flows to iterated maps.

**Displacement field.** Instead of the velocity field
\\\mathbf{f}(\mathbf{x})\\ of an ODE, the displacement
\\\mathbf{F}(\mathbf{x}) - \mathbf{x}\\ is evaluated on a regular grid.
Each arrow indicates the "kick" that the map delivers at that point; its
length is proportional to the displacement magnitude. This replaces the
vector field of the continuous case.

**Isoclines.** For a two-dimensional map \\(x, y) \mapsto (F_1(x, y),\\
F_2(x, y))\\, the isoclines are the curves \\F_1(x, y) = x\\
(x-isocline, where the first coordinate does not change) and \\F_2(x, y)
= y\\ (y-isocline). Their intersections locate fixed points. Isoclines
are extracted as zero-contours of \\F_i - x_i\\ via
[`grDevices::contourLines()`](https://rdrr.io/r/grDevices/contourLines.html),
in direct analogy to nullclines.

**Fixed points** satisfy \\\mathbf{F}(\mathbf{y}^\*) = \mathbf{y}^\*\\,
equivalently \\G(\mathbf{y}) = \mathbf{F}(\mathbf{y}) - \mathbf{y} =
0\\. They are found by distributing `n_fp_grid`\\^{n}\\ initial guesses
across the domain and applying Newton iteration on \\G\\:
\$\$\mathbf{y}\_{k+1} = \mathbf{y}\_k -
\mathbf{J}\_G(\mathbf{y}\_k)^{-1}\\G(\mathbf{y}\_k)\$\$ where
\\\mathbf{J}\_G = \mathbf{J}\_F - \mathbf{I}\\ and \\\mathbf{J}\_F\\ is
estimated by central finite differences. Duplicates closer than
`fp_merge_tol` are merged.

**Fixed point classification** uses the eigenvalue modulus
\\\|\lambda\|\\ of the Jacobian \\\mathbf{J}\_F\\ evaluated at the fixed
point, rather than the real part used for ODEs:

- **Stable node**: all eigenvalues real with \\\|\lambda_i\| \< 1\\

- **Unstable node**: all eigenvalues real with \\\|\lambda_i\| \> 1\\

- **Stable focus**: complex eigenvalues with \\\|\lambda_i\| \< 1\\

- **Unstable focus**: complex eigenvalues with \\\|\lambda_i\| \> 1\\

- **Saddle**: real eigenvalues with mixed \\\|\lambda\| \< 1\\ and
  \\\|\lambda\| \> 1\\

- **Saddle-focus**: complex eigenvalues with mixed moduli (3D or higher
  only)

- **Non-hyperbolic**: any eigenvalue with \\\|\lambda\| = 1\\ (replaces
  "center" from the ODE case)

**Manifolds** are computed for saddle-type fixed points. The unstable
manifold is traced by perturbing along the unstable eigenvector(s) of
\\\mathbf{J}\_F\\ (those with \\\|\lambda\| \> 1\\) and iterating the
map forward. The stable manifold is traced by perturbing along stable
eigenvectors and iterating backward; since maps are not generally
invertible, each preimage is found by Newton's method on
\\\mathbf{F}(\mathbf{y}) - \mathbf{y}\_{\mathrm{target}} = 0\\.

**Orbits** are computed by iterating the map from the specified initial
conditions for `n_iter` steps, analogous to trajectory integration in
the continuous case. When `discard_transient > 0`, the initial fraction
of each orbit is removed. When `orbit_scatter = TRUE`, long orbits are
plotted as scattered points to reveal attractor structure (useful for
chaotic maps such as Henon).

**Cobweb diagram.** When the map has exactly one state variable, a
cobweb (Lamerey) diagram is produced automatically: the function graph
\\y = F(x)\\ and the diagonal \\y = x\\ are drawn, then the iteration
staircase \\(x_n, x_n) \to (x_n, x\_{n+1}) \to (x\_{n+1}, x\_{n+1}) \to
\cdots\\ is traced from each initial condition. This is the canonical
visualization for 1D map dynamics.

## References

Strogatz, S. H. (2015). *Nonlinear dynamics and chaos: with applications
to physics, biology, chemistry, and engineering*. Westview Press. ISBN:
978-0-8133-4910-7.

Alligood, K. T., Sauer, T. D. and Yorke, J. A. (1996). *Chaos: an
introduction to dynamical systems*. Springer. ISBN: 978-0-387-94677-1.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`phase_portrait`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`solver_map`](https://robustecologies.github.io/janos/reference/solver_map.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1D: Logistic map cobweb diagram
logistic <- model_spec(
    map = list(x ~ r * x * (1 - x)),
    state_names = "x",
    parms = list(r = 3.85),
    init = c(x = 0.1)
)
mp <- map_portrait(logistic, xlim = c(0, 1))
print(mp)
summary(mp)
plot(mp)

# 2D: Henon map with orbit scatter
henon <- model_spec(
    map = list(x ~ 1 - a * x^2 + y, y ~ b * x),
    state_names = c("x", "y"),
    parms = list(a = 1.4, b = 0.3),
    init = c(x = 0.1, y = 0.1)
)
mp2 <- map_portrait(henon, orbit_scatter = TRUE,
                    n_iter = 5000, discard_transient = 0.1,
                    manifolds = FALSE)
print(mp2)
summary(mp2)
plot(mp2)
} # }
```
