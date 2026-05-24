# Changelog

## janos 2.1.0

This release introduces a native stiff-ODE toolkit covering four solver
families that complement the existing Rosenbrock (Rodas3) integrator,
together with a slow-fast partition API, stiffness diagnostics and
matrix-exponential utilities. Butcher tableaux are taken verbatim from
Hairer-Wanner *Solving Ordinary Differential Equations II* (Springer
1996); each coefficient set is cross-checked against an alternative
published source.

### Added

- `solver_bdf(order = 1:5)`: fixed-order backward differentiation
  formula with modified-Newton iteration, Jacobian reuse, Gustafsson PI
  step controller and embedded order-(k-1) error estimator. Coefficients
  from Hairer-Wanner Vol. II Section III.1; cross-checked against
  Shampine-Gordon 1975.

- [`solver_esdirk()`](https://robustecologies.github.io/janos/reference/solver_esdirk.md):
  TR-BDF2 explicit-first-stage singly diagonally implicit Runge-Kutta of
  Hosea-Shampine 1996; L-stable order 2 with order-3 embedded estimator.
  Equivalent to MATLAB’s `ode23tb`.

- `solver_imex_ark(method = "ars232")`: additive implicit-explicit
  Runge-Kutta of Ascher-Ruuth-Spiteri 1997, for slow-fast partitioned
  systems where the model supplies separate `rhs_explicit` and
  `rhs_implicit` closures via the new
  [`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
  arguments of the same names.

- `solver_radau(order = 3)`: two-stage fully implicit Radau IIA of order
  3 with block-Newton iteration on the coupled stage system, A- and
  L-stable. Hairer-Wanner Vol. II Section IV.5 Table 5.5.

- [`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
  dispatches to all four new stiff solvers through the standard
  `solver$method` selector and returns the usual `dyn_sim` S3 object.
  Step-level diagnostics (`step_size`, `step_error`, `n_rhs_evals`,
  `n_jac_evals`, `n_lu_decomps`, `solver_name`) are exposed in the
  `metadata` slot. The four C++ engines accept a `verbose` flag that
  periodically prints step-counter and current-time progress.

- [`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
  gains `rhs_explicit`, `rhs_implicit`, and `jac_implicit` arguments for
  additive IMEX-RK splits. When the split is supplied without an
  explicit `rhs`, the combined RHS is synthesised as the sum so that
  non-IMEX solvers can still integrate the model.

- [`slow_fast_partition()`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md):
  detect or accept a partition of the state vector into a fast block and
  a slow block via spectral-gap eigen analysis, Jacobian sparsity graph
  components, or explicit user indices. Returns an S3 object with
  `fast_idx`, `slow_idx`, `tau_fast`, `tau_slow`, `coupling_norm` plus
  the full print/summary/plot triad.

- [`stiffness_ratio()`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md):
  returns `|lambda_max| / |lambda_min|` of the Jacobian at a reference
  state, plus spectral abscissa and field-of-values radius.

- [`gsp_reduce()`](https://robustecologies.github.io/janos/reference/gsp_reduce.md):
  geometric singular perturbation reduction at zeroth order (Tikhonov),
  with a first-order Fenichel slot kept for user-supplied closures.
  Slaves fast variables to slow via Newton on the fast-component
  residual; returns a closure for the reduced right-hand side over slow
  variables.

- `expm_pade(A, t)`: matrix exponential via Higham 2005
  scaling-and-squaring with degree-13 rational Pade approximant.

- `expm_krylov(A, v, t, m)`: Arnoldi-Krylov approximation of `exp(tA) v`
  without forming the matrix exponential, with happy-breakdown
  termination.

### Changed

- `inst/include/janos_stiff.hpp` is the shared C++ header for the stiff
  stack: it declares `ButcherTableau<S>`, `NewtonSolver` (modified
  Newton with Jacobian reuse and W-scaling), `GustafssonPIController`
  (predictive PI step controller) and `SlowFastPartition`. Available to
  any LinkingTo client.

- `src/Makevars` adds `-I../inst/include` to the package compile flags
  so the stiff header resolves at package install time without copying
  it into `src/`.

## janos 2.0.0

### Breaking changes

- `model_spec()` renamed to
  [`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
  and the S3 class `model_spec` renamed to `system_spec`. The rename
  resolves the long-standing name collision with
  [`lucifer::model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.html)
  (a Bayesian model specification) and reflects more accurately what the
  constructor builds: a deterministic or stochastic dynamical system,
  not a statistical model. All S3 methods and helpers follow the rename:
  [`print.system_spec()`](https://robustecologies.github.io/janos/reference/print.system_spec.md),
  [`summary.system_spec()`](https://robustecologies.github.io/janos/reference/summary.system_spec.md),
  [`plot.system_spec()`](https://robustecologies.github.io/janos/reference/plot.system_spec.md),
  [`print.summary.system_spec()`](https://robustecologies.github.io/janos/reference/print.summary.system_spec.md).
  The interoperability helper `model_spec_rhs_cpp()` becomes
  [`system_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/system_spec_rhs_cpp.md)
  and returns objects of class `system_spec_rhs_cpp`. There is no
  deprecation shim; downstream code calling `model_spec()` or
  dispatching on class `"model_spec"` must be updated. The constructor
  signature, semantics and returned object structure are otherwise
  unchanged.

## janos 1.10.1

### New exports

- `model_spec_rhs_cpp()` extracts per-state C++ right-hand-side
  expressions from a formula-based ODE or discrete-map `model_spec`,
  returning the strings that the internal compilers would inline (state
  references as `y[i]`, parameter references as `parms[i]`, both
  0-based). The function is the stable interoperability hook consumed by
  `wadaR::compiled_system(model = ...)` and is available to any
  downstream package that needs to reuse a `model_spec` under a
  different code-generation template. Stochastic, delayed, jump,
  spatial, switched and Markov-chain families are rejected with explicit
  messages, as are models built with `positive_states` clamps or
  function-mode `rhs`.

### Internal

- The `model_spec` object now stores discrete-map formulas in the
  `map_formulas` slot in addition to the compiled C++ source. This
  restores the slot expected by the existing `summary.model_spec()` path
  for maps and enables the symbolic round-trip exposed by
  `model_spec_rhs_cpp()`.

## janos 1.10.0

### Breaking

- The `dynamical_system()` constructor and the `dynamical_system` S3
  class are removed.
  [`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
  now accepts a `model_spec()` object as its first argument, or, as an
  ergonomic shortcut, a square numeric matrix (interpreted as a linear
  system) or a named list with elements `r` and `alpha` (interpreted as
  a generalised Lotka-Volterra system). The first argument has been
  renamed from `system` to `model`. All eight algebraic methods
  (`quadratic`, `goh`, `macarthur`, `gilpin`, `sos`, `rbf`, `massera`,
  `cpa`) consume the unified ingress internally through a plain-list
  envelope.

#### Migration from `dynamical_system()`

The four constructor variants map to the new canonical form as follows.
The three ergonomic shortcuts (matrix, `list(r, alpha)`, `model_spec`)
all produce identical Lyapunov certificates; switch freely according to
which form is most readable for the system at hand.

*Linear system `dot x = A x`*:

``` r

# Old
sys <- dynamical_system(A = A, type = "linear")
lf  <- lyapunov_function(sys)

# New (matrix shortcut; simplest when A is already in hand)
lf <- lyapunov_function(A)

# New (model_spec; preferred when the system is defined symbolically)
m <- model_spec(rhs = list(x ~ -x, y ~ -2 * y),
                state_names = c("x", "y"), parms = list(),
                init = c(x = 1, y = 1))
lf <- lyapunov_function(m)
```

*Generalised Lotka-Volterra `dot x_i = x_i (r_i + sum_j alpha_ij x_j)`*:

``` r

# Old
sys <- dynamical_system(r = r, alpha = alpha, type = "glv")
lf  <- lyapunov_function(sys)

# New (list shortcut; simplest when r and alpha are numeric)
lf <- lyapunov_function(list(r = r, alpha = alpha))

# New (model_spec; preferred when r and alpha are expressed symbolically
# as part of a larger model)
m <- model_spec(rhs = list(x ~ x * (1 - x - 0.3 * y),
                           y ~ y * (0.8 - 0.2 * x - y)),
                state_names = c("x", "y"), parms = list(),
                init = c(x = 0.5, y = 0.5))
lf <- lyapunov_function(m)
```

*Polynomial system*:

``` r

# Old (manual coefficient list, now gone)
sys <- dynamical_system(coefficients = ..., type = "polynomial")

# New: write the polynomial symbolically as a model_spec; the
# structural detector recognises polynomial RHS up to degree 6 and
# routes SOS correctly.
m <- model_spec(rhs = list(x ~ y, y ~ -x + (1 - x^2) * y),
                state_names = c("x", "y"), parms = list(),
                init = c(x = 0.5, y = 0.5))
lf <- lyapunov_function(m, x_star = c(0, 0), method = "sos", degree = 4)
```

*General nonlinear system*:

``` r

# Old (closure-based general specification)
sys <- dynamical_system(
    f = function(x) c(x[2], -x[1] - 0.5 * x[2]),
    n = 2L, type = "general"
)
lf <- lyapunov_function(sys, x_star = c(0, 0), method = "rbf")

# New: express the vector field as a model_spec. The automatic
# Jacobian and the RBF, Massera and CPA methods all consume it.
m <- model_spec(rhs = list(x ~ y, y ~ -x - 0.5 * y),
                state_names = c("x", "y"), parms = list(),
                init = c(x = 1, y = 0))
lf <- lyapunov_function(m, x_star = c(0, 0), method = "rbf")
```

The first argument was renamed `system` -\> `model`. Downstream callers
using named arguments (`lyapunov_function(system = sys)`) must rewrite
to `lyapunov_function(model = sys)` or drop the name. The bare forms
`lyapunov_function(sys)` and `lyapunov_function(A)` work unchanged.

### Added

- Full `print`/`summary`/`plot` trio on every S3 class. New methods:
  [`summary.poincare_section()`](https://robustecologies.github.io/janos/reference/summary.poincare_section.md),
  [`summary.zero_one_test()`](https://robustecologies.github.io/janos/reference/summary.zero_one_test.md),
  [`summary.correlation_dimension()`](https://robustecologies.github.io/janos/reference/summary.correlation_dimension.md),
  [`summary.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/summary.bifurcation_diagram.md),
  [`plot.equilibrium_point()`](https://robustecologies.github.io/janos/reference/plot.equilibrium_point.md),
  `summary.model_spec()`, `plot.model_spec()`,
  [`summary.noise_spec()`](https://robustecologies.github.io/janos/reference/summary.noise_spec.md),
  [`plot.noise_spec()`](https://robustecologies.github.io/janos/reference/plot.noise_spec.md),
  [`summary.solver_spec()`](https://robustecologies.github.io/janos/reference/summary.solver_spec.md)
  and
  [`plot.solver_spec()`](https://robustecologies.github.io/janos/reference/plot.solver_spec.md),
  together with paired `print.summary.*` methods where the summary
  carries a tabular report. `plot.model_spec()` dispatches on the
  detected model type: vector-field grid for 2D ODE/SDE/DDE,
  reaction-graph bipartite diagram for CTMC/RDME, equation card for
  higher-dimensional or spatial models.
  [`plot.noise_spec()`](https://robustecologies.github.io/janos/reference/plot.noise_spec.md)
  renders a two-panel realisation + periodogram (or ACF) diagnostic.
  [`plot.solver_spec()`](https://robustecologies.github.io/janos/reference/plot.solver_spec.md)
  renders the absolute-stability region for RK4, RK4(5) and Rosenbrock
  and a schematic card for event-driven and method-of-lines solvers.

### Changed

- Roxygen documentation uniformised across every exported function per
  the package contract: sentence-case `@title` on every help page,
  complete `@details` with LaTeX mathematical statements, bidirectional
  `@seealso` between constructors and their S3 trio, and full-lifecycle
  `@examples` blocks (constructor -\> print -\> summary -\> plot) for
  every S3-returning function. Added missing `@references` to
  [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
  (Keller 1977, Allgower-Georg 2003) and
  [`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
  (Poincare 1892, Guckenheimer-Holmes 1983). Added missing `@seealso` to
  [`density_landscape_2d()`](https://robustecologies.github.io/janos/reference/density_landscape_2d.md).

## janos 1.9.0

### Changed

- Every S3 `plot.*` method now returns its underlying `ggplot`/`plotly`
  object visibly and no longer calls
  [`print()`](https://rdrr.io/r/base/print.html) internally. At the top
  level, R’s auto-print still renders the figure exactly once, so the
  user-visible behaviour of `plot(x)` is unchanged; but composition with
  ggplot2 layers now works as expected. Expressions such as
  `plot(zero_one_test(y)) + ggplot2::annotate(...)` previously produced
  two figures per call (one bare from the internal
  [`print()`](https://rdrr.io/r/base/print.html), one annotated from the
  top-level auto-print); they now produce the single annotated figure.
  Affects
  [`plot.dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md),
  [`plot.ensemble_sim()`](https://robustecologies.github.io/janos/reference/plot.ensemble_sim.md),
  [`plot.phase_portrait()`](https://robustecologies.github.io/janos/reference/plot.phase_portrait.md),
  [`plot.sde_portrait()`](https://robustecologies.github.io/janos/reference/plot.sde_portrait.md),
  [`plot.dde_portrait()`](https://robustecologies.github.io/janos/reference/plot.dde_portrait.md),
  [`plot.pdmp_portrait()`](https://robustecologies.github.io/janos/reference/plot.pdmp_portrait.md),
  [`plot.map_portrait()`](https://robustecologies.github.io/janos/reference/plot.map_portrait.md),
  [`plot.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/plot.observed_dyn_sim.md),
  [`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md),
  [`plot.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md),
  [`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md),
  `plot.lyapunov_spectrum()`,
  [`plot.correlation_dimension()`](https://robustecologies.github.io/janos/reference/plot.correlation_dimension.md),
  [`plot.poincare_section()`](https://robustecologies.github.io/janos/reference/plot.poincare_section.md),
  [`plot.zero_one_test()`](https://robustecologies.github.io/janos/reference/plot.zero_one_test.md),
  [`plot.bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/plot.bifurcation_diagram.md),
  [`plot.bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/plot.bifurcation_sweep.md),
  [`plot.continuation_result()`](https://robustecologies.github.io/janos/reference/plot.continuation_result.md),
  [`plot.fp_stationary()`](https://robustecologies.github.io/janos/reference/plot.fp_stationary.md),
  [`plot.fp_potential()`](https://robustecologies.github.io/janos/reference/plot.fp_potential.md),
  [`plot.fp_kramers()`](https://robustecologies.github.io/janos/reference/plot.fp_kramers.md),
  [`plot.fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/plot.fp_stationary_2d.md),
  [`plot.qsd_estimate()`](https://robustecologies.github.io/janos/reference/plot.qsd_estimate.md),
  [`plot.rare_event_estimate()`](https://robustecologies.github.io/janos/reference/plot.rare_event_estimate.md),
  [`plot.mlmc_estimate()`](https://robustecologies.github.io/janos/reference/plot.mlmc_estimate.md),
  [`plot.spectral_analysis()`](https://robustecologies.github.io/janos/reference/plot.spectral_analysis.md),
  [`plot.sensitivity_result()`](https://robustecologies.github.io/janos/reference/plot.sensitivity_result.md)
  and
  [`plot.density_landscape_2d()`](https://robustecologies.github.io/janos/reference/plot.density_landscape_2d.md).
  The `.do_print` internal switch in
  [`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
  has been removed as no longer needed.

## janos 1.8.9

### Changed

- Roxygen for
  [`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
  now references the associated `print` and `plot` S3 methods, and those
  methods receive their own roxygen blocks folded into the shared Rd
  page via `@rdname`, so the S3 trio is cross-linked bidirectionally.

### Fixed

- `.janos_parallel_apply()` now saturates every requested core. The
  previous chunk-size default (`ceiling(n_items / (n_cores * 4))`)
  produced chunks smaller than `n_cores`, so `mclapply()` inside each
  chunk only filled a handful of workers while the rest stayed idle. For
  100 parameter values on 19 cores, chunks of 2 items dispatched to 19
  workers meant 17 workers sat empty throughout the scan. Chunk size now
  targets roughly two items per worker per chunk with a floor of
  `n_cores`, which keeps every core busy while preserving responsive Esc
  handling. Affects
  [`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
  [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
  [`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
  [`map_portrait()`](https://robustecologies.github.io/janos/reference/map_portrait.md),
  [`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md),
  [`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md),
  [`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md),
  and
  [`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md).
- Pressing Esc during a parallel scan no longer emits
  `Error in invokeRestart("muffleCondition"): no 'muffleCondition' restart found`.
  R does not register a `muffleCondition` restart for interrupt
  conditions, so the
  [`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) hook
  in the chunk loop could never silence them. The handler has been
  rewritten as a `tryCatch(..., interrupt = ...)` at the chunk boundary,
  which sets the partial-result flag, breaks out cleanly, and returns
  the completed chunks.

## janos 1.8.8

### Fixed

- [`plot.spectral_analysis()`](https://robustecologies.github.io/janos/reference/plot.spectral_analysis.md)
  now renders the log-scale power spectrum correctly. The previous
  implementation clamped non-positive smoothing artefacts from
  [`spec.pgram()`](https://rdrr.io/r/stats/spec.pgram.html) to
  `.Machine$double.xmin`, which produced a y axis spanning 300+ decades
  and collapsed the actual spectrum to a thin band near the top.
  Non-positive bins are now dropped (they are numerical noise from
  Daniell smoothing, not physical content) and the y axis uses
  [`scale_y_log10()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  with `10^k` tick labels. Subtitle and caption are now populated per
  the package plot-method contract.

## janos 1.8.7

### Added

- Cross-platform parallel dispatch for
  [`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
  [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
  [`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
  [`map_portrait()`](https://robustecologies.github.io/janos/reference/map_portrait.md),
  [`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md),
  [`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md),
  [`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md),
  and
  [`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md).
  A single internal dispatcher (`.janos_parallel_apply()`) sits behind
  every parallel call site:
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix, a transparent
  [`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  on Windows. Work is dispatched in chunks; between chunks the master
  polls for Esc. Pressing Esc returns a partial result with
  `$interrupted = TRUE` and `$metadata$n_completed / n_total` instead of
  discarding the scan. Print and plot methods annotate the partial state
  so the object is still informative. The contract replaces the previous
  Windows `n_cores <- 1` fallback in
  [`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md).
- Every parallelised entry point now exposes `parallel` and `n_cores`
  arguments.
- Portrait grid helpers (`compute_vector_field_pp()`,
  `compute_nullclines_pp()`, `find_equilibria_pp()`,
  `compute_displacement_field_mp()`, `compute_isoclines_mp()`,
  `find_fixed_points_mp()`) now parallelise the inner grid/Newton loops
  on the same cross-platform path. Typical speedup on 4 cores is 3-10x
  for `n_nullcline >= 100`.
- Windows support for the compiled OpenMP batch templates
  (`compiler_batch_ssa.R`, `compiler_batch_sde.R`,
  `compiler_batch_tau.R`). A cached probe (`.janos_has_openmp()`)
  compiles a 3-line OMP source once per session; when OMP is available
  on Windows (via Rtools 4.x), the batch path emits `-fopenmp` at
  compile/link time. When OMP is absent,
  `ensemble_sim(backend = "auto")` falls back to the R-level parallel
  path via `.janos_parallel_apply()`.
- `ProgressBar` gains a `check_every` tick that yields control to R’s
  interrupt handler so serial progress-bar loops remain responsive to
  Esc;
  [`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
  inherits this responsiveness even at large `max_points` without
  changing its pseudo-arclength stepping.

### Fixed

- `find_equilibria_pp()` no longer raises
  `no se pudo encontrar la funcion evaluator` when more than 20
  candidate equilibria are found on the search grid; the residual cap
  now calls `evaluator$eval()` instead of the bare evaluator binding.

## janos 1.8.6

### Added

- [`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
  replaces the earlier attractor-level scan. The scan is embarrassingly
  parallel through
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `n_cores = max(1, parallel::detectCores() - 1)` on Unix,
  forced to serial on Windows), turning a 300-point Lorenz sweep from a
  minutes-long wait into seconds on a modern laptop. A new
  `compute_lyapunov` switch evaluates the leading Lyapunov exponent at
  every scan point via
  [`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
  and stores it in the returned S3 object; the plot method then stacks a
  second panel below the attractor scatter through , colouring by sign
  in the janos palette (amber for positive, grey for null, deep purple
  for negative) with Greek-letter labels and axis titles rendered
  through . Family-aware defaults pick a convergent Lyapunov sub-horizon
  automatically (5000 iterations for maps, 120 time units at for flows)
  when is not supplied. An explicit toggle and a compilation warm-up in
  the parent process avoid forked workers racing on the sourceCpp cache
  for expensive ODE scans. The previous name and its S3 class have been
  removed without an alias; existing code using that name must migrate
  to .
- [`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
  now accepts discrete-map models whose underlying `rhs` exposes the
  two-argument signature `function(y, parms)` (every map built through
  `model_spec(map = ...)` was in this state). The variational loop
  detects the arity at `wrap_rhs()` time and dispatches accordingly.
  This is the change that lets
  `bifurcation_diagram(compute_lyapunov = TRUE)` handle the canonical
  Feigenbaum examples (logistic, Henon).

### Fixed

- [`spectral_analysis()`](https://robustecologies.github.io/janos/reference/spectral_analysis.md)
  no longer misclassifies broadband chaotic signals as “periodic”. The
  previous classifier compared each peak to the median of the full
  spectrum; on chaotic attractors with a tail the median underflows by
  many orders of magnitude and the peak-to-background ratio becomes
  astronomical, tripping the periodic threshold. The classifier now uses
  the normalised spectral entropy , calibrated against a pure harmonic
  oscillator (), narrowband Rossler () and broadband Lorenz ().
  Peak-to-background ratios are still reported as descriptive
  diagnostics. Axis labels changed from “cycles / year” to the generic
  “cycles per time unit” now that janos is not tied to ecological
  timescales.
- [`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
  no longer collapses to on oversampled continuous flows. The previous
  stride heuristic () left consecutive samples of a smooth trajectory
  almost perfectly correlated, violating the weak-dependence assumption
  of Gottwald and Melbourne (2009) and producing a vanishing correlation
  between the mean-square displacement and the step count. The stride is
  now chosen as the first lag at which the observable’s autocorrelation
  falls below 0.5, capped so that at least 500 post-subsample points
  survive; a user override is exposed through the new `stride` argument.
  The observable is also centred before the cumulative sums, suppressing
  the deterministic drift term in the translation variables. On Lorenz
  (, ) the test now returns and correctly reports “chaotic” rather than
  “regular”.

### Changed

- `plot(dyn_sim, type = "attractor")` and the delay-coordinate embedding
  renderer now set an explicit plotly camera so that the axis projects
  to the viewer’s left and to the right. This is the janos convention
  for three-dimensional attractors and can still be overridden by
  layering another `plotly::layout()` call on top.
- `plot(lyapunov_spectrum)` now uses the janos accent palette
  (`cs_palette`): amber for positive exponents, grey for null and deep
  purple for negative, with `theme_cs()` for consistency with the rest
  of the package. The convergence panel is likewise coloured by the
  janos palette instead of the ggplot2 defaults.

## janos 1.8.5

### Fixed

- Added a robust wrapper around
  [`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)
  that recovers from corruption of the shared sourceCpp index file under
  the janos cache directory. When many runtime compilations occur in
  rapid succession (for instance during a long vignette render that
  instantiates many `model_spec()` objects back-to-back), the Rcpp
  `cache.rds` index occasionally ends up truncated, and any subsequent
  call raises *“error reading from connection”* / *“error al leer desde
  conexion”* while loading the index. The new `safe_source_cpp()` helper
  catches that exact class of error, deletes the offending index file
  and retries the compilation once. Every compiler entry point in the
  package (`R/compiler.R`, `compiler_map.R`, `compiler_dde.R`,
  `compiler_sde.R`, `compiler_jd.R`, `compiler_implicit.R`,
  `compiler_jacobian.R`, `compiler_ssa.R`, `compiler_tau_leap.R`,
  `compiler_mlmc.R`, `compiler_rdme.R`, `compiler_rdme_graph.R`,
  `compiler_pde.R`, `compiler_pde2d.R`, `compiler_pdmp.R`,
  `compiler_adjoint.R`, `compiler_batch_ssa.R`, `compiler_batch_sde.R`,
  `compiler_batch_tau.R`) now dispatches through the wrapper.
- PDE compiler: boundary blocks now declare finite-difference stencils
  for every state variable that appears anywhere in the PDE system, not
  only for the state whose equation is being written. The old behaviour
  raised an “undeclared identifier” compilation error as soon as a PDE
  had a cross-state spatial reference such as `d2x(v)` inside the
  equation for `u`, which blocked the complex Ginzburg-Landau equation
  and any reaction-diffusion system coupled through derivatives.
- PDE compiler: the fourth-derivative operator `d4x()` is now supported
  through a centred five-point stencil, with modular indexing that
  closes the stencil cleanly on periodic domains. The interior loop
  automatically contracts to when `d4x` is active, and cells are handled
  by the periodic block. Non-periodic boundary conditions combined with
  `d4x` raise an explicit error. This enables the Kuramoto-Sivashinsky
  equation and related fourth-order PDEs from the formula interface.
- PDE compiler: guarded the CFL heuristic in `dyn_sim_pde()` against
  empty or non-numeric parameter lists so parameter-free PDEs (for
  instance the Kuramoto-Sivashinsky equation) no longer raise
  *“non-numeric argument to mathematical function”*.
- New regression tests in `tests/testthat/test-pde_cross_state_d4x.R`
  cover the complex Ginzburg-Landau cross-state configuration, the
  Kuramoto-Sivashinsky `d4x` compilation and integration, and the error
  thrown when `d4x` is combined with a Dirichlet boundary.
- `vignettes/chaotic-systems.Rmd`: Kuramoto-Sivashinsky and complex
  Ginzburg-Landau chunks re-enabled with `eval = TRUE` now that the
  compiler supports both constructs natively.

## janos 1.8.4

### Added

- New chaos-analysis module `R/analysis_chaos.R` introducing five
  exported functions usable on any `model_spec` or `dyn_sim` object.
  - [`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md):
    full Lyapunov spectrum via QR renormalisation of the variational
    flow (Benettin-Shimada), finite-difference Jacobian, returns
    exponents sorted descending, running-estimate convergence,
    Kaplan-Yorke dimension and S3
    [`print()`](https://rdrr.io/r/base/print.html),
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (types
    `"convergence"`, `"spectrum"`).
  - [`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md):
    Grassberger-Procaccia correlation integral with Theiler window,
    linear fit on the scaling region, S3 with plot of log C versus log
    epsilon.
  - [`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md):
    interpolated intersections of a trajectory with a hyperplane of the
    form ; S3 with plot of the projected scatter on the remaining
    coordinates.
  - [`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md):
    Gottwald-Melbourne 0-1 test for chaos; median K over random
    frequencies with S3 print and plot.
  - `orbit_diagram()`: attractor-level bifurcation diagram by direct
    simulation across a scanned parameter (local maxima for flows,
    direct iterates for maps), complementing
    [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
    which traces equilibrium branches.
- `plot.dyn_sim(type = "dygraph")`: interactive `dygraphs` rendering of
  time series with range selector, highlight-on-hover and epsilon
  overlay. Requires the optional `dygraphs` and `zoo` packages (added to
  Suggests).
- `plot.dyn_sim(type = "delay_embedding")`: 3D Takens delay-coordinate
  reconstruction of a scalar observable; accepts an explicit `lag` or
  infers one from the model’s declared DDE delay. `type = "attractor"`
  on a DDE with fewer than three state variables now automatically
  delegates to `"delay_embedding"`, so Mackey-Glass and similar scalar
  DDEs produce a proper attractor plot out of the box.
- `plot.dyn_sim(type = "timeseries" | "dygraph")` accepts a `vars`
  argument to restrict the plot to a subset of state variables, useful
  for forced oscillators whose auxiliary phase state is otherwise
  rendered as a trivial straight ramp.
- `vignettes/chaotic-systems.Rmd` rewritten and expanded from 16 to 29
  canonical chaotic dynamical systems with per-system enrichment
  (spectrum, dimension, Poincare, 0-1 test, orbit diagram, stroboscopic
  section, delay embedding). New systems: Chen, Lu, Shimizu-Morioka,
  Nose-Hoover thermostat, Sprott minimal jerk, Lorenz-84, Lorenz-96 in
  five dimensions, logistic map, Chirikov standard map, Lozi map, Ikeda
  map, Ikeda delay DDE, Lang-Kobayashi laser DDE, Rossler hyperchaotic
  four-dimensional, seasonally forced SIR, Lotka-Volterra with Allee
  effect. May-Leonard is reframed as the heteroclinic intransitive cycle
  it is rather than a chaotic attractor. Newton-Leipnik now shows both
  coexisting strange attractors overlaid in a single 3D figure.

## janos 1.8.3

### Added

- `model_spec()` gains a `positive_states` argument (default `FALSE`)
  that selects which state variables the compiled RK4, RK45 and
  Rosenbrock integrators clamp to non-negative after each accepted step.
  Accepts `FALSE` (no clamp; correct for any sign-free ODE), `TRUE` (all
  states; correct for population densities on the carrying simplex), a
  logical vector of length `n_states`, or a character subset of
  `state_names` for mixed systems. The resolved flag participates in the
  compilation hash, so the source cache distinguishes models that differ
  only in their clamp configuration. The `print.model_spec()` method
  reports the active clamp.

### Fixed

- The compiled RK4 and RK45 templates used to clamp every component of
  the state vector to non-negative after every accepted step,
  irrespective of the mathematical nature of the system. That silently
  destroyed the dynamics of every sign-free chaotic flow (Lorenz,
  Rossler, Thomas, Halvorsen, Aizawa, Rabinovich-Fabrikant, Van der Pol,
  Duffing, Newton-Leipnik, Chua) and invalidated the RK45 FSAL cache on
  every clamped step near the origin, producing attractors collapsed
  onto the non-negative orthant and pathological step-size dynamics. The
  clamp is now opt-in through the `positive_states` argument of
  `model_spec()`.
- `plot.dyn_sim(type = "phase")` now renders discrete maps
  (`is_map == TRUE`) as a point cloud via
  [`geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  instead of joining successive iterates with
  [`geom_path()`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  The fractal structure of the Henon attractor and of companion
  two-dimensional strange sets is now visible in the default plot.
- `plot.fp_kramers(type = "summary")` no longer aborts with *argument of
  length zero* when the input model is 2D. The 1D summary view
  dereferences `potential$x`, `potential$V` and `potential$wells`,
  fields that do not exist in the 2D pathway where the potential is
  returned as a gridded surface. The method now detects the 2D schema
  and renders a compact transition-table panel (one row per transition
  with , rate and MFPT); full 2D landscape visualisations remain the
  responsibility of
  [`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md).

## janos 1.8.2

### Fixed

- `plot.fp_kramers(type = "summary")` no longer aborts with *argument of
  length zero* when the input model is 2D. The 1D summary view
  dereferences `potential$x`, `potential$V` and `potential$wells`,
  fields that do not exist in the 2D pathway where the potential is
  returned as a gridded surface. The method now detects the 2D schema
  and renders a compact transition-table panel (one row per transition
  with , rate and MFPT); full 2D landscape visualisations remain the
  responsibility of
  [`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md).

## janos 1.8.1

### Added

- Extended
  [`plot.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md)
  with two additional views: `type = "detector_scores"` renders a
  horizontal bar plot of the structural detector scores with evidence
  annotations, and `type = "radar"` draws a radial confidence polygon
  over the seven supported families.
- Extended
  [`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md)
  with two report-level views: `type = "certificate_stack"` shows a
  three-panel summary of the algebraic certificate (positivity of V,
  negativity of dV/dt, residual), and `type = "timeline"` renders the
  wall-clock sequence of the dispatcher phases (advise, classify,
  construct, verify).
- Family-specific plot types accessible via
  [`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
  or
  [`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md):
  `"cobweb"` and `"iterate_decay"` for discrete maps;
  `"generator_field"`, `"ensemble_decay"` and `"noise_ellipse"` for
  SDEs; `"lmi_spectrum"` and `"delay_margin"` for DDEs; `"regime_lmi"`
  and `"switched_trajectory"` for PDMPs; `"drift_grid"` and
  `"fluid_vs_ctmc"` for CTMC Foster lifts; `"energy_decay"`,
  `"gradient_field_check"` and `"profile"` for reaction-diffusion PDEs.

## janos 1.8.0

### Added

- `analyse_lyapunov(model)` unified Lyapunov entry point over
  `model_spec` objects. Dispatches by system family to seven specialised
  constructors and returns an S3 `lyapunov_report` with `feasible`,
  `certificate_type` (algebraic, local_algebraic, numerical, none) and a
  reasoned rationale.
- `lyapunov_advisor(model)` pure diagnostic inspector. Classifies the
  family and subtype, runs structural detectors (linearity, gLV,
  polynomial degree, symmetry, additive noise, gradient field), emits a
  ranked plan and narrates the reasoning verbosely. S3 class with
  `print`, `summary`, `plot` methods.
- Family-specific constructors, each with its own certificate theory and
  verbose narration:
  - [`lyapunov_discrete()`](https://robustecologies.github.io/janos/reference/lyapunov_discrete.md)
    discrete Lyapunov equation for maps.
  - [`lyapunov_stochastic()`](https://robustecologies.github.io/janos/reference/lyapunov_stochastic.md)
    continuous Lyapunov equation plus the Khasminskii generator for
    SDEs.
  - [`lyapunov_krasovskii()`](https://robustecologies.github.io/janos/reference/lyapunov_krasovskii.md)
    LMI feasibility for linear(ised) DDEs via CVXR.
  - [`lyapunov_pdmp()`](https://robustecologies.github.io/janos/reference/lyapunov_pdmp.md)
    common-quadratic switched-system LMI via CVXR.
  - [`lyapunov_foster()`](https://robustecologies.github.io/janos/reference/lyapunov_foster.md)
    Foster-Lyapunov drift lifted from the fluid-limit ODE for CTMC
    reaction networks.
  - [`lyapunov_functional()`](https://robustecologies.github.io/janos/reference/lyapunov_functional.md)
    energy functional for reaction-diffusion PDEs with gradient
    reaction.
- [`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
  now accepts a `model_spec` argument and delegates to
  [`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
  transparently.
- `model_spec()` preserves the raw `pde` formulas in a new
  `pde_formulas` field for downstream PDE tooling.

## janos 1.7.0

### Added

- Lyapunov function construction module with 8 computational methods:
  quadratic (Bartels-Stewart), Goh logarithmic (1977), MacArthur
  Q-function (1969), Gilpin line-integral (1974), sum-of-squares (SOS)
  via CVXR, RBF collocation, Massera converse, and CPA piecewise affine.
  Unified entry point
  [`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
  with auto-dispatch based on system type and dimension.
- `dynamical_system()` S3 constructor for ODE specifications (linear,
  gLV, polynomial, general) with auto-generated vector fields and
  analytic Jacobians.
- C++ Bartels-Stewart solver via `arma::sylvester()` for O(n^3)
  continuous and discrete Lyapunov equations.
- Six plot types for
  [`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md):
  level_sets, landscape (3D), gradient_field, phase_portrait, V_dot, and
  comparison.
- Self-contained VL-stability check (heuristic + CVXR) for the Goh
  method.

## janos 1.6.5

### Added

- Responsive Esc/Ctrl+C interrupt handling across all 19 C++ compiler
  templates. Non-throwing `pending_interrupt()` via `R_ToplevelExec`
  replaces the old `Rcpp::checkUserInterrupt()`. Check interval reduced
  from 10,000-100,000 steps to 1,024 (bitwise AND). Partial results
  preserved and returned as valid `dyn_sim` objects with
  `interrupted = TRUE` flag. Batch/OpenMP solvers use chunked
  parallelism for safe inter-chunk interrupt checks. S3
  print/summary/plot methods indicate interrupted state.
- `TEMPLATE_VERSION` salt in compiler hash computations to invalidate
  cached C++ binaries on template changes.
- [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
  for robust parameter-sweep bifurcation analysis with Newton
  multi-guess, branch tracking, and bifurcation detection. S3 triad
  (print/summary/plot) producing publication-quality line diagrams with
  stability coloring and bifurcation point markers. Replaces fragile
  pseudo-arclength continuation for most applications.
- [`density_landscape_2d()`](https://robustecologies.github.io/janos/reference/density_landscape_2d.md)
  for empirical 2D probability density estimation from SDE ensemble
  terminal states via kernel density estimation. S3 triad with
  filled-contour and 3D plotly surface plot types.
- [`fp_stationary()`](https://robustecologies.github.io/janos/reference/fp_stationary.md)
  for computing the stationary density of the 1D Fokker-Planck equation
  via Chang-Cooper discretization, with Boltzmann comparison for
  gradient systems and spectral gap estimation.
- [`fp_potential()`](https://robustecologies.github.io/janos/reference/fp_potential.md)
  for reconstructing the drift potential and quasi-potential landscape,
  identifying wells, barriers, and barrier heights.
- [`fp_kramers_rate()`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md)
  for computing Kramers escape rates and mean first-passage times
  between metastable states, with Arrhenius scaling validation.
- S3 triads (print/summary/plot) for `fp_stationary`, `fp_potential`,
  and `fp_kramers` classes.
- Internal Fokker-Planck helpers: `build_1d_drift_evaluator()`,
  `chang_cooper_weight()`, `build_fp_operator()`, `find_extrema_1d()`.
- [`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md)
  for computing and visualizing the qualitative structure of ODE
  systems. 2D portraits (ggplot2) with vector fields, nullclines,
  equilibrium points (7 types with distinct shape/color encoding),
  streamlines, stable/unstable manifolds and trajectory overlays; 3D
  portraits (plotly) with cone vector fields, equilibrium markers and
  streamline/manifold curves. Automatic multi-equilibrium search via
  grid-seeded Newton’s method with compiled Jacobian. Manifold
  computation for saddle-type equilibria via eigenvector perturbation
  and RK4 integration. Nullcline extraction via
  [`contourLines()`](https://rdrr.io/r/grDevices/contourLines.html) on
  RHS zero-level sets. Support for high-dimensional models via state
  projection. `phase_portrait` S3 class with print, summary and plot
  methods.
- [`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md)
  for SDE systems: drift/diffusion fields, Lyapunov covariance ellipses,
  Euler-Maruyama sample paths, full deterministic skeleton analysis.
- [`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md)
  for DDE systems: frozen-system analysis (vector field, nullclines,
  equilibria, manifolds), Chebyshev collocation delay spectrum, actual
  DDE trajectory overlay.
- [`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md)
  for PDMP systems: per-regime ODE portraits (fields, nullclines,
  equilibria), faceted visualization, switching trajectory overlay with
  regime coloring.
- Shared portrait utilities in `R/analysis_portrait_utils.R`: extracted
  helpers from `phase_portrait.R` and `map_portrait.R` to avoid code
  duplication across five portrait types. Includes
  `make_drift_evaluator()` unified evaluator supporting ODE, SDE (via
  `drift_parsed`), and DDE (frozen system), `reconstruct_formulas()` to
  recover formula objects from parsed representations for Jacobian
  compilation, `fixed_point_visual_encoding()` and
  `equilibrium_plotly_symbol()` shared across portrait types.
- 67 new tests in `test-phase-portrait.R`; 16 tests for SDE portrait, 11
  for DDE portrait, 10 for PDMP portrait; 75 tests in
  `test-fokker-planck.R` covering input validation, OU exact density,
  double-well bimodal density, spectral gap, potential reconstruction,
  Kramers rate, S3 methods, Chang-Cooper weight properties, operator
  matrix properties.

### Changed

- `classify_equilibrium()`: added configurable `zero_tol` parameter and
  2D trace-determinant fallback via optional `jacobian` argument.
  Reduces “unclassified” equilibria for models with near-zero
  eigenvalues.
- `find_equilibria_pp()`: added `max_equilibria = 20` cap to prevent
  spurious equilibria from dominating phase portraits. Passes Jacobian
  to classifier for 2D fallback.
- [`plot.fp_potential()`](https://robustecologies.github.io/janos/reference/plot.fp_potential.md):
  improved all plot types with informative subtitles (well/barrier
  locations), captions explaining symbols, and grey annotation styling.
- [`plot.fp_kramers()`](https://robustecologies.github.io/janos/reference/plot.fp_kramers.md):
  transition labels now show direction arrows, rate k, and MFPT with
  vertical staggering to prevent overlap.
- Phase portrait aspect ratio: replaced `coord_fixed(ratio=1)` with
  conditional `coord_cartesian` when axis scales differ by more than 5x,
  across all five portrait types.
- `plot.dyn_sim(type = "phase")`: increased linewidth from 0.3 to 0.5,
  added title, rounded line joins.
- Refactored `phase_portrait.R`: extracted ~800 lines of helper
  functions to `analysis_portrait_utils.R`. `make_rhs_evaluator()` now
  delegates to `make_drift_evaluator()` from shared utils.
- Refactored `map_portrait.R`: removed duplicate
  `fixed_point_visual_encoding()`.
- All examples previously using `spom()` now use a Lotka-Volterra model
  (intentionally using `alpha` as a parameter name to verify the bug
  fix).
- `src/init.cpp` now contains proper `R_init_janos()` with
  `R_registerRoutines()` (previously just a header include).
- DESCRIPTION title simplified (removed “Autonomous and
  Non-Autonomous”).
- Compiled dispatch for `compiled_id` starting with `"usr_"` routes
  directly to `dyn_sim_user_compiled()` without the intermediate
  `dyn_sim_compiled()` function.

### Fixed

- [`plot.phase_portrait()`](https://robustecologies.github.io/janos/reference/plot.phase_portrait.md):
  crash when one nullcline has no visible contours in the plotting
  domain (override.aes color vector length mismatch).
- SDE portrait: fill scale conflict between diffusion heatmap
  (continuous) and equilibria markers (discrete); equilibria fill now
  set as fixed aesthetic rather than mapped scale.
- [`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md):
  fixed 2D Chang-Cooper operator (was producing uniform density due to
  incorrect assembly). Rewritten row-by-row following the working 1D
  pattern. Validated against analytical OU Gaussian.
- Bug where any model with a parameter named `alpha` was misidentified
  as “non-autonomous (seasonal forcing)” due to hardcoded detection at
  `dyn_sim.R:473`. All models now correctly report as autonomous.
- Rosenbrock (Rodas3) solver: complete rewrite with verified
  Hairer-Wanner coefficients. Previous version had wrong coefficients,
  wrong C-matrix divisor (gamma\*h instead of h), wrong update formula
  (extra h factor), and wrong stage evaluation points.
- RK45/Rosenbrock output: added linear interpolation between adaptive
  steps. Previously produced staircase (piecewise-constant) output at
  output grid points.
- Legacy \*365 year-to-day time conversion removed from all compilers
  (`compiler.R`, `compiler_implicit.R`, `dyn_sim.R`,
  `analysis_sensitivity.R`).
- 1D PDE state references: added y\[k\] to y\[k\*N+i\] replacement in
  reaction terms so multi-state PDEs read from the correct spatial grid
  point.
- Levy CMS Cauchy case: removed spurious M_PI_2 factor inside log
  argument that biased asymmetric Cauchy variates.
- Hybrid SSA/CLE timing: time now advances by full dt_step (was only
  tau_slow, missing the CLE remainder portion).
- Hybrid SSA/CLE fallback reaction: now finds the last slow reaction
  (was falling back to nr-1 which may be a fast reaction).
- fBm noise scaling: pre-scale noise matrix by dt^(H-0.5) so C++ pregen
  template correctly produces dt^H total scaling (was using sqrt(dt) for
  all Hurst values).
- RK45 FSAL invalidation: detect positivity clamp activation and force
  k1 re-evaluation on next step (clamped state differs from the state
  where k7 was computed).
- Correlated Milstein correction: use (dW\[j\]^2 - Sigma\[j,j\]\*dt)
  instead of (dW\[j\]^2 - dt) to account for non-unit diagonal of
  covariance matrix.
- PDMP Lewis-Shedler thinning: added 2x safety factor to instantaneous
  rate bound, reducing probability of invalid thinning when rates
  increase during ODE integration.
- 2D PDE: reject periodic BCs with informative error (accepted in
  validation but never generated correct C++).
- Periodicity estimation in summary: removed restrictive unique-diff
  check that failed for adaptive solvers; now uses median dt.

### Removed

- `spom()` function and all SPOM-specific code (`dyn_sim_compiled()`,
  `dyn_sim_spom()`, `is_spom()`).
- `seasonal_forcing()` and `validate_forcing_params()` from
  `R/forcing_functions.R`.
- `recurrence_analysis()` and RQA S3 methods (`print.rqa_result`,
  `summary.rqa_result`, `plot.rqa_result`).
- C++ SPOM integrator (`src/spom_ode.cpp`, 765 lines) and RQA engine
  (`src/rqa.cpp`, 294 lines).
- OU noise parameters (`stochastic`, `beta`, `noise_phi`, `noise_sigma`)
  and demographic noise parameters (`demographic`, `sigma_demo`) from
  [`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
- Non-autonomous system type detection (hardcoded `alpha` parameter
  check).
- SPOM-specific derived quantities (`B_total`, `R`) from print/summary
  methods.
- SPOM and RQA references from README.
- 27 SPOM/RQA-specific tests.

## janos 1.6.0

### Added

- Network-structured graph RDME on arbitrary topologies.
- Graph constructors:
  [`lattice_graph()`](https://robustecologies.github.io/janos/reference/lattice_graph.md),
  [`ring_graph()`](https://robustecologies.github.io/janos/reference/ring_graph.md),
  [`star_graph()`](https://robustecologies.github.io/janos/reference/star_graph.md),
  [`complete_graph()`](https://robustecologies.github.io/janos/reference/complete_graph.md),
  [`random_graph()`](https://robustecologies.github.io/janos/reference/random_graph.md).
- Graph RDME C++ compilation with Gillespie direct over nodes and edges
  and CSR offsets.
- Test suite grown to 1047 tests.

## janos 1.5.0

### Added

- Correlated noise via Cholesky decomposition.
- Levy alpha-stable SDEs with the Chambers-Mallows-Stuck algorithm.
- Fractional Brownian motion via circulant embedding (Wood-Chan) with
  Hosking fallback.
- Colored (1/f^beta) noise via FFT spectral synthesis.
- Piecewise deterministic Markov processes (PDMP) with Lewis-Shedler
  thinning.
- Noise specification layer:
  [`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
  [`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md),
  [`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
  [`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md).

## janos 1.3.1

### Changed

- Consolidation release: bug fixes, numerical stability improvements,
  compiler extensions, code deduplication.
- FSAL optimization in RK45, sparse RQA storage.
- Test suite grown to 803 tests.

## janos 1.3.0

### Added

- Ensemble simulation with two-tier parallelization: C++ OpenMP batch
  plus R-level mclapply fallback.
- Batch compilers for SSA Direct, SDE Euler-Maruyama, and tau-leap.
- Fan charts, extinction tracking, convergence diagnostics.
- Test suite grown to 718 tests.

## janos 1.2.0

### Added

- Continuous adjoint sensitivity analysis for ODE parameter gradients.
- Jump-diffusion processes (Ito SDE + Poisson jumps).
- Reaction-diffusion master equation (RDME) on 1D voxel grids.
- Test suite grown to 645 tests.

## janos 1.1.0

### Added

- Adaptive tau-leap step rejection.
- Rosenbrock (Rodas3) implicit solver with symbolic Jacobian generation.
- Numerical bifurcation continuation with fold and Hopf detection.
- 2D PDE via method of lines.
- Multi-level Monte Carlo variance reduction.
- Test suite grown to 548 tests.

## janos 1.0.0

### Added

- PDE method of lines with RK4.
- Full test suite (419 tests).
- Stable API.

## janos 0.9.0

### Added

- Discrete maps, delay differential equations, Milstein and
  Euler-Maruyama SDE solvers.
- Recurrence quantification analysis.

## janos 0.7.0

### Added

- Tau-leap, midpoint tau-leap, hybrid SSA/CLE.
- Quasi-stationary distribution estimation and rare-event probability
  estimation.

## janos 0.5.0

### Added

- Exact stochastic simulation: Direct, Next-Reaction, Modified
  Next-Reaction.
- Stoichiometry-based model specification.

## janos 0.3.0

### Added

- Formula-to-C++ expression compiler.
- Observation noise layer.

## janos 0.2.0

### Added

- `model_spec()` and
  [`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
  new API.
- Adaptive Dormand-Prince RK4(5) solver.
- Spectral analysis.

## janos 0.1.0

### Added

- Initial development release. Core simulation engine: deterministic ODE
  integration, basic stochastic forcing, and S3 infrastructure for
  simulation objects. Superseded in 0.2.0 by the `model_spec()` /
  [`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
  API.
