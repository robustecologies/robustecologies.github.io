# Technical architecture of the janos simulation framework

  

## Design philosophy

janos is a general-purpose simulation framework for spatio-temporal
stochastic dynamical systems. A single model specification layer feeds a
catalogue of seventeen solver backends covering deterministic ODEs,
discrete maps, delay differential equations, partial differential
equations in one and two spatial dimensions, stochastic differential
equations with correlated, Levy alpha-stable, fractional Brownian, and
colored noise, jump-diffusion processes, piecewise deterministic Markov
processes, continuous-time Markov chains, and the reaction-diffusion
master equation on regular grids and on arbitrary graph topologies. The
architecture is designed so that a user writes the model once and swaps
solvers by changing a single argument; the downstream analysis pipeline,
ensemble simulation, observation layer, and Lyapunov subsystem all
consume the same object without modification. This section lays out the
four architectural principles that make this possible.

### Separation of model and solver

Every simulation follows a three-step pattern: define the model with
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md),
choose an integration scheme with one of the solver constructors, and
hand both to
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
The model is a declarative object; the solver is a configuration object;
the simulator is the dispatcher that binds them. A Lotka-Volterra system
integrated with classical RK4, adaptive Dormand-Prince RK4(5), implicit
Rosenbrock, or Euler-Maruyama with additive noise differs only in the
`solver` argument. The model equations never change. This separation has
two consequences: analysis tools are agnostic to how the trajectory was
generated, and a family of solvers can be compared on a fixed model
without rewriting anything.

### Formula-to-C++ compilation

The default entry point is a list of R formulas where each left-hand
side is a state variable and each right-hand side is a derivative,
propensity, map, or diffusion expression. The
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
function walks these expressions with a restricted abstract-syntax-tree
translator, validates every symbol against the known state names,
parameter names, and a whitelist of allowed functions (`exp`, `log`,
`sqrt`, trigonometric and hyperbolic families, `pow`, `abs`, `sign`,
`ifelse`, and comparison operators), and emits equivalent C++ source.
The generated source is compiled at runtime through
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html),
cached on disk under `tools::R_user_dir("janos", "cache")`, and loaded
as a native function pointer. Subsequent calls with the same formulas
reuse the cache. Compilation is therefore amortized across sessions; the
user writes R, but the inner loop runs at native speed.

### Progressive analysis

janos is designed for an iterative analysis workflow: simulate, inspect,
escalate. A trajectory carries enough metadata for downstream tools to
reconstruct the dynamics without recomputing. Phase portraits sit on the
drift field, Fokker-Planck tools sit on the SDE generator, the
bifurcation engine sits on the symbolic Jacobian, the Lyapunov advisor
sits on the model structure, and the quasi-stationary distribution and
rare-event engines sit on the CTMC propensities. Every analysis tool
returns an S3 object equipped with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods with a
uniform grammar. Escalation between tools is a matter of feeding one
output into the next, without format negotiation.

### Two-tier parallelism

Large ensembles use two parallelization tiers. The primary tier is an
OpenMP batch compiled for SSA direct, Euler-Maruyama, and adaptive
tau-leap solvers, where the entire replicate loop runs inside a single
C++ shared object with per-thread xoshiro256+ state. The secondary tier
is an R-level fallback through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
Unix or
[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)
on Windows, which covers every solver type including the implicit, PDE,
DDE, jump-diffusion, PDMP, and RDME families. The backend selection is
automatic; when OpenMP is unavailable on the build (detected once per
session by `.janos_has_openmp()`) `ensemble_sim(backend = "auto")` falls
back to the R-level path. Reproducibility is guaranteed regardless of
thread count: each replicate i consumes seed = master_seed + i - 1,
decoupling results from the number of workers.

### Parallelism and early termination

Every analysis entry point that performs independent-per-iteration work
([`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
[`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md),
[`map_portrait()`](https://robustecologies.github.io/janos/reference/map_portrait.md),
[`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md),
[`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md),
[`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md),
[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md))
accepts the pair `parallel = TRUE, n_cores = NULL` and dispatches
through a single internal helper (`.janos_parallel_apply()`). On Unix
the helper invokes
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) so
that forked children inherit compiled C++ environments; on Windows it
opens a transparent
[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)
with the model, parameters and helpers exported once. There is no
separate Windows branch in user-facing code, and no Unix-only assumption
in the implementation. The helper chunks the work so that, between
chunks, the master polls R’s interrupt flag; pressing Esc terminates the
scan and returns a partial result with `$interrupted = TRUE` and
`$metadata$n_completed / n_total` set, rather than discarding the run.
Print methods annotate the partial state, and plot methods append
`(interrupted at N/M)` to the figure caption.

For the compiled OpenMP batch templates the same contract holds at C++
level: a short non-longjmp helper built around `R_ToplevelExec` polls
for the interrupt flag between replicate chunks, so a long Lorenz
ensemble or a large tau-leap scan yields to Esc in bounded time (a few
hundred milliseconds). On Windows the batch templates link against the
Rtools OpenMP runtime when available; the session-cached probe decides
this automatically.

Portrait grids benefit the most from parallelism in practice. A 200x200
nullcline computation at `n_nullcline = 200` runs roughly 8-10x faster
on 4 cores than serially, since the inner loop is a few million RHS
evaluations. The `n_eq_grid` Newton search and the `n_grid` vector-field
evaluation share the same dispatcher and are similarly accelerated.
Bifurcation scans at 300-500 parameter values on Unix see near-linear
speedup up to the number of physical cores;
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
in parallel mode intentionally drops the warm-start between parameter
values so workers stay independent, so a run of 200 points may miss the
occasional branch reachable only by continuation from a nearby solution.
In that case, run `bifurcation_sweep(..., parallel = FALSE)` for the
last-mile branch tracing.

Inside `R CMD check` the helper silently forces `n_cores = 1` (honoring
`_R_CHECK_LIMIT_CORES_`). The ensemble OpenMP path honors the same flag
through its `n_threads` argument. This means test suites remain
single-threaded by default without changing call sites.

Layered architecture of the janos framework. Each row is a subsystem;
arrows show the data flow from user-facing definitions through the
compiler and solver engine down to the C++ runtime.

  

## Architecture overview

The package is organized around five subsystems. The **specification
layer**
([`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md),
`noise_spec()`, and the seventeen `solver_*()` constructors) provides
the declarative vocabulary. The **compiler** (`R/compiler*.R`)
translates formulas and stoichiometries into C++ templates, generates
symbolic Jacobians when needed, and routes the source through a
persistent cache with self-repair against index corruption. The
**simulation engine**
([`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md),
[`observe()`](https://robustecologies.github.io/janos/reference/observe.md))
dispatches to the appropriate backend and returns an S3 object with a
uniform shape. The **analysis toolkit** (phase portraits, Fokker-Planck,
sensitivity, bifurcation, QSD, rare events, spectral, MLMC, chaos
diagnostics) operates on model objects or on trajectories. The
**Lyapunov subsystem**
([`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
and fourteen family builders) is an orthogonal layer that extracts
algebraic or numerical stability certificates from any model.

### The model contract

A `model_spec` object is a list with a stable set of slots that
downstream tools read. The core fields are `type` (one of `ode`, `map`,
`dde`, `sde`, `jump_diffusion`, `pdmp`, `ssa`, `tau_leap`, `pde`,
`pde2d`, `rdme`), `state_names` (character vector of state variable
names), `parms` (named list of numeric parameters), `init` (named
numeric vector of initial conditions), `rhs` or `propensities` or `map`
(the formula list or function), and `compiled_id` (an MD5 hash that keys
the compilation cache). SDEs add `diffusion`; DDEs add `delays`; PDMPs
add `regimes` and `switching`; SSAs add `stoichiometry`; PDEs add
`mesh`, `bc_left`, `bc_right`, and the spatial-operator registry. The
[`print.model_spec()`](https://robustecologies.github.io/janos/reference/print.model_spec.md)
method renders a one-screen synopsis naming the type, the state
dimension, the parameter roster, the compilation status, and any
auxiliary structure (delays, diffusion, jumps, mesh).

The contract is what makes the framework composable.
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
reads `compiled_id` and dispatches to the compiled backend when present;
[`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md)
evaluates `rhs` on a spatial grid;
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
matches the `type` field against its advisor table;
[`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
calls the symbolic-Jacobian generator on the `rhs` list;
[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md)
checks whether the compiled backend supports its OpenMP batch path. None
of these tools manipulate model internals directly; they ask the
contract.

Component relationship map. Orange nodes are infrastructure, yellow
nodes are user-level objects, blue nodes are engines, green nodes are
orchestration and analysis. Dashed arrows indicate dispatch; solid
arrows indicate data flow.

### Subsystem summary

| Subsystem | Source files | Exported API | Role |
|:---|:---|:---|:---|
| Specification | model_spec.R + noise_spec.R + solvers.R | 19 (1 spec + 4 noise + 17 solvers - overlap) | Declare dynamics and noise |
| Compilation | compiler\*.R (22 files) + safe_source_cpp() | 0 (internal) | Translate formulas to native code |
| Simulation | dyn_sim.R + ensemble.R | 3 (dyn_sim, ensemble_sim, observe) | Integrate the system and return S3 objects |
| Analysis | analysis\_\*.R (15 files) | 20 (portraits, FP, continuation, QSD, rare event, sensitivity, spectral, MLMC, chaos diagnostics) | Extract structure from models and trajectories |
| Lyapunov | lyapunov\_\*.R (18 files) | 9 (analyse + advisor + 7 family constructors) | Prove stability algebraically or numerically |
| Observation | observe.R | 1 (observe) | Add measurement noise after integration |
| Utilities | graph_utils.R + utils.R + VonNeumanns_elephant.R | 5 (graph constructors, elephant curve) | Support types, helpers, demos |

janos subsystem summary. The five core subsystems account for roughly 50
R files and the entirety of the exported API surface.

  

## Model specification

A `model_spec` is the unique integration point between user intent and
runtime machinery. The constructor is polymorphic: the same function
handles eleven model types through a dispatch on which optional
arguments are supplied. Each type has its own validation path and its
own downstream compilation template.

### Formula grammar

The ODE, SDE, DDE, jump-diffusion, map, and PDE specifications use R
formulas of the form `state ~ expression`. The parser in `R/compiler.R`
enforces a closed expression language: the left-hand side must be a
single symbol matching a state name, the right-hand side is walked
symbol by symbol, and every identifier must resolve to a state name, a
parameter name, one of the reserved constants `t` or `pi`, one of the
delay variables (DDE), or one of the spatial operators `d1x`, `d2x`,
`d1y`, `d2y`, `lap` (PDE and PDE2D). Supported operators and functions
are listed in the table below. Any unknown symbol aborts the
specification with an explicit error citing the allowed set, so
misspelled parameters fail fast rather than producing silent zeros at
runtime.

| Category       | Accepted tokens                     |
|:---------------|:------------------------------------|
| Arithmetic     | \+ - \* /                           |
| Power          | ^ pow                               |
| Exponential    | exp log log2 log10 sqrt             |
| Trigonometric  | sin cos tan atan2                   |
| Inverse trig   | asin acos atan                      |
| Hyperbolic     | sinh cosh tanh                      |
| Absolute value | abs sign                            |
| Rounding       | floor ceiling                       |
| Comparison     | \< \> \<= \>= == !=                 |
| Conditional    | ifelse min max                      |
| Constants      | t pi                                |
| Spatial (PDE)  | d1x d2x d1y d2y lap                 |
| Delay (DDE)    | lag(state, tau) via delays argument |

Closed expression language accepted by the janos compiler. Every
identifier in a formula must resolve to a state, parameter, constant,
reserved name, or one of these tokens.

### The eleven model types

The `type` field of a `model_spec` is determined by inspecting which
arguments are supplied. The matrix below gives the diagnostic rules and
the downstream solvers.

| Type | Trigger argument | Default solver | Cache prefix |
|:---|:---|:---|:---|
| ode | rhs supplied | solver_rk4() | ode\_ |
| map | map supplied | solver_map() | map\_ |
| dde | rhs + delays | solver_dde() | dde\_ |
| sde | rhs + diffusion | solver_euler_maruyama() | sde\_ |
| jump_diffusion | rhs + diffusion (optional) + jumps | solver_jump_diffusion() | jd\_ |
| pdmp | regimes + switching | solver_pdmp() | pdmp\_ |
| ssa | stoichiometry + propensities | solver_ssa_direct() | ssa\_ |
| tau_leap | stoichiometry + propensities + hor | solver_tau_leap() | tau\_ |
| pde | pde (1D) + mesh | solver_mol() | pde\_ |
| pde2d | pde2d + mesh_x + mesh_y | solver_mol2d() | pde2d\_ |
| rdme | rdme (reactions + diffusivities + voxel structure) | solver_rdme() | rdme\_ |

Model type dispatch. The first argument combination that matches drives
type detection, validation, and the compilation template choice.

### Noise specifications

Stochastic models with exotic noise use `noise_spec` helpers that
compose with the base SDE type. The constructors are
[`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md)
for cross-correlated Gaussian noise with a Cholesky factor of the
covariance matrix,
[`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md)
for alpha-stable processes sampled through the Chambers-Mallows-Stuck
algorithm with optional Cholesky rotation,
[`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md)
for fractional Brownian motion generated by Wood-Chan circulant
embedding with a Hosking fallback when embedding fails, and
[`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md)
for 1/f^beta power-law noise synthesized in the frequency domain by FFT
filtering. Each noise spec attaches to a model through the `noise`
argument of
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md);
the compiler then routes to the corresponding SDE template
(`compiler_levy.R`, `compiler_fbm.R`, `compiler_colored.R`). Standalone
fBm and colored noise are restricted to single-replicate use because
their generation is non-Markovian; ensembles must use the R-level
parallel fallback.

``` r
# Deterministic Lotka-Volterra
lv <- model_spec(
    rhs = list(
        prey     ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init  = c(prey = 40, predator = 9)
)

# Same model with multiplicative Levy alpha-stable noise
lv_levy <- model_spec(
    rhs = list(
        prey     ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    diffusion = list(
        prey     ~ 0.1 * prey,
        predator ~ 0.1 * predator
    ),
    noise = levy_noise(alpha = 1.7, beta = 0),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init  = c(prey = 40, predator = 9)
)
```

  

## Compilation pipeline

Every compiled model passes through the same three-stage pipeline. The
**walker** (`expr_to_cpp()` in `R/compiler.R`) performs a depth-first
traversal of the formula’s abstract syntax tree and emits C++ tokens,
remapping R infix operators to their C equivalents, translating `^` into
`std::pow`, and substituting state and parameter references with indexed
array accesses. The **template renderer** selects a per-type C++
skeleton (one of thirteen in `R/compiler_*.R`) and splices the generated
expressions into the skeleton through `__HASH__` placeholders replaced
by `gsub(fixed = TRUE)`. The **runtime compiler** calls
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) with
the synthesized source, captures the resulting function pointers into an
environment, and attaches them to the `model_spec` object via
`compiled_id`. The entire pipeline is wrapped in a disk cache keyed by
an MD5 of the source: if the same formulas have been seen before, the
cached shared library is loaded without recompilation.

Compilation pipeline. Formulas enter at the top, pass through the AST
walker and a per-type template, and exit as cached function pointers in
the model_spec object.

### Symbolic Jacobian generation

Stiff integrators, bifurcation continuation, and several Lyapunov-family
builders need the Jacobian matrix J(y) = partial f / partial y of the
drift. Rather than asking the user for an analytic derivative or
resorting to finite differences, `compiler_jacobian.R` performs symbolic
differentiation of the parsed formula tree. The differentiation rules
cover the full closed grammar listed earlier, including `abs`, `sign`,
`ifelse`, the trigonometric and hyperbolic families, and the `pow`
operator with both integer and symbolic exponents. The output is a
matrix of C++ expressions that feeds directly into the Rosenbrock
template (`compiler_implicit.R`) for implicit time stepping, into the
continuation engine for fold and Hopf detection, and into the quadratic
Lyapunov builder when an equilibrium-point linearization is requested.
Because differentiation is symbolic, the result is exact up to floating
point round-off; there is no finite-difference step size to choose.

### Adjoint sensitivity

For parameter gradients of observables integrated along an ODE
trajectory, `compiler_adjoint.R` generates a continuous adjoint system
following the Cao-Li-Petzold (2003) formulation. The adjoint ODE is
integrated backward in time with the same RK4 scheme used for the
forward pass; the parameter gradient is accumulated through an inner
product against the pre-generated state trajectory. The generator emits
two extra C++ functions alongside the forward drift: one evaluates the
Jacobian-vector product J^T lambda, the other evaluates the
parameter-sensitivity product (partial f / partial p)^T lambda. This
route is materially faster than forward sensitivity analysis when the
number of parameters exceeds the number of observables, which is the
common case in calibration.

### Cache semantics

The cache is strictly content-addressed: the key is an MD5 of the
generated C++ source, not of the formula text. Two models that differ in
whitespace or in the order of algebraically equivalent terms but that
compile to identical source share a cache entry. The cache directory is
resolved through `tools::R_user_dir("janos", "cache")` so it persists
across sessions and respects platform conventions (XDG on Linux,
`~/Library/Caches` on macOS, `%LOCALAPPDATA%` on Windows). A cache entry
is a `.so`, `.dylib`, or `.dll` plus an index record. Invalidation is
automatic: changing any formula, parameter name, or reserved constant
flips the hash and triggers recompilation.

  

## Solver catalogue

Seventeen solver constructors cover deterministic, stochastic, spatial,
event-driven, and stiff regimes. Each solver is a lightweight
configuration object of class `solver_spec` carrying the integration
method name, the step size or tolerances, subsampling and
output-thinning settings, and solver-specific flags. The dispatcher
reads the method string and routes the compiled model through the
matching C++ runtime template.

| Constructor | Regime | Integration method | Order | Cost per trajectory |
|:---|:---|:---|:---|:---|
| solver_rk4() | ODE | Fixed-step RK4 | 4 | O(N_t \* d) |
| solver_rk45() | ODE | Adaptive Dormand-Prince 4(5) | 4(5) | O(N_t \* d) with step rejection |
| solver_rosenbrock() | ODE (stiff) | Rosenbrock Rodas3 + symbolic Jacobian | 3 | O(N_t \* d^3) dense / O(N_t \* d \* bw) banded |
| solver_map() | Map | Direct iteration | – | O(N_t \* d) |
| solver_dde() | DDE | RK4 + Hermite delay buffer | 4 | O(N_t \* d \* N_lag) |
| solver_euler_maruyama() | SDE | Euler-Maruyama | 0.5 (strong) / 1 (weak) | O(N_t \* d) |
| solver_milstein() | SDE | Milstein with central-FD g’(X) | 1 (strong) | O(N_t \* d^2) |
| solver_jump_diffusion() | SDE with jumps | Euler-Maruyama + Poisson jump channels | 0.5 / 1 | O(N_t \* d + N_jumps) |
| solver_pdmp() | PDMP | Lewis-Shedler thinning + per-regime drift | exact | O(N_t \* d + N_switches) |
| solver_ssa_direct() | CTMC | Gillespie direct | exact | O(N_events \* M) |
| solver_ssa_nrm() | CTMC | Gibson-Bruck next-reaction | exact | O(N_events \* log M) |
| solver_ssa_mnrm() | CTMC | Anderson modified next-reaction | exact | O(N_events \* log M) |
| solver_tau_leap() | CTMC approx. | Adaptive tau with step rejection | 1 (weak) | O(N_leaps \* M) |
| solver_tau_leap_midpoint() | CTMC approx. | Midpoint tau-leap | 1 (weak) | O(N_leaps \* M) |
| solver_hybrid() | CTMC / CLE | SSA / CLE switch by HOR threshold | exact / 1 | variable |
| solver_mol() | 1D PDE | Method of lines + RK4 | 4 (time) / 2 (space) | O(N_t \* N_mesh \* d) |
| solver_mol2d() | 2D PDE | Method of lines (tensor grid) + RK4 | 4 (time) / 2 (space) | O(N_t \* Nx \* Ny \* d) |
| solver_rdme() | Spatial CTMC | Gillespie over voxels or graph nodes | exact | O(N_events \* V) |

Solver catalogue. d = state dimension, N_t = number of time steps, M =
number of reaction channels, V = number of voxels or graph nodes, bw =
Jacobian bandwidth. For SSA, N_events is the number of firings drawn
adaptively from the propensity distribution.

### Decision tree for solver choice

The choice of solver is driven by the model type, by stiffness
considerations, and by whether the quantity of interest is a single
trajectory, an ensemble statistic, or an invariant distribution. The
flowchart below is the triage used internally by the documentation and
by the Lyapunov advisor; it also serves as a cheat sheet for users.

Solver selection flowchart. Diamond nodes are decision points,
rectangular nodes are recommended constructors. Stiffness is a dominant
driver for ODEs; for CTMCs the choice is between exactness and
computational cost.

### Runtime internals

The runtime layer lives in the per-type compiler files and in the fixed
stubs of `src/init.cpp`. All C++ compilation is deferred to runtime
through
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html); the
static shared library shipped with the package contains only the
Bartels-Stewart routine used by the Lyapunov subsystem
(`src/lyapunov_bartels_stewart.cpp`) and a minimal `R_init_janos()`
registration stub. Every generated C++ file uses `arma::vec` for state
vectors, `arma::ivec` for integer populations in SSA, `arma::mat` and
`arma::cube` for trajectory outputs, and the thread-safe xoshiro256+
PRNG bundled in the SSA and SDE templates when OpenMP is available.
State layout for spatial problems is contiguous: 1D PDEs use
`y[k * N + i]` for species k at voxel i, 2D PDEs use
`y[k * Nx * Ny + j * Nx + i]`, and RDME uses
`arma::imat state(n_species, n_voxels)` with `arma::icube` trajectory
output.

  

## Simulation engine

[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
is the single simulation entry point. It accepts a model, a solver, and
optional overrides for `t_max`, `init`, `parms`, and
`discard_transient`, and returns a `dyn_sim` S3 object. The object
carries the full trajectory, the attractor segment with transients
removed, a pointer to the source model, the effective solver
configuration, the elapsed wall-clock time, and any solver-specific
diagnostics (step rejection rate for RK45, leap rejection rate for
adaptive tau, jump counts for jump-diffusion, regime-switch counts for
PDMP).

### The dispatch order

The dispatch inside
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
is deterministic. The function first validates the model and the solver,
then checks whether the solver method is compatible with the model type
(an SSA solver on an ODE model, for example, fails with an explicit
error before any integration is attempted). It then reads `compiled_id`
and, when present, loads the cached shared object and binds the native
function pointers. If no `compiled_id` is available (rare, only in
legacy user-supplied `rhs` closures), the dispatcher falls back to
interpreted R evaluation with a non-trivial per-step cost that is
signaled to the user. The integration is then performed by the matching
runtime template; output is trimmed by `discard_transient`; an S3 object
is assembled with `class = "dyn_sim"`.

### The dyn_sim object

| Slot | Type | Purpose |
|:---|:---|:---|
| trajectory | matrix \[N_out x d\] | Full sampled trajectory including transient |
| attractor | matrix \[N_kept x d\] | Trajectory after discard_transient |
| t_grid | numeric \[N_out\] | Time grid corresponding to trajectory rows |
| model | model_spec | Back-pointer to the model spec that produced it |
| solver | solver_spec | Effective solver configuration (after defaults merged) |
| init_used | numeric \[d\] | Initial conditions actually used (after overrides) |
| parms_used | list | Parameters actually used (after overrides) |
| elapsed | numeric scalar | Wall-clock time in seconds |
| diagnostics | list or NULL | Per-solver diagnostics: step/leap rejections, jumps, switches |
| call | matched call | The call that generated the object (for reproducibility) |
| seed | integer scalar or NULL | Master seed, if any, for stochastic runs |
| interrupted | logical | Whether the run was interrupted and returned partial output |

Anatomy of a dyn_sim S3 object. Every analysis tool consumes a subset of
these slots; none of them are considered internal.

### Ensemble simulation

[`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md)
runs many independent replicates of the same model. Two backends
coexist. The **OpenMP batch backend** compiles a specialized C++
template that holds the entire replicate loop inside a single parallel
region, using per-thread xoshiro256+ state and per-thread state arrays.
It is available when the solver is
[`solver_ssa_direct()`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_euler_maruyama()`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
or
[`solver_tau_leap()`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md),
and when `vary` is NULL (uniform replicates). The **R-level fallback**
uses [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
on Linux and macOS or
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
across platforms. It covers every solver type and supports per-replicate
parameter variation via the `vary` argument. Automatic backend selection
picks the batch path when both conditions are met; otherwise the R
fallback is used.

Reproducibility is independent of the backend and of the number of
threads. Replicate i consumes seed = master_seed + i - 1, so the
sequence of outputs depends only on the master seed and the replicate
index; rerunning with a different thread count returns identical
trajectories. The output `ensemble_sim` object carries a 3D array of
trajectories, a matrix of terminal states, the running mean and standard
deviation over replicates, the number of extinct replicates for CTMC
models, and a diagnostics list.

Ensemble dispatch. The batch path runs all replicates inside a single
OpenMP parallel region; the fallback path spawns R workers through
mclapply or future.apply. Both paths share the seed protocol that
guarantees reproducibility.

  

## Observation layer

Real data is almost never the true state.
[`observe()`](https://robustecologies.github.io/janos/reference/observe.md)
applies a post-processing observation error model to the attractor
portion of a completed `dyn_sim`, producing an `observed_dyn_sim` object
that retains the noise-free trajectory alongside the noisy observation.
Keeping observation separate from integration is a deliberate design
choice: the same simulation can be observed under several noise regimes
without re-integrating, and the signal-to-noise diagnostics reported by
[`summary.observed_dyn_sim()`](https://robustecologies.github.io/janos/reference/summary.observed_dyn_sim.md)
compare the two copies directly. Six noise families are supported:
Gaussian, lognormal, Poisson, negative binomial, Student-t, and
multiplicative Beta for bounded states. Each family has a documented
parameterization so that users report identifiable quantities rather
than internal reparametrizations.

  

## Analysis toolkit

The analysis layer provides nineteen tools that extract structure from
models or from trajectories. The common pattern is that each tool
returns an S3 object with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods, uses
the model’s declarative structure directly when possible (never
differentiating finite-difference approximations when the symbolic
Jacobian is available), and produces publication-ready visualizations
with technical subtitles and captions.

### Phase portraits

A phase portrait draws the drift vector field on a two-dimensional slice
of state space together with a selection of trajectories. janos has five
phase-portrait constructors, each tailored to a model type:
[`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md)
for autonomous ODEs,
[`map_portrait()`](https://robustecologies.github.io/janos/reference/map_portrait.md)
for discrete maps with cobweb or basin-of-attraction visualizations,
[`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md)
which overlays the deterministic drift on Lyapunov ellipses computed
from the linearization and a fan of sample paths,
[`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md)
which combines the frozen-system flow with the Chebyshev spectrum of the
delay operator, and
[`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md)
which shows one flow per regime and traces the stochastic switching
history. Shared helpers live in `analysis_portrait_utils.R` so that grid
evaluation, nullcline detection, and ggplot theming are consistent
across variants.

### Fokker-Planck tools

For SDEs, the stationary probability density is computed through one of
three routes.
[`fp_potential()`](https://robustecologies.github.io/janos/reference/fp_potential.md)
reconstructs the one-dimensional quasi-potential from the drift and
diffusion fields,
[`fp_stationary()`](https://robustecologies.github.io/janos/reference/fp_stationary.md)
and
[`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md)
solve the stationary Fokker-Planck equation through finite-volume
discretization and sparse eigenvalue extraction, and
[`fp_kramers_rate()`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md)
computes Kramers escape rates between basins using the potential
curvature at the saddle and at the wells.
[`density_landscape_2d()`](https://robustecologies.github.io/janos/reference/density_landscape_2d.md)
produces a dimensional-reduction projection of a high-dimensional
stationary distribution onto a chosen plane. All Fokker-Planck tools
respect the noise specification: correlated, Levy, fBm, and colored
noise receive the appropriate generator.

### Bifurcation analysis

[`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
implements pseudo-arclength continuation (Keller, 1977) over a single
parameter, with fold detection through the determinant of the Jacobian
and Hopf detection through eigenvalue crossings of the imaginary axis.
The symbolic Jacobian generator feeds the continuation engine directly,
so bifurcation points are located to the precision of the underlying
linear algebra rather than to the resolution of a finite-difference
scheme.
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
is the higher-level wrapper that performs a multi-parameter sweep with
automatic refinement near detected codimension-one points.

### Quasi-stationary distribution and rare events

For CTMC models with absorbing states (extinction being the motivating
example), the dynamics conditional on non-absorption converge to a
quasi-stationary distribution.
[`estimate_qsd()`](https://robustecologies.github.io/janos/reference/estimate_qsd.md)
reports this distribution through iterated conditioning, with an
effective-sample-size diagnostic on the reweighting step.
[`estimate_extinction()`](https://robustecologies.github.io/janos/reference/estimate_extinction.md)
estimates mean absorption times and survival curves.
[`mlmc_estimate()`](https://robustecologies.github.io/janos/reference/mlmc_estimate.md)
provides multi-level Monte Carlo variance reduction for expensive
stochastic observables, coupling coarse and fine tau-leap trajectories
through the `compiler_mlmc.R` template to cancel most of the stochastic
variance analytically.

### Spectral and sensitivity analysis

[`spectral_analysis()`](https://robustecologies.github.io/janos/reference/spectral_analysis.md)
computes the Koopman or Perron-Frobenius spectrum of the deterministic
flow on a user-defined dictionary of observables using extended dynamic
mode decomposition, returning eigenvalues, eigenfunctions, and
reconstruction residuals.
[`adjoint_sensitivity()`](https://robustecologies.github.io/janos/reference/adjoint_sensitivity.md)
computes gradients of scalar observables with respect to parameters
along an ODE trajectory through the continuous adjoint, using the
adjoint C++ template described earlier. Together these tools provide the
building blocks for parameter calibration, model reduction, and
structural inference.

### Chaos diagnostics

For systems suspected of chaotic behaviour, janos ships a five-function
suite that reports quantitative fingerprints of a trajectory or of the
flow itself.
[`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
computes the full spectrum by propagating an orthonormal tangent basis
alongside the reference trajectory and restoring orthogonality through
periodic QR decomposition (Benettin-Shimada algorithm); it returns
running estimates of each exponent, the Kaplan-Yorke dimension, and
convergence plots.
[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
estimates the Grassberger-Procaccia correlation dimension from pair
statistics of attractor points with Theiler-window correction.
[`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
produces interpolated intersections of a trajectory with a user-chosen
hyperplane, collapsing a three-dimensional flow to a two-dimensional
return map.
[`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
implements the Gottwald-Melbourne 0-1 test, returning a median
correlation that separates regular from chaotic dynamics without
embedding; the stride is selected automatically from the autocorrelation
decay so that oversampled flows are not misdiagnosed as regular.
[`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
complements
[`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
by sampling the attractor (local maxima for flows, direct iterates for
maps) across a scanned parameter, revealing Feigenbaum cascades and
periodic windows on the true attractor rather than on equilibrium
branches; the scan is embarrassingly parallel through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) and
can optionally overlay the leading Lyapunov exponent at every scan point
on a second stacked panel.

The chaos toolkit shares the common S3 interface: every function returns
an object with [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) (where applicable),
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
For DDEs the dedicated `plot.dyn_sim(type = "delay_embedding")`
reconstructs a three-dimensional Takens embedding of a scalar
observable, and `plot.dyn_sim(type = "dygraph")` produces an interactive
time-series viewer through `dygraphs`.

| Tool | Input type | Output class | Dominant cost |
|:---|:---|:---|:---|
| phase_portrait() | ode model | phase_portrait | O(G \* d) |
| map_portrait() | map model | map_portrait | O(G \* d) |
| sde_portrait() | sde model | sde_portrait | O(G \* d + N_paths \* N_t \* d) |
| dde_portrait() | dde model | dde_portrait | O(G \* d + d^3) |
| pdmp_portrait() | pdmp model | pdmp_portrait | O(G \* R \* d) |
| fp_potential() | 1D sde model | fp_potential | O(G) |
| fp_stationary() | 1D sde model | fp_stationary | O(G^d_min) |
| fp_stationary_2d() | 2D sde model | fp_stationary_2d | O(Nx \* Ny) |
| fp_kramers_rate() | 1D sde model with bistable potential | fp_kramers | O(1) after fp_potential |
| density_landscape_2d() | dyn_sim of SDE or ensemble | density_landscape_2d | O(N_rep \* d \* N_t) |
| continuation() | model + parameter handle | continuation_result | O(N_steps \* d^3) |
| bifurcation_sweep() | model + two parameter handles | bifurcation_sweep | O(N_par^2 \* d^3) |
| estimate_qsd() | CTMC model | qsd_estimate | O(N_iter \* N_events \* M) |
| estimate_extinction() | CTMC model | rare_event_estimate | O(N_rep \* N_events \* M) |
| mlmc_estimate() | CTMC model (coupled tau-leap) | mlmc_estimate | O(N_levels \* N_rep \* M) |
| spectral_analysis() | ode or map with observables | spectral_analysis | O(N_obs^3) |
| adjoint_sensitivity() | ode model + observable | sensitivity_result | O(N_t \* d^2) |
| lyapunov_spectrum() | ode or map model | lyapunov_spectrum | O(N_renorm \* d^3) |
| correlation_dimension() | dyn_sim of chaotic flow | correlation_dimension | O(N_pts^2) |
| poincare_section() | dyn_sim of 3D flow | poincare_section | O(N_t) |
| zero_one_test() | dyn_sim trajectory | zero_one_test | O(n_c \* N_t) |
| bifurcation_diagram() | ode or map + parameter handle | bifurcation_diagram | O(N_par \* N_t) |

Analysis toolkit. G = number of grid points for spatial discretizations,
R = number of PDMP regimes, N_par = number of bifurcation parameters,
N_levels = number of MLMC levels, N_obs = size of the observable
dictionary, N_renorm = number of QR renormalisation intervals in the
Lyapunov spectrum, n_c = number of random frequencies in the 0-1 test.

  

## Lyapunov subsystem

Beyond simulation and analysis, janos carries a dedicated framework for
constructing and verifying Lyapunov functions. A Lyapunov function V(x)
is a scalar function that is zero at the equilibrium, strictly positive
elsewhere, and decreasing along trajectories; its existence certifies
local or global asymptotic stability of the equilibrium. The subsystem
is built around a single entry point,
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
which reads the model type and delegates to one of fourteen family
builders. A companion advisor,
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
inspects the model and recommends which families are theoretically
applicable before any construction is attempted.

### Family catalogue

| Family | Intended regime | Source file | Mathematical backend |
|:---|:---|:---|:---|
| quadratic | Linear / local ODE | lyapunov_quadratic | Bartels-Stewart (src/) |
| goh | Lotka-Volterra ODE | lyapunov_goh | Algebraic conservation laws |
| macarthur | MacArthur consumer-resource | lyapunov_macarthur | Algebraic |
| gilpin | Gilpin generalization | lyapunov_gilpin | Algebraic |
| sos | Polynomial ODE | lyapunov_sos | SOS programming (CVXR) |
| rbf | Semi-algebraic ODE | lyapunov_rbf | RBF + quadratic programming (quadprog) |
| massera | Nonlinear ODE (constructive) | lyapunov_massera | Massera series |
| cpa | Continuous PWA ODE | lyapunov_cpa | Continuous PWA via geometry |
| discrete | Discrete map | lyapunov_discrete | Discrete-time Lyapunov equation |
| stochastic | SDE | lyapunov_stochastic | Ito formula + PDE residual |
| krasovskii | DDE | lyapunov_krasovskii | Razumikhin / Krasovskii |
| pdmp | PDMP | lyapunov_pdmp | Per-regime quadratic + coupling |
| foster | Discrete CTMC (drift condition) | lyapunov_foster | Foster-Lyapunov drift test |
| functional | Delay functional | lyapunov_functional | Delay-functional construction |

Lyapunov family catalogue. Each family targets a specific mathematical
regime; the advisor selects among them based on the model type and
symbolic structure.

### The advisor

[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md)
walks the model, extracts a handful of structural features (type,
dimensionality, presence of a quadratic Jacobian at the equilibrium,
polynomial degree of the drift, conservation laws detectable through
null-space analysis, presence of stochastic or delay components, number
of PDMP regimes), and scores each family on a fitness metric combining
theoretical applicability and expected computational cost. The
top-ranked family is the default dispatch target. The advisor output is
a `lyapunov_advisor` S3 object that can be printed, inspected, and
plotted; users who want to override the automatic choice pass
`method = "quadratic"` (or any other family name) directly to
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md).

Lyapunov subsystem dispatch. The advisor ranks applicable families;
analyse_lyapunov() then calls the corresponding builder and wraps the
result in a lyapunov_report with a plot method that chooses the
informative visualization for the selected family.

### The Bartels-Stewart routine

The quadratic family (`lyapunov_quadratic`) requires the solution of the
continuous Lyapunov equation A^T P + P A = -Q for the symmetric
positive-definite P. This is the only C++ code shipped in the static
library; `src/lyapunov_bartels_stewart.cpp` implements the
Bartels-Stewart algorithm (1972) by real Schur decomposition of A,
substitution through the triangular quasi-factors, and reassembly. The
routine handles stable and unstable spectra correctly, returning an
error flag when the equation has no positive-definite solution. Because
Bartels-Stewart runs in O(d^3) and is exact up to floating-point
round-off, quadratic Lyapunov analysis scales to state dimensions well
beyond what symbolic approaches can handle.

### Certificate types

Each builder returns a certificate. For the quadratic, Goh, MacArthur,
and Gilpin families the certificate is algebraic: the Lyapunov function
is expressed in closed form and the decrease condition is verified
symbolically. For the SOS family the certificate is the SDP feasibility
witness returned by CVXR. For the RBF family the certificate is a set of
Gram-matrix inequalities verified by a QP. For the Massera construction
the certificate is convergence of the series truncation. For PDMP and
stochastic families the certificate combines algebraic decrease of the
deterministic part with numerical bounds on the stochastic generator.
The `lyapunov_report` S3 wrapper unifies these into a uniform interface:
[`print()`](https://rdrr.io/r/base/print.html) reports the family and
the decrease bound, [`summary()`](https://rdrr.io/r/base/summary.html)
gives the full verification trace, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows level
sets, decrease heatmaps, or per-regime sublevel sets depending on the
family.

  

## S3 class system

janos exports thirty-one S3 classes covering models, simulations,
analyses, Lyapunov reports, and chaos diagnostics. Every class
implements at least the `print` / `plot` pair; most also carry a
`summary` method. [`print()`](https://rdrr.io/r/base/print.html) is a
one-screen synopsis returning the object invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) extends the synopsis
with diagnostics and returns a dedicated summary class with its own
print method. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
defaults to the most informative single visualization for the class and
accepts a `type` argument for alternatives. Every plot method enforces a
uniform visual grammar: a technical subtitle naming the quantity being
displayed (for example, “quasi-stationary distribution, conditioning on
non-absorption”) and a grey caption citing the generating function and
any numerical conventions.

| Class                 | print | summary | plot |
|:----------------------|:-----:|:-------:|:----:|
| model_spec            |  yes  |   yes   |  no  |
| noise_spec            |  yes  |   no    |  no  |
| solver_spec           |  yes  |   no    |  no  |
| dyn_sim               |  yes  |   yes   | yes  |
| ensemble_sim          |  yes  |   yes   | yes  |
| observed_dyn_sim      |  yes  |   yes   | yes  |
| dynamical_system      |  yes  |   yes   |  no  |
| equilibrium_point     |  yes  |   yes   |  no  |
| phase_portrait        |  yes  |   yes   | yes  |
| map_portrait          |  yes  |   yes   | yes  |
| sde_portrait          |  yes  |   yes   | yes  |
| dde_portrait          |  yes  |   yes   | yes  |
| pdmp_portrait         |  yes  |   yes   | yes  |
| fp_potential          |  yes  |   yes   | yes  |
| fp_stationary         |  yes  |   yes   | yes  |
| fp_stationary_2d      |  yes  |   yes   | yes  |
| fp_kramers            |  yes  |   yes   | yes  |
| density_landscape_2d  |  yes  |   yes   | yes  |
| continuation_result   |  yes  |   yes   | yes  |
| bifurcation_sweep     |  yes  |   yes   | yes  |
| qsd_estimate          |  yes  |   yes   | yes  |
| rare_event_estimate   |  yes  |   yes   | yes  |
| mlmc_estimate         |  yes  |   yes   | yes  |
| spectral_analysis     |  yes  |   yes   | yes  |
| sensitivity_result    |  yes  |   yes   | yes  |
| lyapunov_function     |  yes  |   yes   | yes  |
| lyapunov_advisor      |  yes  |   yes   | yes  |
| lyapunov_report       |  yes  |   yes   | yes  |
| lyapunov_spectrum     |  yes  |   yes   | yes  |
| correlation_dimension |  yes  |   no    | yes  |
| poincare_section      |  yes  |   no    | yes  |
| zero_one_test         |  yes  |   no    | yes  |
| bifurcation_diagram   |  yes  |   no    | yes  |

S3 class method matrix. Thirty-one of thirty-three classes implement the
full print/summary/plot triad; specification classes omit plotting
because they are pre-simulation descriptors, and a few chaos diagnostics
return a scalar summary through print alone.

  

## Parallelism

janos exposes three levels of parallelism.

**C++ OpenMP inside batch solvers.** The ensemble batch templates for
SSA direct, Euler-Maruyama, and tau-leap execute the replicate loop
inside a single `#pragma omp parallel for`. Each thread owns a private
xoshiro256+ state initialized deterministically from the master seed so
that replicate assignments are stable across thread counts. Distribution
evaluations, propensity recomputation, and state updates all run without
R callbacks.

**R-level parallelism through mclapply or future.** The fallback path
spawns R workers for per-replicate
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
calls. It supports every solver type, including PDE and DDE which are
not in the OpenMP batch catalogue, and it supports per-replicate
parameter or initial-condition variation through the `vary` argument. On
Linux and macOS,
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) is
used when available; on Windows, a PSOCK cluster is used instead. When
`future.apply` is installed, users can override the backend through the
`future` plan.

**OpenMP inside per-step primitives.** A few analysis tools internally
parallelize over evaluation points:
[`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md)
vectorizes grid evaluation,
[`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md)
parallelizes the sparse matrix assembly over voxels, and the Lyapunov
RBF builder parallelizes kernel-matrix assembly over center pairs.

### Reproducibility

The seed protocol is universal. A single master seed reproduces every
stochastic run regardless of the backend, regardless of the number of
threads, and regardless of whether the ensemble is run in one call or
split across sessions. For replicates i = 1, …, N, the effective seed is
master + i - 1. This convention means that incrementally extending an
ensemble (rerunning with a larger N) preserves the first N trajectories
bit-for-bit.

### CRAN compliance

The parallel tests in `tests/testthat/` respect `_R_CHECK_LIMIT_CORES_`.
All parallel examples are wrapped in `\dontrun{}` when they exceed five
seconds; default test configurations set `n_cores = 1`. OpenMP is
guarded by `#ifdef _OPENMP` in every template, so compilation without
OpenMP (macOS without libomp, R CMD check on CRAN build farms) degrades
gracefully to single-threaded execution.

  

## Graceful interruption

Long simulations and stiff continuation runs benefit from cooperative
interruption. Every major engine in janos checks for user interrupts at
regular intervals: the forward integrators at every output step, the
bifurcation continuation at every corrector pass, the ensemble dispatch
between replicates, and the Lyapunov SOS solver at every CVXR iteration.
When Esc (or Ctrl+C) is detected the current run returns a partial
`dyn_sim` or analysis object with `interrupted = TRUE` set; the
[`print()`](https://rdrr.io/r/base/print.html) methods surface this flag
so that the user sees at a glance that the result is incomplete. For
multi-process ensembles the master process sends SIGINT to each worker,
waits briefly for graceful shutdown, and force-kills survivors;
collected partial replicates are returned rather than discarded.

  

## Complexity reference

For the common operating regimes, the dominant runtime and memory costs
are summarized below. These figures are measured empirically on the
included benchmarks and agree with the textbook analyses of the
underlying algorithms.

| Operation | Time complexity | Memory |
|:---|:---|:---|
| ODE integration, fixed step | O(N_t d) | O(N_out d) |
| ODE integration, adaptive step | O(N_t d) amortized, more with rejections | O(N_out d) |
| Implicit ODE integration | O(N_t d^3) dense | O(d^2) for Jacobian + O(N_out d) |
| DDE integration with k delays | O(N_t d k) | O(N_lag d) |
| 1D PDE, method of lines, M mesh points | O(N_t M d) | O(M d) + O(N_out M d) |
| 2D PDE, tensor mesh (Nx, Ny) | O(N_t Nx Ny d) | O(Nx Ny d) |
| SSA direct | O(N_events M) | O(M + d + N_out d) |
| SSA Gibson-Bruck NRM | O(N_events log M) | O(M) heap + O(N_out d) |
| Tau-leap, adaptive | O(N_leaps M) | O(M + N_out d) |
| Euler-Maruyama SDE | O(N_t d) with noise generation O(d) | O(N_out d) |
| Milstein SDE with multiplicative noise | O(N_t d^2) | O(d^2) + O(N_out d) |
| Jump-diffusion with N_J channels | O(N_t d + N_events N_J) | O(N_out d + N_J) |
| PDMP with R regimes | O(N_t d + N_sw R) | O(R d) + O(N_out d) |
| RDME grid, V voxels, M reactions | O(N_events (V + M)) | O(V M) + O(N_out V M) |
| RDME graph, V nodes + E edges | O(N_events (V + E)) | O(V + E) + O(N_out V) |
| Ensemble batch, R replicates, T threads | O(R / T \* cost_single) | O(R d) steady, O(d) per thread |
| Quadratic Lyapunov (Bartels-Stewart) | O(d^3) | O(d^2) |
| SOS Lyapunov of degree p | O(d^(3p)) SDP dominated | O(d^(2p)) |
| Continuation step with corrector | O(d^3) per Newton iter | O(d^2) |
| Stationary FP 1D, G cells | O(G) direct, O(G^2) eig | O(G) |
| Stationary FP 2D, Nx by Ny cells | O(Nx Ny) assembly, sparse eig | O(Nx Ny) sparse |
| Symbolic Jacobian generation | O(\|formulas\| d) | O(\|formulas\| d) |
| Formula-to-C++ compilation | O(\|C++ source\|) | O(source) |
| Cache hit | O(1) modulo dynamic loader | O(1) |

Complexity reference. d = state dimension, N_t = number of integration
steps, N_out = number of stored output points, M = number of reaction
channels, V = voxels or nodes, E = edges, N_J = jump channels, R =
replicates or regimes, T = threads, G = FP grid cells, p = SOS
polynomial degree.

  

## Extending janos

Because the compiler, the solver dispatch, and the analysis layer all
consume the `model_spec` contract, adding new capability is a matter of
plugging into one of three integration points.

### Adding a solver

A new solver is a template in `R/compiler_*.R` plus a constructor in
`R/solvers.R`. The template must emit C++ that conforms to the runtime
ABI (named arguments, state layout, output array shape) expected by
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
for the model type it serves. The constructor validates user-level
arguments and attaches a `method` string that the dispatcher reads. No
modification of `dyn_sim.R` is required for new solvers targeting
existing model types; the dispatcher discovers them through the method
string.

### Adding a model type

A new model type needs a branch in
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
that validates its specific arguments and sets the `type` field, a
per-type template in `R/compiler_*.R`, a default solver, and a dispatch
branch in
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md).
Downstream analysis tools and the Lyapunov advisor discover the new type
automatically when they encounter a model with the new `type` field;
either they support it through existing code paths (as for closely
related types) or they report that the type is outside their domain.

### Adding an analysis tool

A new analysis tool is a file in `R/analysis_*.R` implementing the tool
function, the three S3 methods (`print`, `summary`, `plot`) for its
output class, and the appropriate cross-references in roxygen
`@seealso`. The shared grid-evaluation and plotting utilities in
`analysis_portrait_utils.R` and `lyapunov_plot_utils.R` should be reused
wherever possible; reimplementing them breaks the visual consistency
guarantee of the plot layer.

### Adding a Lyapunov family

A new Lyapunov family is a file in `R/lyapunov_*.R` containing the
family builder, a registration entry in
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md)
to make it visible to automatic dispatch, and an entry in the
family-to-plot routing inside `lyapunov_plot_family.R`. The family
builder must return a `lyapunov_function` S3 object with the standard
slots (`V_fun`, `V_decrease_fun`, `certificate`, `family`, `model`) so
that the report wrapper and the plot method work without modification.

Extension points. Green nodes indicate where new code attaches; orange
nodes indicate the central contracts that must be respected.

  

## Quick reference

| Function | Purpose | Output class |
|:---|:---|:---|
| model_spec() | Declare a dynamical system (11 types) | model_spec |
| noise_spec helpers | Exotic noise: correlated, Levy, fBm, colored | noise_spec |
| solver_rk4() / rk45() / rosenbrock() | Deterministic ODE integrators | solver_spec |
| solver_map() / dde() | Discrete map and delay integrator | solver_spec |
| solver_euler_maruyama() / milstein() | SDE integrators (Euler-Maruyama, Milstein) | solver_spec |
| solver_jump_diffusion() / pdmp() | Event-driven solvers (JD, PDMP) | solver_spec |
| solver_ssa_direct() / nrm() / mnrm() | Exact CTMC simulation | solver_spec |
| solver_tau_leap() / tau_leap_midpoint() / hybrid() | CTMC approximations (tau-leap, hybrid SSA/CLE) | solver_spec |
| solver_mol() / mol2d() / rdme() | Spatial solvers (MOL 1D/2D, RDME) | solver_spec |
| dyn_sim() | Run a single simulation | dyn_sim |
| ensemble_sim() | Run many replicates with two-tier parallelism | ensemble_sim |
| observe() | Overlay observation noise on a simulation | observed_dyn_sim |
| phase_portrait() / map_portrait() | Deterministic phase portraits | phase_portrait / map_portrait |
| sde_portrait() / dde_portrait() / pdmp_portrait() | Stochastic/delay/regime-switching portraits | sde_portrait / dde_portrait / pdmp_portrait |
| fp_potential() / fp_stationary() / fp_stationary_2d() | Fokker-Planck stationary densities | fp_stationary(\_2d) / fp_potential |
| fp_kramers_rate() / density_landscape_2d() | Escape rates, density landscapes | fp_kramers / density_landscape_2d |
| continuation() / bifurcation_sweep() | Bifurcation continuation and sweeps | continuation_result / bifurcation_sweep |
| estimate_qsd() / estimate_extinction() | Quasi-stationary and extinction analysis | qsd_estimate / rare_event_estimate |
| mlmc_estimate() / spectral_analysis() | Variance reduction and spectral decomposition | mlmc_estimate / spectral_analysis |
| adjoint_sensitivity() | Continuous adjoint parameter gradients | sensitivity_result |
| analyse_lyapunov() / lyapunov_advisor() | Lyapunov construction and advisor | lyapunov_report / lyapunov_advisor |
| lyapunov_function() + family builders | 14 family builders for stability certificates | lyapunov_function |
| lyapunov_spectrum() / correlation_dimension() | Lyapunov spectrum and correlation dimension | lyapunov_spectrum / correlation_dimension |
| poincare_section() / zero_one_test() | Poincare section and 0-1 test for chaos | poincare_section / zero_one_test |
| bifurcation_diagram() | Attractor-level bifurcation diagram | bifurcation_diagram |

Primary API reference.

### End-to-end example

``` r
library(janos)

# 1. Declare a bistable SDE
ms <- model_spec(
    rhs = list(x ~ x - x^3 + alpha),
    diffusion = list(x ~ sigma),
    state_names = "x",
    parms = list(alpha = 0.0, sigma = 0.35),
    init  = c(x = -1)
)

# 2. Integrate
sim <- dyn_sim(ms, t_max = 400,
               solver = solver_euler_maruyama(dt = 0.005),
               discard_transient = 50)

# 3. Inspect
plot(sim)
summary(sim)

# 4. Stationary Fokker-Planck density
fp <- fp_stationary(ms, grid = seq(-2.5, 2.5, length.out = 401))
plot(fp)

# 5. Kramers escape rate between the two wells
rate <- fp_kramers_rate(fp)
print(rate)

# 6. Ensemble with the OpenMP batch backend
ens <- ensemble_sim(ms, solver = solver_euler_maruyama(dt = 0.005),
                    t_max = 200, n_rep = 1024, seed = 1L)
summary(ens)

# 7. Lyapunov analysis (falls back to stochastic family)
lya <- analyse_lyapunov(ms)
plot(lya)

# 8. Add observation noise and compare
obs <- observe(sim, noise = "gaussian", sigma = 0.1)
plot(obs)
```

This eight-step pipeline touches every subsystem of the package:
declarative specification, runtime compilation, deterministic and
stochastic integration, analytical post-processing, ensemble
parallelism, stability analysis, and observation modeling. None of the
intermediate objects require format translation; each step consumes what
the previous step produced.
