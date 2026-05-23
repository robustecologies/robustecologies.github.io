# Changelog

## lucifer 3.10.0

### Breaking changes

  

The S3 classes assigned to the outputs of
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md),
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md)
and
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
are renamed with a `lucifer_` prefix to avoid colliding with the
homonymous classes exported by the CRAN package `loo`. The class mapping
is `loo` -\> `lucifer_loo`, `waic` -\> `lucifer_waic`, `kfold` -\>
`lucifer_kfold`, `lfo` -\> `lucifer_lfo`, `stacking_weights` -\>
`lucifer_stacking_weights`, `loo_comparison` -\>
`lucifer_loo_comparison`. The S3 methods `print`, `plot` and `summary`
follow the rename: `print.lucifer_loo`, `plot.lucifer_loo`,
`summary.lucifer_loo`, etc. The functions themselves
([`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
etc.) keep their names, signatures, return-value field structure and
computation; the change is invisible to code that calls these functions
and prints, plots or accesses fields of the result via the standard S3
generics. The change is visible only to code that does explicit class
introspection: any `inherits(x, "loo")`, `inherits(x, "waic")` and
similar must be updated to the `lucifer_*` class names. The CRAN package
`loo` is now declared in `Suggests` so that downstream code can use the
loo package’s tools alongside lucifer without the startup
`Registered S3 methods overwritten by 'lucifer'` warning previously
emitted on load.

  

## lucifer 3.9.0

### API surface reduction

  

Moved 120 dispatcher-internal sampler and optimiser implementations to
the package private namespace. These dot-prefix functions (`.la*`,
`.mcmc*`, `.vb*`, `.iq*`, `.colVars`, `.rowVars`,
`.run_one_pathfinder_path`) were never part of the documented API; they
are dispatched internally by
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
and
[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md).
The user-facing export count is now 446 rather than 566. Any code that
called these symbols directly must switch to the public dispatcher with
the corresponding `Algorithm` argument.

  

### Documentation completeness

  

All exported functions now carry full roxygen2 blocks compliant with the
project documentation standard: `@title`, `@description`, `@param`,
`@return`, `@details` with mathematical or algorithmic exposition,
`@references` (primary literature, verified DOIs), `@examples` wrapped
in `\dontrun{}`, and `@seealso` with bidirectional S3 family
cross-references. No user-visible behaviour is affected.

------------------------------------------------------------------------

## lucifer 3.8.0

### Bayesian mode inference

New
[`BayesMode()`](https://robustecologies.github.io/lucifer/reference/BayesMode.md)
function for fully Bayesian inference on the number and location of
modes in any distribution, implementing the Sparse Finite Mixture (SFM)
framework of Cross, Hoogerheide, Ulker and van Dijk (2024, *Economics
Letters*) extended to 11 distribution families. The method fits a
mixture model via Gibbs MCMC with a sparsity-inducing Dirichlet
hyperprior that automatically prunes unnecessary components, then
detects modes from the fitted mixture at each MCMC draw to produce a
posterior distribution over both mode count and mode locations. The
result is a `bayes_mode` S3 object with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

The C++ engine (`src/mode_sfm.cpp`, ~750 lines) implements
family-specific Gibbs update steps for all 11 distributions: Normal,
Skew-Normal, Student-t, Gamma, Log-Normal, Beta, Poisson,
Shifted-Poisson, Negative Binomial, Geometric, and Binomial. A companion
C++ file (`src/mode_detect.cpp`, ~400 lines) provides three mode
detection algorithms: fixed-point iteration for Normal mixtures
(Carreira-Perpinán 2000), discrete enumeration for count distributions,
and grid search with golden-section refinement for general continuous
families. Both files support OpenMP parallelization.

[`BayesMode()`](https://robustecologies.github.io/lucifer/reference/BayesMode.md)
accepts raw numeric vectors or any lucifer posterior object (`demonoid`,
`laplace`, `iterquad`, `vb`, `pmc`); when given a posterior object, it
returns a `bayes_mode_multi` container with per-parameter mode inference
and a compact one-line-per-parameter summary table.

[`print.bayes_mode()`](https://robustecologies.github.io/lucifer/reference/print.bayes_mode.md)
produces bayestestR-inspired pipe-separated tabular output showing the
posterior distribution over mode count, P(multimodal), and mode
locations with 95% credible intervals.
[`summary.bayes_mode()`](https://robustecologies.github.io/lucifer/reference/summary.bayes_mode.md)
adds MCMC diagnostics (ESS for alpha, effective K, log-likelihood).
[`plot.bayes_mode()`](https://robustecologies.github.io/lucifer/reference/plot.bayes_mode.md)
provides six plot types: `"modes"` (data density with posterior mode
bands), `"n_modes"` (bar chart of P(M modes)), `"mixture"` (fitted
mixture density with posterior draws), `"trace"` (MCMC diagnostics),
`"posterior"` (violin plots of mode locations), and `"convergence"`
(running mean diagnostics).

The family auto-detection logic selects Poisson for non-negative integer
data, Gamma for positive continuous data, and Normal for general
continuous data. Users can override with `family = "shifted_poisson"`
(or any of the 11 families) and tune the sparsity prior via
`alpha.prior`.

### Mode estimation overhaul

[`Mode()`](https://robustecologies.github.io/lucifer/reference/Mode.md)
and
[`Modes()`](https://robustecologies.github.io/lucifer/reference/Mode.md)
now return objects of class `mode_estimate` with dedicated
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
return structure is fully backward-compatible: `Mode(x)[1]`,
`length(Modes(x)[[1]])`, and `Modes(x)$modes` all work identically to
previous versions.

Both functions gain a `method` argument supporting four estimation
strategies: `"kde"` (kernel density, the existing default), `"hsm"`
(half-sample mode; Bickel & Fruhwirth 2006), `"venter"` (Venter’s mode,
shortest interval containing half the data), and `"shorth"` (shortest
half-range). The HSM and Venter estimators have C++ backends for
performance.
[`summary.mode_estimate()`](https://robustecologies.github.io/lucifer/reference/summary.mode_estimate.md)
with `compare = TRUE` tabulates mode estimates across all four methods.

[`plot.mode_estimate()`](https://robustecologies.github.io/lucifer/reference/plot.mode_estimate.md)
produces ggplot2 visualizations with
[`theme_relab()`](https://robustecologies.github.io/lucifer/reference/theme_relab.md):
`"density"` (default) overlays vertical mode markers on a kernel density
curve, `"histogram"` marks modes on a histogram, and `"comparison"`
shows a panel comparing all four methods on the same data.

Added `.print_pipe_table()` reusable utility for bayestestR-style
pipe-separated tabular console output, used across all new print
methods.

------------------------------------------------------------------------

## lucifer 3.7.1

### Burn-in overhaul

[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
gains a `BurnIn` argument that gives the user full control over burn-in
(warmup) removal. When `NULL` (default), the BMK stationarity diagnostic
estimates the burn-in boundary automatically as before. When a
non-negative integer, that many thinned samples are discarded from the
start, bypassing BMK entirely. The returned object includes a
`BurnIn.Method` field (`"user"` or `"BMK"`) shown in
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html). The same parameter
threads through
[`Combine()`](https://robustecologies.github.io/lucifer/reference/Combine.md)
for multi-chain runs.

Fixed a critical degeneracy in `Posterior2` construction. When BMK
detected no burn-in (all samples stationary), `Posterior2` collapsed to
a single-row vector instead of equaling `Posterior1`. This caused
`Summary2` and `DIC2` to be `NA`,
[`predict.demonoid()`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md)
to fall back to `Posterior1` with a warning, and `LML` estimation to
silently fail. `Posterior2` is now always a proper matrix. When BMK
finds no stationary region at all, the first 50% of samples are
discarded as a conservative fallback (with a warning) instead of
producing a degenerate single-row result.

[`Combine()`](https://robustecologies.github.io/lucifer/reference/Combine.md)
gains a `BurnIn` argument matching
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and its stationarity assessment mirrors the same fix: proper Posterior2
matrix, 50% fallback when BMK fails.

### Interrupt handling and progress reporting

Added `check_interrupt_r()` calls inside the stepping-out and shrinkage
loops of `.slice_sample_1d()` (auto-FC Gibbs), fixing unresponsiveness
to Esc during slice sampling. Previously, pressing Esc had no effect
until the slice sample completed (up to 300 Model evaluations without an
interrupt check).

Fixed a crash (“subscript out of bounds”) when pressing Esc twice during
parallel chain result collection. The callr interrupt handler now
catches secondary interrupts via `interrupt = function(e) NULL` in the
inner `tryCatch`.

Added stall detection to the parallel chain progress poller: if a chain
shows no progress for 120 seconds, a warning is printed advising the
user to press Esc.

Replaced inline status printing with `print_status()` in six C++
samplers that were bypassing the progress file mechanism: DEMC, AIES,
twalk, DREAM (`sampler_ensemble.cpp`) and HARM branches 2/4
(`sampler_other.cpp`). These samplers now write progress to the callr
polling file and respond to interrupt checks, fixing the “0% -\> done”
progress display in parallel runs. Added `.write_progress()` to the
R-level flowMC and NeuTra samplers.

### SMC console output

[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
now prints professional per-stage status showing the temperature
schedule (beta), ESS with percentage, running log marginal likelihood
estimate, and flags for resampling and rejuvenation acceptance rate. A
header summarizes the configuration (particles, rejuvenation method,
schedule type) and a footer reports the total stage count and wall time.
IBIS mode prints analogous per-batch status. The previous implementation
only printed at `stage %% Status` intervals, which for adaptive
schedules with fewer than 100 stages meant no output was ever shown.

### Spike-and-slab (SSVS) visualization module

New
[`ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/ssvs_summary.md)
constructor and S3 class providing dedicated tooling for spike-and-slab
variable selection models. The constructor extracts posterior inclusion
probabilities (PIPs) from any lucifer fit (`demonoid`, `vb`, `pmc`,
`smc`, `abc`, `sbi`, `laplace`, `iterquad`, `bayesquad`), accepts
indicator and coefficient parameter patterns either by name (regex) or
column index, and optionally takes ground truth indicators and array
dimensions for matrix-valued indicators (e.g. VAR adjacency matrices).
[`detect_indicators()`](https://robustecologies.github.io/lucifer/reference/detect_indicators.md)
auto-discovers binary `{0, 1}` columns in a posterior matrix when
explicit patterns are not supplied.

[`print.ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/print.ssvs_summary.md)
reports the median-probability model with coloured selection symbols and
optional ground truth annotations.
[`summary.ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/summary.ssvs_summary.md)
produces tabular output including conditional coefficient means and
standard deviations given inclusion, plus Bayesian false discovery rate
(FDR) and false omission rate (FOR) at the chosen threshold when ground
truth is available.

[`plot.ssvs_summary()`](https://robustecologies.github.io/lucifer/reference/plot.ssvs_summary.md)
dispatches to seven specialized visualizations: a PIP bar chart with
prior and threshold reference lines (`type = "pip"`); spike vs slab
prior densities overlaid against the marginal posterior of each
coefficient (`type = "spike_slab"`); coefficient density split by
inclusion status (`type = "conditional"`); a matrix heatmap reshaped via
`dims` for VAR-style adjacency layouts (`type = "heatmap"`);
cumulative-PIP trajectories of the most uncertain indicators across MCMC
iterations (`type = "trajectory"`); a directed network graph where edge
opacity and width are proportional to PIP, with self- loops absorbed
into node markers (`type = "network"`); and Bayesian FDR/FOR curves as
functions of the inclusion threshold when ground truth is supplied
(`type = "fdr"`).

### Fitted PPC plot enhancements

The `Style = "Fitted"` posterior predictive plot in `.plot_ppc` and its
multivariate variants `"Fitted, Multivariate, C"` and
`"Fitted, Multivariate, R"` have been redesigned without breaking the
default look. Axis labels now use proper mathematical typography via
[`expression()`](https://rdrr.io/r/base/expression.html) and
[`bquote()`](https://rdrr.io/r/base/bquote.html) rendering and rather
than literal “y” / “yhat” strings. A dashed grey 1:1 identity line has
been added so departures from perfect calibration are immediately
visible, and the plotting region is forced equal via
[`coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
so the diagonal is geometrically meaningful. Points are now coloured by
the absolute residual via a `scale_color_gradient`, providing a visual
cue for outlier observations. A subtle grey caption in the bottom-right
summarizes the elements of the plot (predictive intervals, identity
line, LOESS) so the figure is interpretable without external annotation.
The original posterior-predictive interval bars and LOESS smoother are
preserved.

### Posterior imputation visualization

New
[`plot_imputed()`](https://robustecologies.github.io/lucifer/reference/plot_imputed.md)
function for publication-quality visualization of Bayesian missing-data
imputation in univariate and multivariate time series. The function
accepts any lucifer fit (`demonoid`, `vb`, `pmc`, `smc`, `abc`, `sbi`,
`laplace`, `iterquad`, `bayesquad`) together with the original data
containing `NA` values, locates the imputed-value parameters in the
posterior matrix (auto-detecting columns whose names match `imp` or
`miss`), matches them positionally to the `NA` cells in column-major
order, and renders the time series with observed values, posterior
median imputations, credible intervals, and (optionally) ground truth.

Three `Style` options are available. `"default"` draws vertical credible
interval bars and hollow diamonds at each missing position; `"ribbon"`
converts the per-position intervals into a continuous shaded band best
suited to structured or contiguous missingness; `"draws"` overlays a
spaghetti of posterior samples threading through the missing positions,
revealing the joint correlation structure of nearby imputations that
marginal credible intervals collapse away. The multivariate path facets
by series via `facet_wrap(~ series, scales = "free_y")`. A grey caption
in the bottom-right documents the visual encoding of the chosen `Style`,
and the subtitle reports the count of observed and imputed cells with
the missing fraction. Custom time indices (numeric, `Date`, or
`POSIXct`), explicit `indicators`/`positions` overrides, custom
`series_names` for facet labels, and a `level` argument for the credible
interval are all supported.

### Categorical posterior predictive summary fix

Fixed a critical bug in the `Categorical = TRUE` branch of the
`summary.*.ppc()` family (affecting `summary.demonoid.ppc`,
`summary.laplace.ppc`, `summary.pmc.ppc`, `summary.vb.ppc`, and
`summary.iterquad.ppc`). The lift and `p(yhat[i,] != y[i])` discrepancy
columns used `grep(Summ[i, 1], names(catcounts))` to map each
observation to its category, which is a regex match against the
stringified category labels. With Poisson-style count data the category
names are `"0"`, `"1"`, `"2"`, …, `"10"`, `"11"`, so for example
`grep("1", c("0","1","2","10","11"))` returns three positions (matching
`"1"`, `"10"`, and `"11"`). The multi-element vector this produced could
not be assigned into a single matrix cell, raising “número de elementos
para sustituir no es un múltiplo de la longitud del reemplazo” / “number
of items to replace is not a multiple of replacement length” and
breaking `plot(., Style = "Covariates, Categorical DV", ...)` for any
count outcome with multi-digit categories. Replaced with
`match(as.character(Summ[, 1]), names(catcounts))`, which performs the
correct exact lookup in a single vectorized call before the loop.

### Expanded `plot.demonoid.ppc` documentation

The roxygen block for
[`plot.demonoid.ppc()`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)
has been substantially expanded. The `Style` argument is now documented
in a dedicated `@details` section that organizes the 40 plot styles into
seven thematic groups (univariate predictive distributions, discrepancy
and goodness-of-fit, fitted-vs-observed, residual diagnostics,
predictive intervals and counts, covariate-aware checks, spatial and
time-series checks, calibration, and grouped diagnostics) and lists
exactly which styles require the `Data` argument. A dedicated subsection
on the `"Fitted"` family explains the 1:1 identity line, the
equal-aspect coordinate region, the LOESS smoother, the
residual-magnitude colour gradient, and the `"Multivariate, C"` /
`"Multivariate, R"` variants. The `...` argument now lists the most
useful pass-through options (`Group`, `stat_fun`, `stat_fun2`, `loo`,
`forecast_start`, `col`). Added a runnable `@examples` block
demonstrating `Style = "Density"`, `"Fitted"`, `"Density Overlay"`,
`"Intervals"`, `"Ribbon"`, `"Covariates"`, and the convenience wrappers
[`ppc_dens_overlay()`](https://robustecologies.github.io/lucifer/reference/ppc_dens_overlay.md)
and
[`ppc_stat()`](https://robustecologies.github.io/lucifer/reference/ppc_stat.md).
References to Gelman, Meng & Stern (1996) on PPC and Vehtari, Gelman &
Gabry (2017) on PSIS-LOO have been added, and the `@seealso` block now
cross-links `summary.demonoid.ppc`, `ppc_dens_overlay`, `ppc_stat`,
`ppc_rootogram`, and the new `plot_imputed`.

------------------------------------------------------------------------

## lucifer 3.7.0

### SDE infrastructure audit

Exact transition density likelihoods for Ornstein-Uhlenbeck, geometric
Brownian motion, and Cox-Ingersoll-Ross processes with observation noise
have been replaced with affine Kalman filter implementations
(`loglik_ou_kalman`, `loglik_gbm_kalman`, `loglik_cir_kalman` in
`src/sde_engines.cpp`). The previous code used the noisy observation
y\[k-1\] as a proxy for the previous latent state in the transition
density; with observation noise present this produces biased posteriors:
kappa overestimated by a factor of four and obs_sd underestimated. The
new implementations track the filtered state E\[X_k \| Y\_{1:k}\] via
proper Kalman recursion.
[`simulate.sde_model()`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md)
now adds observation noise by default (`add.obs.noise = TRUE`),
correcting posterior predictive coverage from approximately 91% to 95%.
The EKF covariance update uses the Joseph form P = (I-KH)P(I-KH)’ + KRK’
throughout; the initial covariance was changed from 1e-6·I to 1e-2·I and
the finite-difference step to sqrt(machine_epsilon) for optimal
central-difference accuracy. A new C++ function
`sde_filtered_states_exact_cpp()` returns Kalman-filtered state
estimates for exact families.
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
now warns when Euler-Maruyama or Milstein is selected with observation
noise present.

SDE trajectory and prediction fixes. The trajectory plot now uses
Kalman-filtered state estimates via `.filtered_posterior_paths()`,
achieving correlation with observed data of 0.97+ (previously near zero
due to unconditional forward simulation). `predict.sde_fit` had a
double-simulation bug that discarded the x0 override on the second call;
this is fixed. `sde_simulate_exact_cpp` now accepts an optional x0
initial state to start paths from the first observation rather than the
stationary distribution.

### Robust interrupt handling

Multi-chain parallel execution responds to Esc and Ctrl+C gracefully:
callr children receive SIGINT via `handle$interrupt()`, and PSOCK
workers check a stop-signal file written by C++
`set_stop_signal_cpp()`/`clear_stop_signal_cpp()`. After a grace period,
the parent collects whatever partial results are available rather than
crashing. R-level MCMC samplers (Gibbs, GG, AGG, LSS, QSS, GPSS, NeuTra,
flowMC, MCMCMC) cache partial state at thinning points via
`.cache_partial()`/`.get_cached_partial()` in the new
`R/interrupt_utils.R`. VB methods return the accumulated partial ELBO
trace on interrupt rather than a two-row zeros stub. PMC,
LaplaceApproximation, and IterativeQuadrature gain `tryCatch(interrupt)`
wrappers. A C++-level `maybe_check_interrupt()` function in
`sampler_common.h` uses `std::chrono::steady_clock` at 200 ms intervals
to check for user interrupts during computation-intensive inner loops.
All output objects gain an `Interrupted` field; all `print` methods
display a warning banner when results come from an interrupted run.

### Algorithm selection system expansion

`algo_registry.R` gains six new metadata fields per entry:
`multimodal_affinity`, `constraint_affinity`, `requires_torch`,
`eval_cost_multiplier`, `dim_range`, and `quality_tier`. Algorithm
subcategories have been corrected: Simulated Tempering, NRST,
Wang-Landau, and NRPT moved to the multimodal subcategory; BPS,
Boomerang, ZigZag, and Randomised HMC moved to pdmp; Projected Langevin
and ProxMCMC moved to constraint; NeuTra and flowMC moved to flow;
Delayed Acceptance moved to surrogate.
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
replaces seven hard-coded MCMC scoring blocks with a unified
`.prescribe_mcmc_factors()` dispatching on 13 subcategories.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
and
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
dispatch updated to cover all 82+ algorithms. New exported function
[`algo_info()`](https://robustecologies.github.io/lucifer/reference/algo_info.md)
supports registry introspection with filtering by category, subcategory,
and gradient requirement.

------------------------------------------------------------------------

## lucifer 3.6.0

### Multimodal samplers

Three sampling algorithms for multimodal posteriors added to
`src/sampler_gradient.cpp`. **Simulated Tempering** (Marinari & Parisi
1992): augmented state space with a temperature variable; Gibbs updates
alternate between the parameter draw at the current temperature and a
temperature transition proposal. **NRST** (Non-Reversible Simulated
Tempering; Sakai & Hukushima 2016): deterministic even-odd sweep over
temperature levels with direction variables, improving mixing relative
to reversible tempering. **Wang-Landau** (Wang & Landau 2001):
flat-histogram adaptive estimation of the density of states, converging
to uniform visits across the temperature ladder.

### Flow-enhanced MCMC

**NeuTra** (neural transport MCMC; Hoffman et al. 2019) reparameterises
the target via a trained normalising flow, running standard HMC in the
latent space where the geometry is approximately Gaussian. **flowMC**
(Gabrie et al. 2022) adapts the flow online using MCMC samples,
improving flow quality iteratively. `R/normalizing_flow.R` implements a
RealNVP module via torch `nn_module` for both methods. Both samplers are
implemented in pure R with `tryCatch(interrupt)` wrappers for partial
result recovery.

### Coupling diagnostics

[`coupling_diagnostic()`](https://robustecologies.github.io/lucifer/reference/coupling_diagnostic.md)
added to `R/coupled_mcmc.R` computes L-lag total variation upper bounds
following Biswas, Jacob and Vanetti (2019). The function accepts two
chains run from a common coupling kernel and produces TV bound traces
along with meeting time statistics. S3 `print` and `plot` methods
included.

### Constraint-handling samplers

**Projected Langevin**: Euler-Maruyama discretisation of the Langevin
diffusion with hard projection onto box or simplex constraints at each
step. **ProxMCMC** (Durmus et al. 2022): Moreau-Yosida envelope
smoothing for L1 and total variation penalties, enabling gradient-based
sampling of non-smooth posteriors. Both implemented in
`src/sampler_gradient.cpp`.

### Additional algorithms

**Relativistic MC** (Lu, Niu, Lan & Liu 2017): replaces the standard
kinetic energy with a relativistic form, producing heavier-tailed
momentum proposals that improve exploration of funnel-shaped posteriors.
**NR Overdamped Langevin** (Duncan, Lelièvre & Pavliotis 2016): adds a
non-reversible antisymmetric perturbation to the overdamped Langevin
diffusion, improving convergence rate without changing the invariant
distribution. **GIST-NUTS** (Bou-Rabee, Carpenter & Marsden 2024): at
each iteration, K candidate step sizes are drawn on a log-uniform grid;
the grid is scored via softmax(−α\|ΔH\|) and a step size is
Gibbs-sampled from the resulting distribution; a full NUTS tree is then
built at the selected ε. This design achieves approximately 87%
acceptance rate compared to approximately 70% for standard NUTS on test
models.

Total MCMC algorithm count after Phases 51–58 and the additional
proposals: 82.

------------------------------------------------------------------------

## lucifer 3.5.0

### PDMP samplers

Three piecewise deterministic Markov process (PDMP) samplers added to
`src/sampler_gradient.cpp`. **Bouncy Particle Sampler** (Bouchard-Côté,
Vollmer & Doucet 2018): velocity refreshment at a Poisson rate and
reflection at equipotential boundaries, producing a non-reversible
continuous-time trajectory. **Boomerang Sampler** (Bierkens et
al. 2020): replaces the global reference measure with a Gaussian
reference trajectory, achieving lower event rates for log-concave
targets. **Randomised HMC** (Bou-Rabee & Sanz-Serna 2017): stochastic
trajectory length drawn uniformly from \[0, 2πL\], breaking the
resonance artefacts of fixed-length HMC.

### Modern slice sampling

Three new slice samplers implemented in pure R. **QSS** (Quantile Slice
Sampler; Heiner et al. 2024): introduces a unit-uniform auxiliary
variable and constructs proposals via quantile evaluations, reducing the
number of likelihood evaluations compared to standard stepping-out.
**LSS** (Latent Slice Sampler; Li & Walker 2023): augments the target
with a latent variable whose marginal characterises the slice, enabling
efficient sampling in high dimensions. **GPSS** (Gibbsian Polar Slice
Sampler; Biron-Lattes, Surjanovic, Campbell & Bouchard-Côté 2023): uses
polar coordinates in the augmented space to avoid the difficulties of
anisotropic targets in the standard slice framework. All three include
`tryCatch(interrupt)` wrappers for partial result recovery.

------------------------------------------------------------------------

## lucifer 3.4.0

### Riemannian and geometric MCMC

Three samplers exploiting the local curvature of the posterior added in
`src/sampler_geometric.cpp` (~400 lines) with shared infrastructure in
`src/metric_tensor.h` and `src/metric_tensor.cpp`. **RMHMC** (Riemannian
Manifold HMC; Girolami & Calderhead 2011): uses a position-dependent
Fisher information metric via SoftAbs regularisation of the Hessian,
adapting the leapfrog geometry to local curvature without divergences.
**LMC** (Lagrangian Monte Carlo; Lan, Streets & Shahbaba 2015): casts
the dynamics as a Lagrangian system with explicit velocity variables and
a geometric integrator, producing a sampler well-suited to non-Euclidean
target geometries. **MHMC** (Magnetic HMC; Tripuraneni, Rowland,
Ghahramani & Turner 2017): adds an antisymmetric component to the
symplectic dynamics, breaking detailed balance locally while preserving
the target as the invariant distribution and improving mixing. Three
metric modes are supported: identity (standard HMC), Hessian (numerical
Fisher information), and SoftAbs (regularised absolute-value spectrum).

### Locally adaptive gradient methods

Three gradient-based algorithms with automatic step-size control
appended to `src/sampler_gradient.cpp`. **autoMALA** (Biron-Lattes,
Surjanovic, Campbell & Bouchard-Côté 2024): searches for a locally
optimal Metropolis-adjusted Langevin step size using a geometric binary
search, adapting to the local curvature at each iteration without a
global warm-up phase. **MALT** (Metropolis Adjusted Langevin
Trajectories; Riou-Durand & Vogrinc 2022): applies BAOAB operator
splitting to the underdamped Langevin dynamics, combining the mixing
speed of HMC with the momentum persistence of SGLD. **AAPS**
(Apogee-to-Apogee Path Sampler; Sherlock, Thiery & Vogrinc 2023):
constructs HMC trajectories from apogee to apogee of the Hamiltonian,
using multinomial trajectory sampling to select the next state.

------------------------------------------------------------------------

## lucifer 3.3.0

### Parallelisation overhaul

[`lucifer_threads()`](https://robustecologies.github.io/lucifer/reference/lucifer_threads.md)
provides a user-facing API for thread budget management with four
strategies: auto (detect available cores and distribute optimally),
chains_only (all threads to chain-level parallelism, one OpenMP thread
per chain), omp_only (all threads to intra-chain OpenMP operations), and
balanced (even split between chain and OpenMP parallelism).
[`lucifer_parallel_info()`](https://robustecologies.github.io/lucifer/reference/lucifer_parallel_info.md)
produces a diagnostic report covering platform type, core counts, OpenMP
availability, and platform-specific guidance.
`.lucifer_parallel_apply()` is a cross-platform parallel dispatch helper
using PSOCK clusters on all operating systems. `.make_cluster()`
provides managed cluster creation.
[`set_omp_threads_cpp()`](https://robustecologies.github.io/lucifer/reference/set_omp_threads_cpp.md)
enables runtime OpenMP thread count adjustment from R.

`inst/include/lucifer/thread_rng.h` provides a shared xoshiro256+
pseudo-random number generator header for use across OpenMP regions in
C++ code, ensuring independent per-thread streams. SMC particle
rejuvenation is conditionally parallelised via OpenMP when compiled
models are used. VB Pathfinder multi-path execution and ABC rejection
sampling now run in parallel when `CPUs > 1`.

Critical fix: all child chain processes previously had `OMP_NUM_THREADS`
hardcoded to 1, eliminating OpenMP parallelism for distribution
evaluations, diagnostics, and PSIS computations within each chain. The
thread budget is now inherited correctly from the parent process. SBI
simulation via the `n_cores` parameter and Kfold cross-validation via
the `CPUs` parameter now work on Windows through PSOCK dispatch,
replacing the previous Unix-only `mclapply` implementation. Explicit
variable exports replace `ls(.GlobalEnv)` cluster export patterns across
14 files, eliminating a source of unpredictable failures.

[`server_Listening()`](https://robustecologies.github.io/lucifer/reference/server_Listening.md),
the legacy HPC socket server, is deprecated. Use
`lucifer(Chains = n, CPUs = n)` instead.

------------------------------------------------------------------------

## lucifer 3.2.0

### Neural ODE inference via Bayesian gradient matching

[`NODE()`](https://robustecologies.github.io/lucifer/reference/NODE.md)
implements neural ordinary differential equation inference following the
Bayesian gradient matching approach of Bonnaffé (2023). A sinusoidal
observation network captures smooth temporal trajectories; a
dual-pathway single-layer perceptron process model approximates the
latent dynamics. BFGS optimisation with OpenMP-parallelised ensemble
fitting produces a posterior ensemble of neural network weights that
characterises uncertainty in the inferred dynamics. The Geber criterion
drives sparse network architecture selection; RK4 integration performs
forward trajectory generation. The C++ backend in `src/node_bngm.cpp`
(~620 lines, 8 exported functions) achieves approximately 0.2 seconds
for two-species Lotka-Volterra inference, compared to approximately 5
minutes for the pure R reference implementation.
[`NODE_predict()`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md)
and
[`predict.node_fit()`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md)
generate forward forecasts from the posterior ensemble.
[`plot.node_fit()`](https://robustecologies.github.io/lucifer/reference/plot.node_fit.md)
supports eight visualisation types: interpolation, dynamics, effects,
contributions, network diagram, phase portrait, cross-validation
performance, and forecast.

------------------------------------------------------------------------

## lucifer 3.1.0

### Simulation-based inference overhaul

[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
gains five new parameters. `log_prior_fn` accepts a user-supplied exact
log-prior function for NLE and NRE methods, eliminating the multivariate
Gaussian prior approximation that was used previously and produced
biased posteriors for non-Gaussian prior shapes. `max_grad_norm` makes
gradient clipping configurable (was hardcoded at 10.0). `seed` enables
reproducible neural network initialisation. `n_cores` controls parallel
simulation dispatch. `pretrained_network` enables warm-starting from a
previously trained SBI network for transfer learning and sequential
refinement.

Three new diagnostic functions:
[`C2ST()`](https://robustecologies.github.io/lucifer/reference/C2ST.md)
implements the classifier two-sample test, which trains a binary
classifier to distinguish approximate from exact posterior samples and
reports the accuracy as a quality score (0.5 = perfect, 1.0 = maximally
wrong).
[`TARP()`](https://robustecologies.github.io/lucifer/reference/TARP.md)
implements tests of accuracy with random points (Lemos et al. 2023),
which is more statistically powerful than standard simulation-based
calibration for detecting posterior miscoverage.
[`sbi_ppc()`](https://robustecologies.github.io/lucifer/reference/sbi_ppc.md)
is now exported for direct use. Plot types `"tarp"` and `"c2st"` added
to
[`plot.sbi()`](https://robustecologies.github.io/lucifer/reference/plot.sbi.md).
[`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md)
and
[`expected_coverage()`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md)
extended to NLE and NRE via lightweight MCMC per calibration sample.

NPE posterior sampling now uses rejection sampling within the prior
support, eliminating tail artefacts from MDN density extrapolation
outside the prior. SNPE uses the TSNPE truncation scheme with IQR-based
robust bounds, making it immune to heavy-tail bias from outlier
simulations. The prior approximation uses a full multivariate Gaussian
including covariance structure and support bounds, replacing the
previous independent normal marginals.

Critical bug fixed: `.sbi_eval_log_prior` used
[`environment()`](https://rdrr.io/r/base/environment.html) to cache
prior statistics, which returned the ephemeral execution environment
rather than the model closure environment. The cache never persisted
between calls, causing 10,000 prior draws to be generated at every NLE
and NRE model evaluation.

------------------------------------------------------------------------

## lucifer 3.0.1

### Algorithm correctness fixes

Comprehensive testing across 141 algorithm, model, and chain
combinations identified nine algorithm failures requiring correction.
SAMWG, SMWG, USAMWG, and USMWG: the `Specs$Dyn` validator now
auto-generates a character matrix from `parm.names` when Dyn is NULL or
non-character, matching the documented default behaviour. MEADS: a NULL
epsilon guard prevents crashes during startup; NaN guards added at six
points in the dual-averaging adaptation (initial gradients, avg_alpha,
Hbar, epsilon, and epsilonbar). RJ: `bin.p` is coerced to a scalar
before the if() condition that was receiving a vector. SVGD: particle
initialisation widened to `pmax(abs(parm), 1.0)` with a retry loop; the
convergence condition at iter == 1 relaxed to allow genuine fast
convergence.

VB Salimans1, Salimans2, and NGD had an infinite loop when the variance
matrix V collapsed to near-zero: the `while(identical(xstar, m))` loop
would hang indefinitely because Cholesky of V produces zero
displacements. The loop now breaks after 100 attempts with a
noise-escape perturbation of magnitude 1e-4. INCA crashed in
single-machine mode with `get("con")` when the HPC socket connection
objects were absent; a tryCatch fallback with an `.hpc_mode` flag now
skips interchain serialisation in non-HPC mode. Stochastic yhat (from
`rnorm` in the Model function) was inflating sigma in
`robust_influence.R` by approximately sqrt(2), making all observations
appear equally likely and producing uniform PSIS importance weights;
fixed via leave-one-out log-likelihood decomposition:
`ll_i = total_LL - LL_{-i}`. BMK diagnostic gains an NA guard in the Ind
matrix for small sample sizes.

------------------------------------------------------------------------

## lucifer 3.0.0

### Eight new MCMC algorithms

The MCMC algorithm count increases from 54 to 62. **Barker Proposal**
(Livingstone & Zanella 2022): gradient-based proposals from the Barker
balancing function, which is the unique locally balanced proposal
satisfying detailed balance with Laplace increments; dual-averaging
step-size adaptation targets 0.40 acceptance rate. **Ensemble Slice
Sampler / Zeus** (Karamanis & Beutler 2021): gradient-free
affine-invariant ensemble sampler using the doubling-out slice
procedure; immune to linear degeneracies and highly effective for
correlated posteriors. **DREAM** (Vrugt, ter Braak, Diks, Robinson,
Hyman & Higdon 2009): multi-pair differential evolution mutation,
crossover, and snooker updates; well-suited to multi-modal and
high-dimensional problems. **Non-Reversible Parallel Tempering** (Syed,
Bouchard-Côté, Deligiannidis & Doucet 2022): deterministic even-odd
sweep over the replica ladder with direction variables, outperforming
reversible parallel tempering at the same computational cost. **Delayed
Acceptance MCMC** (Christen & Fox 2005): two-stage acceptance using a
cheap surrogate model to pre-screen proposals, with the exact likelihood
evaluated only for pre-accepted proposals; an auto subset-data surrogate
is provided when no custom surrogate is supplied. **Zanella Locally
Balanced** (Zanella 2020): optimal discrete proposals via sqrt- or
min-balancing functions on the proposal weight;
`src/sampler_discrete.cpp` new. **Polya-Gamma Data Augmentation**
(Polson, Scott & Windle 2013): exact Gibbs sampler for logistic and
binomial regression models via Polya-Gamma auxiliary variable
augmentation, producing conjugate Gaussian updates for the regression
coefficients. **Zig-Zag Sampler** (Bierkens, Fearnhead & Roberts 2019):
piecewise deterministic Markov process with componentwise velocity flips
at Poisson event times determined by the gradient of the log-target.

### SSM module and C++ infrastructure

A Kalman filter C++ module (`src/kalman.h` + `src/kalman.cpp`) provides
forward filter, RTS smoother, and backward simulation (Carter-Kohn) with
Joseph-form covariance updates, handling missing observations and
time-varying system matrices. A Polya-Gamma random variate generator
(`src/polyagamma.h` + `src/polyagamma.cpp`) implements the exact Devroye
(2009) method for PG(1,z), integer summation for PG(b,z), and gamma
approximation for non-integer b. A conditional particle filter with
ancestor sampling (`src/pgas.cpp`) extends the bootstrap particle filter
with reference trajectory anchoring.

[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md) is
exported as a top-level dispatcher for state-space model inference,
auto-selecting FFBS, PGAS, SMC2, or KSC based on model structure. The
`ssm_fit` S3 class gains `print`, `summary`, `plot`, and
`as.demonoid.ssm_fit` methods. IBIS mode in
[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
via `Schedule = "ibis"` (Chopin 2002) enables online sequential Bayesian
learning with incremental data batches.
[`coupled_mcmc()`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md)
implements reflection-maximal coupling of any MCMC kernel for unbiased
posterior estimation via the Jacob-O’Leary-Atchadé (2020) debiasing
formula.

### C++ porting (Phase 44b)

Seven existing samplers ported to C++ in a single pass: PG Gibbs
(`mcmc_pg_cpp` in `src/sampler_other.cpp`), Zig-Zag (`mcmc_zigzag_cpp`
in `src/sampler_gradient.cpp`), KSC stochastic volatility
(`ssm_ksc_cpp`), SMC² (`ssm_smc2_cpp`), coupled MCMC
(`coupled_mcmc_cpp`), IBIS batch weights (`smc_ibis_batch_weights_cpp`),
INCA non-HPC path (`mcmc_inca_cpp`), and USMWG/USAMWG (`mcmc_usmwg_cpp`,
`mcmc_usamwg_cpp`). A new dedicated file `src/ssm_engines.cpp` (1,036
lines) houses the KSC and SMC² C++ backends. Total C++ source reaches
approximately 26,262 lines.

------------------------------------------------------------------------

## lucifer 2.12.0

### Major SDE extension

The standalone SDE family functions (`sde_ou()`, `sde_gbm()`, etc.) have
been replaced by a unified family registry: `SDE(family = "ou")`
dispatches over eleven families (ou, vasicek, gbm, cir, lotka_volterra,
sir, fitzhugh_nagumo, merton, kou, heston, and double_well). Each family
entry carries a native C++ identifier enabling OpenMP-parallelised
particle filter dispatch without R callbacks when the native family is
recognised.

The EKF likelihood engine is available via `method = "ekf"`, computing a
deterministic approximate log-likelihood through continuous-discrete
Kalman filtering. A symbolic differentiation engine
(`.sde_symbolic_deriv()`) computes analytical drift and diffusion
derivatives by recursive AST traversal; the companion
`.sde_stratonovich_correction()` computes the Ito-Stratonovich
noise-induced drift g(x)·g’(x)/2 symbolically. Setting
`interpretation = "stratonovich"` in
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
triggers automatic correction of the drift when the model is specified
in Stratonovich form.

[`compile.sde_model()`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md)
compiles drift and diffusion R expressions to native C++ via
[`sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) with
hash-based caching. `R/sde_compile.R` provides the recursive AST parser
and C++ code generator. `inst/include/lucifer/sde_native.h` provides an
xoshiro256+ thread-safe PRNG for compiled particle propagation.

Additional particle filter variants: auxiliary particle filter
(`pf.type = "auxiliary"`, Pitt & Shephard 1999) applies look-ahead
weighting for better ESS with informative observations; bridge particle
filter (`pf.type = "bridge"`, Delyon & Hu 2006) uses guided proposals
toward the next observation for sparse data.
[`smooth.sde_fit()`](https://robustecologies.github.io/lucifer/reference/smooth.sde_fit.md)
runs the forward-filtering backward-sampling smoother for joint
smoothing.
[`predict.sde_fit()`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md)
generates posterior predictive forecasts.
[`SDE.fit.pmmh()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.pmmh.md)
is a convenience wrapper for pseudo-marginal MCMC.
[`log_lik.sde_fit()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
extracts pointwise log-likelihoods for LOO-PSIS model comparison. The
N-species generalised Lotka-Volterra family auto-detects the species
count from `ncol(data)`, removing the two-species constraint. A
community ecology family with hierarchical priors and configurable
interaction structure is added. SDE `yhat` now returns predicted or
filtered state means; SDE `Monitor` carries log-likelihood, last
filtered state per variable, and minimum ESS.

### SSM module

[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
is an exported constructor for discrete-time state-space model
specification without requiring
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md).
Five builder types are registered: `local_level`, `local_linear_trend`,
`seasonal` (configurable period), `bsm` (Harvey 1989 basic structural
model), and `var_p` (VAR(p) in companion form with Minnesota prior). Six
non-Gaussian observation families registered in `R/ssm_obs_families.R`:
Gaussian, Poisson, negative-binomial, Student-t, zero-inflated Poisson,
binomial. The UKF (`src/ukf.h` + `src/ukf.cpp`, ~300 C++ lines)
propagates sigma points through nonlinear dynamics without Jacobians.
The EnKF (`src/enkf.cpp`, ~350 lines) is a stochastic ensemble Kalman
filter with perturbed observations and covariance inflation. The RBPF
(`src/rbpf.cpp`, ~260 lines) maintains per-particle Kalman filters for
conditionally linear substates. MS-FFBS implements K-regime
Markov-switching with the Hamilton filter in C++ and R Gibbs updates.
Total SSM algorithm count: 8.

------------------------------------------------------------------------

## lucifer 2.11.0

### Stan/JAGS interoperability bridge

[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
is an S3 generic that imports stanfit, brmsfit, CmdStanMCMC, mcmc.list,
runjags, and posterior::draws objects into the demonoid class, making
the full lucifer post-processing pipeline available for fits from
external Bayesian engines.
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
and
[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)
are backend fitting functions that fit Stan or JAGS models and return
demonoid objects with full compatibility for Kfold and LFO
cross-validation.
[`to_mcmc_list()`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md),
[`to_draws_array()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md),
[`to_draws_matrix()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md),
and
[`to_draws_df()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md)
export lucifer fits to coda and posterior package formats.
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
is a new S3 generic for uniform pointwise log-likelihood extraction from
any fit object.
[`stan_to_spec()`](https://robustecologies.github.io/lucifer/reference/stan_to_spec.md)
and
[`jags_to_spec()`](https://robustecologies.github.io/lucifer/reference/jags_to_spec.md)
are experimental translators mapping Stan/JAGS model code to lucifer
`model_spec` IR.

Namespace masking by brms and loo (which export `log_lik`, `LOO`,
`loo_compare`, `Rhat`) is handled via `.onLoad`/`.onAttach` hooks that
register `log_lik.demonoid` on brms and loo generics through
[`registerS3method()`](https://rdrr.io/r/base/ns-internal.html) and
`setHook(packageEvent())`. Generated quantities arrays such as
`log_lik[1]...log_lik[N]` are filtered from the draws before building
the demonoid posterior matrix.

The `inst/include/lucifer/dist_gradients.h` and
`inst/include/lucifer/transforms.h` header-only libraries are finalised
with 24 distributions carrying analytical log-pdf gradients and five
bijective constraint transforms with log-Jacobian and gradient
functions; these headers are included in generated C++ model code at
zero overhead.

------------------------------------------------------------------------

## lucifer 2.10.0

### Frequentist bridge

[`freq_summary()`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)
extracts a frequentist coefficient table (estimate, standard error, z or
t statistic, p-value, and confidence interval) from any lucifer fit
object (demonoid, laplace, iterquad, or data_cloning) and returns a S3
class with `print` (significance stars) and `plot` (forest plot with
theme_relab) methods.
[`wald_test()`](https://robustecologies.github.io/lucifer/reference/wald_test.md)
performs Wald tests for individual parameters or joint hypotheses
specified as contrast matrices.
[`lr_test()`](https://robustecologies.github.io/lucifer/reference/lr_test.md)
computes the likelihood ratio test between nested models using the
maximised log-likelihood from the Dev component.
[`score_test()`](https://robustecologies.github.io/lucifer/reference/score_test.md)
provides a score (Lagrange multiplier) test via finite-difference
gradient at the restricted estimate.
[`freq_ic()`](https://robustecologies.github.io/lucifer/reference/freq_ic.md)
extracts AIC, AICc, and BIC.
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
produces a forest plot comparing Bayesian credible intervals and
frequentist confidence intervals from multiple fit objects.
[`coverage_sim()`](https://robustecologies.github.io/lucifer/reference/coverage_sim.md)
runs a Monte Carlo simulation to assess empirical coverage. All
functions are implemented in `R/freq_tools.R` (~1,000 lines) with S3
methods.

[`profile_likelihood()`](https://robustecologies.github.io/lucifer/reference/profile_likelihood.md)
computes profile log-likelihood curves for any lucifer fit object: each
parameter is fixed at a grid of values, nuisance parameters are
optimised via BFGS, and the profile confidence interval is determined
from the chi-squared threshold. An asymmetry ratio diagnoses the
adequacy of the Wald approximation. The S3 class has `print` (CI
comparison table), `summary` (per-parameter asymmetry diagnostics), and
`plot` (faceted ggplot2 with profile and Wald bounds) methods.

[`freq_residuals()`](https://robustecologies.github.io/lucifer/reference/freq_residuals.md)
produces four-panel residual diagnostics (residuals vs fitted, Q-Q plot,
scale-location, leverage and Cook’s distance) from any lucifer fit
object; the `print` method reports a Shapiro-Wilk normality test. The RJ
between-model birth proposal is fixed to use `nonzero_post(v_change)`
when constructing proposals for dead parameters, rather than always
proposing from the zero value, so that birth proposals are centred on
the last known nonzero value as documented.

------------------------------------------------------------------------

## lucifer 2.9.0

### Cross-validation framework

[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
provides unified model comparison across LOO, WAIC, K-fold, and LFO
criteria, ranking models by ELPD and computing correlated standard
errors of pairwise differences, with S3 `print` and `plot` (forest plot)
methods.
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
performs exact K-fold cross-validation by refitting the model K times,
supporting custom fold assignments including blocked, group-level, and
spatial partitions.
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
implements leave-future-out cross-validation for time series (Bürkner,
Gabry & Vehtari 2020) with PSIS approximation and automatic refitting
when Pareto k exceeds a threshold; the S3 class includes two-panel plots
combining cumulative ELPD and Pareto k diagnostics.
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md)
computes Bayesian model stacking and pseudo-BMA weights (Yao, Vehtari,
Simpson & Gelman 2018) via three methods: constrained softmax
optimisation, pseudo-BMA, and pseudo-BMA+.
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
now returns a `waic` class with a `pointwise` data frame preserving
backward compatibility.

### Algorithm orchestration

[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
profiles a Bayesian model across six dimensions (dimensionality, speed,
gradient availability, parameter constraints, conditioning, and
multimodality) and recommends an inference strategy across MCMC, VB,
Laplace, IQ, PMC, SMC, ABC, and SBI.
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
benchmarks any collection of fit objects across ESS/second, marginal
KLD, Wasserstein-1 distance, and a Pareto frontier, with eight plot
types.
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
orchestrates the full Prescribe → Fit → Consort → Refine → Arena
pipeline with configurable rounds and families.
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
v2 returns an S3 `consort` object with structured `next_specs` for
programmatic use; appeasement thresholds upgraded to ESS bulk ≥ 400,
tail ESS ≥ 200, and split Rhat \< 1.01. Algorithm registry
`R/algo_registry.R` covers 90+ methods with uniform property fields.
`Juxtapose()` is soft-deprecated in favour of
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md).

[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
implements maximum likelihood estimation via Bayesian MCMC (Lele, Dennis
& Lutscher 2007): the log-likelihood is multiplied by K to concentrate
the posterior on the MLE, with the asymptotic covariance estimated from
K times the posterior variance. The Campbell-Lele (2014) ANOVA
estimability test detects structural non-identifiability. Ten plot types
are provided via S3 `plot.data_cloning`.

### Critical bug fixes

The C++ RNG starvation bug is fixed: `PutRNGstate()` is now called
before every `Rcpp::Function` callback to prevent the C-level RNG state
from being overwritten, which caused all C++ samplers to reuse identical
random draws within each iteration. NUTS acceptance tracking is
corrected from a conflation of tree-doubling counts with iteration-level
acceptance. Non-finite `log_alpha` is handled correctly across all C++
sampler files (26 instances): NaN is now mapped to negative infinity
(reject) rather than zero (always accept), preserving correct ±Inf
semantics. The numerical gradient is switched from first-order forward
differences to second-order central differences. ADVI full-rank Cholesky
gains a diagonal floor of 1e-7. QMC Lagrange basis gains a
division-by-zero guard for coincident knots. Bayesian quadrature
Cholesky rank-1 update diagonal floor raised from 1e-12 to 1e-8.
[`CovEstim()`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
Hessian method fixed from `-as.inverse(H)` to `as.inverse(-H)` so the
positive-definite `-H` is inverted directly. Hessian Richardson d
parameter corrected from 0.0001 to 0.1 for second derivatives.

------------------------------------------------------------------------

## lucifer 2.8.0

### Model specification DSL

[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
is a declarative probabilistic model specification language that
compiles plain-text or LaTeX notation into
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)-compatible
`Model` functions. The compiler is split across five files:
`R/model_spec.R` (main function and S3 methods), `R/model_spec_parse.R`
(parser), `R/model_spec_ir.R` (intermediate representation builder),
`R/model_spec_compile.R` (code generator), and `R/model_spec_registry.R`
(distribution registry covering 86 distributions including all 76
lucifer density functions as aliases). 81 yhat templates are registered
for posterior predictive generation. Cholesky and precision
parameterisation variants are registered for all multivariate
distributions. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
generic extracts generated source code. `model_spec` objects are
auto-unwrapped by
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

### C++ model compilation backend

`inst/include/lucifer/dist_gradients.h` provides a header-only
analytical log-pdf gradient library for 24 distributions (Normal,
NormalV, NormalP, HalfCauchy, HalfNormal, HalfT, Gamma, InvGamma, Beta,
Exponential, Uniform, Poisson, Bernoulli, Binomial, NegBinomial,
StudentT, LogNormal, Laplace, Weibull, Pareto, InvChiSq, InvGaussian,
LaplaceP, LogNormalP). Each function returns the log-pdf and partial
derivatives with respect to all arguments.
`inst/include/lucifer/transforms.h` provides five bijective constraint
transforms (log, logit, scaled logit, stick-breaking simplex, Cholesky
factor) with log-Jacobian and gradient functions.

`R/model_spec_compile_cpp.R` generates self-contained C++ source from
model_spec IR, supporting 24 distributions, parameter transforms through
all five bijections, and chain rule composition through exp, invlogit,
log, and sqrt. `R/model_spec_runtime.R` compiles generated C++ via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) with
hash-based caching, graceful R fallback on compilation failure, a
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
S3 generic, and a `compiled_model` class. Parallel chains support for
compiled models: serialisable metadata extracted via
`.compiled_model_metadata()` and recompiled in child processes via
`.recompile_in_child()` using the cached source. `src/Makevars` and
`src/Makevars.win` gain `-I../inst/include` for header access.

### NUTS diagnostics and dual number AD

[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
(`R/nuts_diagnostics.R`) provides Stan-like NUTS diagnostics: divergent
transition count, percentage, and iteration indices; E-BFMI with a 0.3
threshold; and tree depth saturation detection. NUTS in
`src/sampler_gradient.cpp` now tracks per-iteration energy, tree depth,
divergent flag, and leapfrog count, returning them as named list
elements. `inst/include/lucifer/dual.h` implements forward-mode
automatic differentiation via dual numbers (~170 lines) with arithmetic
operators and 20+ elementary functions including exp, log, sqrt, pow,
sin, cos, tanh, invlogit, logit, and softplus.

------------------------------------------------------------------------

## lucifer 2.7.0

### Extended visualisation system

`R/plot_mcmc.R` provides nine algorithm-aware MCMC diagnostic functions:
rank histograms, R-hat bar plot, ESS ratio plot, energy diagnostic,
divergence scatter plot, NUTS diagnostics, pairs plot, parallel
coordinates, and ridgeline density areas. `R/plot_prior.R` provides
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
(four types),
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md),
and
[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md).
`R/plot_dag.R` provides
[`plot_dag()`](https://robustecologies.github.io/lucifer/reference/plot_dag.md)
for directed acyclic graph visualisation of model structure via igraph
and ggraph. `R/plot_convenience.R` provides 16 bayesplot-compatible
convenience wrappers: `mcmc_trace`, `mcmc_dens`, `mcmc_intervals`,
`mcmc_areas`, `mcmc_pairs`, `mcmc_parcoord`, `mcmc_rank`, `mcmc_rhat`,
`mcmc_neff`, `mcmc_energy`, `ppc_dens_overlay`, `ppc_intervals`,
`ppc_ribbon`, `ppc_rootogram`, `ppc_stat`, and `ppc_loo_pit`.

Fourteen new PPC styles added to the shared `.plot_ppc()`
implementation: density overlay, histogram overlay, ECDF overlay, stat,
stat 2D, scatter average, error scatter, error histogram, intervals,
ribbon, bars, rootogram, LOO-PIT, and violin grouped. Group faceting is
handled by `.apply_group_facet()`. The `type` argument with
[`match.arg()`](https://rdrr.io/r/base/match.arg.html) dispatch has been
added to all S3 plot methods: `plot.demonoid` gains 14 types, `plot.vb`
6 types, and `plot.pmc`, `plot.laplace`, `plot.iterquad`, and
`plot.bayesquad` are extended. A `ground_truth` overlay argument is
added across all plot methods, drawing vertical lines for densities,
diamonds for interval plots, and horizontal lines for trace plots.
Informative subtitle generation is centralised via `.build_subtitle()`.

------------------------------------------------------------------------

## lucifer 2.6.0

### ProductSpace model selection

[`ProductSpace()`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md)
implements Bayesian model selection via the product space method (Carlin
& Chib 1995). The hybrid sampler uses block RAM (Vihola 2012) for the
active model parameters, exact enumeration over the model index M, and
direct sampling from pseudopriors for inactive model parameters,
reducing per-iteration cost from approximately 10·Σdₖ + K to K + 1 model
evaluations. Bisection calibration for pseudoprior estimation uses the
same hybrid RAM sampler, maintaining cost efficiency throughout
calibration. Multi-chain support is provided via callr with independent
adaptation states per chain. The S3 class `product_space` includes
`print`, `summary`, `plot` (types: probabilities, trace, transition,
bayes_factors, posteriors), and `is.product_space`.

### RobustBayes sensitivity analysis

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
combines four complementary robustness analyses. Power-scaling
sensitivity (Kallioinen, Paananen, Bürkner & Vehtari 2024) measures how
the posterior changes under perturbation of the prior or likelihood.
PSIS-LOO observation influence identifies observations with
disproportionate impact on the posterior. The Evans-Moshonov test
detects prior-data conflict by comparing the prior predictive
distribution with the observed marginal likelihood. An aggregate
robustness score synthesises the four assessments.
`src/robust_bayes.cpp` implements the C++ backend for power-scaling
divergence and calibration. The S3 class `robust_bayes` has `print`,
`summary`, and `plot` methods covering power, influence, conflict,
divergence, aggregate, and power_density visualisation types.

### PSIS bug fixes

Three critical bugs in `src/psis_internal.h` are corrected. The GPD
theta grid used only positive values, preventing detection of heavy
tails in the Pareto k estimator; the correct Zhang-Stephens formula
`theta = 1/x[N] + (1 - sqrt(M/(j-0.5)))/3/xstar` allowing negative theta
is ported from loo v2.9. GPD fitting was applied to log-scale
exceedances rather than natural-scale exceedances; the correct procedure
exponentiates the tail, fits GPD in natural scale, smooths, and logs
back. The profile likelihood formula in the GPD fit produced −Inf for
all valid heavy-tail theta values; replaced with the correct `lx`
function from the loo package.

------------------------------------------------------------------------

## lucifer 2.5.0

### Pathfinder variational inference

`VariationalBayes(Method = "Pathfinder")` implements Pathfinder (Zhang,
Carpenter, Gelman & Vehtari 2022): quasi-Newton variational inference
that constructs approximations along an L-BFGS trajectory and selects
the best approximation using Pareto smoothed importance sampling.
`src/pathfinder.cpp` provides `lbfgs_inverse_hessian_cpp()` implementing
the two-loop L-BFGS recursion with explicit inverse Hessian recovery and
`psis_resample_cpp()` for PSIS smoothing followed by systematic
resampling. Both `src/psis.cpp` and `src/pathfinder.cpp` include the
shared `src/psis_internal.h` header, eliminating the duplicate static
function definitions that previously existed.
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
now supports eight methods. Pathfinder diagnostics, including ELBO along
the L-BFGS trajectory per path and Pareto k assessment, are included in
`print.vb`, `summary.vb`, and `plot.vb`.

------------------------------------------------------------------------

## lucifer 2.4.0

### QMC ensemble methods

Eight quadratic Monte Carlo ensemble MCMC methods from Militzer (2023,
2025) are added: QMC, QMC3, QMCN, DQMC, SQMC, MAMC, SAMC, and WMC.
`src/sampler_qmc.cpp` (~650 lines) implements C++ backends for all
eight, providing Lagrange interpolation helpers, ensemble
initialisation, and guide selection procedures. `R/sampler_qmc.R`
provides thin R wrappers. Specs validation entries and dispatcher
branches are added for all eight methods.

### Bayesian quadrature

[`BayesianQuadrature()`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)
implements GP-based probabilistic numerical integration with four
algorithms. Vanilla BQ (O’Hagan 1991) treats integration as a Bayesian
inference problem over the integrand under a squared-exponential GP
prior. WSABI-L (Gunter, Osborne, Garnett, Hennig & Roberts 2014) warps
the integrand by a square-root transformation and linearises for
tractable GP updates. FWBQ (Bach 2017) uses Frank-Wolfe optimisation for
batch acquisition of integration nodes. BatchBQ extends FWBQ with greedy
batch node selection via OpenMP parallelisation.
`src/bayesian_quadrature.cpp` implements SE kernel matrices, kernel mean
embeddings, GP posterior updates, Cholesky rank-1 updates, and
acquisition functions. The S3 class `bayesquad` has `print`, `summary`,
`plot`, and `predict` methods.

### C++ ports and bug fix

Pseudo-marginal MCMC (`mcmc_pmcmc_cpp` in `src/sampler_metropolis.cpp`),
SGLD (`mcmc_sgld_cpp` in `src/sampler_gradient.cpp`), and
reversible-jump MCMC (`mcmc_rj_cpp` in `src/sampler_other.cpp`) are
ported from R to C++. A `LaplaceApproximation` parallel-chain crash is
fixed: when all initial values are zero, each callr subprocess
independently attempted LA and chain 1 crashed, killing the entire
parallel run. LA now runs once in the parent process before chain
dispatch, with a `tryCatch` wrapper on the single-chain path as well.

------------------------------------------------------------------------

## lucifer 2.3.0

### Automatic full-conditional Gibbs sampling

`lucifer(Algorithm = "Gibbs", Specs = NULL)` now constructs univariate
sampling kernels automatically, eliminating the requirement for
user-supplied closed-form full-conditional distributions. `R/auto_fc.R`
implements Neal (2003) stepping-out slice sampling for continuous
parameters and exhaustive enumeration for discrete parameters specified
via `Data$dparm` (parameter indices) and `Data$dsupport` (named list of
support values). The factory `.make_auto_fc()` constructs per-parameter
kernels at runtime by inspecting the model structure.

### Simulation-based inference

[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
provides simulation-based inference with three method families. NPE
(neural posterior estimation) trains a mixture density network
(`src/mdn.cpp`, ~560 lines) on simulated prior-parameter pairs to
directly estimate the posterior. NLE (neural likelihood estimation)
trains a density network on simulated parameter-data pairs and plugs the
trained network as a likelihood into any of the 62 MCMC algorithms. NRE
(neural ratio estimation) trains a binary classifier
(`src/nre_classifier.cpp`, ~350 lines) to estimate the
likelihood-to-prior ratio. Sequential variants SNPE, SNLE, and SNRE
refine the approximation through rounds of targeted simulation.
[`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md)
performs simulation-based calibration with rank histograms and KS test;
[`expected_coverage()`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md)
computes HPD coverage diagnostics. S3 class `sbi` with `print`,
`summary`, and `plot` methods covering posterior, pairs, training, sbc,
coverage, and ppc plot types.

### SA-ABC and callr parallel chains

[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
gains a fourth algorithm: simulated annealing ABC (SA-ABC). A Boltzmann
soft-kernel annealing phase (30% of iterations) is followed by a
hard-threshold sampling phase (70%). T0 is auto-calibrated from the
median pilot distance; the cooling rate is computed so the temperature
reaches epsilon at the end of Phase 1. Phase 2 reinitialises at a point
within epsilon to avoid frozen chains.
[`plot.abc()`](https://robustecologies.github.io/lucifer/reference/plot.abc.md)
gains `type = "temperature"` for SA-ABC diagnostics.

Parallel chain execution is rewritten with
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html) as the
primary dispatch (`.run_chains_callr()`), with PSOCK fallback
(`.run_chains_psock()`). Per-chain progress bars use a C++ progress file
mechanism (`set_progress_file_cpp()`/`clear_progress_file_cpp()`), with
terminal-aware rendering: ANSI multi-line block bars with cursor-up
overwrite for real terminals and a single-line aggregate summary for
RStudio. The previous `.clusterApplyLB_progress()` implementation, which
used unexported `parallel` internals `sendCall`/`recvOneData`, is
removed. Console output across the package is unified via `.log_msg()`
and `get_colored_symbol()`.
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
defaults changed to `Chains = 3L`, `CPUs = 3L`, `Type = NULL`
(auto-detects FORK on Unix, PSOCK on Windows).

------------------------------------------------------------------------

## lucifer 2.2.0

### Variational Bayes expansion

[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
now supports seven general-purpose variational inference algorithms via
the `Method` argument, up from the single Salimans2 method. The new
methods span the modern VI landscape from fixed-form stochastic
approximation to particle-based inference.

**Salimans1** implements the first algorithm of Salimans and Knowles
(2013), using only gradient information. The precision matrix is updated
via the outer product of gradients, reducing per-iteration cost from
O(J²/2) to O(J+1) model evaluations. **ADVI.mf** implements mean-field
automatic differentiation variational inference (Kucukelbir et al. 2017)
with the reparameterisation trick and Adam optimisation. **ADVI.fr** is
the full-rank variant, parameterising the covariance through a
lower-triangular Cholesky factor; a warning is issued when J \> 50.
**BBVI** implements black-box variational inference (Ranganath et al.
2014) using the score function estimator with per-component control
variate baselines. **SVGD** implements Stein variational gradient
descent (Liu & Wang 2016), a particle-based method that updates K
particles via the functional gradient of KL divergence in a reproducing
kernel Hilbert space. **NGD** implements natural gradient variational
inference (Amari 1998; Khan & Lin 2017) with a Robbins-Monro schedule.

Shared infrastructure: numerically stable softplus and sigmoid
functions, an Adam optimiser, and `.clip_grad()` (max_norm = 10). Four
bug fixes: the reparameterisation trick no longer overwrites theta with
the model-evaluation parm, which broke gradient consistency; the Adam
v-hat is no longer permanently corrupted by gradient explosions from
constraint boundaries; false early convergence from noisy ELBO is
prevented by a patience mechanism requiring 20 consecutive passes with
min_iter; LP \< −1e100 at the initial parameter value is now caught and
reported.

------------------------------------------------------------------------

## lucifer 2.1.0

### LOO-PSIS cross-validation

[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
computes approximate leave-one-out cross-validation via Pareto smoothed
importance sampling, following Vehtari, Gelman and Gabry (2017) and
Vehtari, Simpson, Gelman, Yao and Gabry (2024). The function accepts the
same N × S pointwise log-likelihood matrix as
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
and returns `elpd_loo`, `p_loo`, and `looic` with standard errors, plus
per-observation Pareto k shape parameters that diagnose reliability. GPD
fitting uses the Zhang and Stephens (2009) method with the weakly
informative prior adjustment. `src/psis.cpp` (~230 lines) parallelises
over observations via OpenMP.
[`PSIS()`](https://robustecologies.github.io/lucifer/reference/PSIS.md)
provides standalone Pareto smoothed importance sampling for general
importance ratio stabilisation. `print.loo` displays the LOO table with
Pareto k threshold counts. `plot.loo` produces a ggplot2 scatter of
Pareto k values coloured by reliability category. `summary.demonoid` and
`print.demonoid` report LOO-PSIS when Monitor contains pointwise
log-likelihoods.

------------------------------------------------------------------------

## lucifer 2.0.0

### Package renamed from LaplacesDemon to lucifer

The primary inference function is renamed from `LaplacesDemon()` to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
`LaplacesDemon.hpc()` becomes `lucifer.hpc()` and `LaplacesDemon.RAM()`
becomes
[`lucifer.RAM()`](https://robustecologies.github.io/lucifer/reference/lucifer.RAM.md).
The `demonoid` S3 class and
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
retain their names for backward compatibility.

### Stochastic differential equation interface

[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
provides structured model specification for stochastic dynamical
systems. Users specify drift, diffusion, observation model, and priors;
the package generates a
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)-compatible
`Model` function backed by C++ likelihood engines. Four likelihood
engines are available: exact transition densities for
Ornstein-Uhlenbeck, geometric Brownian motion, Cox-Ingersoll-Ross, and
Vasicek processes; Euler-Maruyama discretisation for general nonlinear
SDEs; the Milstein scheme with higher-order correction; and a bootstrap
particle filter producing unbiased log-marginal-likelihood estimates for
pseudo-marginal MCMC.

Observation models include Gaussian, Poisson, binomial, and user-defined
likelihoods. Jump-diffusion processes are supported via optional
`jump.rate` and `jump.size` arguments. Eight pre-built families provide
convenience wrappers: `sde_ou()`, `sde_gbm()`, `sde_cir()`,
`sde_lotka_volterra()`, `sde_sir()`, `sde_fitzhugh_nagumo()`,
`sde_merton()`, and `sde_kou()`.
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
bundles model specification with MCMC inference, returning an `sde_fit`
S3 class with trajectory, predictive, phase portrait, residuals,
posterior, and filtered state plot types.
[`simulate.sde_model()`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md)
provides forward simulation using exact transition densities or
Euler-Maruyama.

New C++ files: `src/sde_engines.cpp` (~370 lines, 6 exported functions)
and `src/particle_filter.cpp` (~260 lines, 2 exported functions). An SMC
bug is fixed: `arma::vec` returned from C++ was silently promoted to a
1-column matrix, causing
[`cov.wt()`](https://rdrr.io/r/stats/cov.wt.html) to fail with
“non-conformable arrays” when the tempering schedule jumped to beta = 1
in a single stage.

------------------------------------------------------------------------

## lucifer 1.4.0

### Sequential Monte Carlo and Approximate Bayesian Computation

[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
implements adaptive tempering sequential Monte Carlo with systematic
resampling, random-walk Metropolis rejuvenation, and online covariance
adaptation. The C++ backend in `src/smc.cpp` handles resampling and
weight computation.
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
implements approximate Bayesian computation with three algorithms:
rejection sampling, MCMC-ABC (Marjoram, Molitor, Plagnol & Tavaré 2003),
and SMC-ABC (Sisson, Fan & Tanaka 2007; Beaumont, Cornuet, Marin &
Robert 2009). Both return S3 classes with `print`, `summary`, and `plot`
methods.

Four new MCMC algorithms are added, bringing the total to 47. Stochastic
gradient HMC (SGHMC; Chen, Fox & Guestrin 2014) extends SGLD with
momentum. Microcanonical HMC (MCHMC; Robnik, De Luca & Ramirez-Ruiz
2022) uses an energy-conserving microcanonical ensemble. Dual-chain
adaptive HMC (MEADS; Hoffman & Sountsov 2022) estimates the mass matrix
from dual chains running in parallel. Pseudo-marginal MCMC uses an
unbiased likelihood estimator (typically a particle filter) as a noisy
oracle, producing exact posterior samples for intractable likelihoods.

Two NUTS and HMCDA bug fixes: momentum initialisation changed from
`R::runif` to `R::rnorm`; dual-averaging now correctly accumulates
`alphaprime` and `nalphaprime` across tree-doubling iterations.

------------------------------------------------------------------------

## lucifer 1.3.0

### Extended distributions

Nine new distribution families are added with C++ implementations in
`src/distributions_extended.cpp` (~650 lines, 23 C++ functions) and R
wrappers in `R/distributions.R` (~700 lines, 28 R functions).

Tier 1 covers high-demand univariate distributions. The ex-Gaussian
(`dexgaussian`, `pexgaussian`, `rexgaussian`) models reaction-time data
as the convolution of a normal and an exponential. The von Mises
(`dvonmises`, `pvonmises`, `rvonmises`) models circular data using the
Best-Fisher exact sampling algorithm. The stable distribution
(`dstable`, `pstable`, `rstable`) uses Chambers-Mallows-Stuck random
generation and numerical characteristic function inversion for density
evaluation. The Tweedie (`dtweedie`, `rtweedie`) is computed via the
Dunn-Smyth (2005) series expansion for density and compound
Poisson-gamma simulation.

Tier 2 covers copulas and specialised distributions. Clayton
(`dclayton`, `rclayton`), Gumbel (`dgumbel.copula`, `rgumbel.copula`),
and Frank (`dfrank`, `rfrank`) copulas for dependence modelling; the
simplex distribution (`dsimplex`, `rsimplex`) for bounded continuous
data on (0,1); the matrix-t (`dmatrixt`, `rmatrixt`) for matrix-variate
heavy-tailed modelling.

Tier 3 covers spatial distributions. The intrinsic CAR prior (`dicar`,
`ricar`) for spatially structured random effects on irregular lattices,
using eigendecomposition-based pseudo-determinant computation. Matern
correlation and covariance functions (`matern_corr`, `matern_cov`) with
OpenMP-parallelised covariance matrix construction.

------------------------------------------------------------------------

## lucifer 1.2.0

### User-supplied gradients

[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
now accepts a `Grad` argument for gradient-based MCMC algorithms (HMC,
AHMC, HMCDA, MALA, NUTS, Refractive, SGLD, THMC). Three modes are
supported: `Grad = NULL` (default) uses the built-in forward
finite-difference approximation via `partial_cpp`;
`Grad = function(parm, Data)` accepts a user-supplied analytical
gradient function that is validated at startup; and `Grad = "numDeriv"`
wraps [`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html)
with Richardson extrapolation for higher-accuracy gradients at adaptive
step size. Non-gradient algorithms ignore the `Grad` parameter. The
dispatch function `compute_gradient()` in `src/sampler_common.cpp`
routes to C++ finite differences or the R gradient function depending on
the mode. All seven gradient-based C++ samplers accept a
`Rcpp::Nullable<Rcpp::Function> GradFn` parameter. The gradient function
is forwarded through multi-chain execution when `Chains > 1`.

------------------------------------------------------------------------

## lucifer 1.1.0

### Modern convergence diagnostics

[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
implements split R-hat (Vehtari, Gelman, Simpson, Betancourt & Gabry
2021) with Blom rank-normalisation, dispatching on vectors, matrices,
lists, and `demonoid` objects.
[`ESS.bulk()`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md)
computes rank-normalised bulk effective sample size.
[`ESS.tail()`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md)
computes folded rank-normalised tail ESS.
[`MCSE.quantile()`](https://robustecologies.github.io/lucifer/reference/MCSE.quantile.md)
provides MCSE for quantiles via the Flegal and Jones (2010) batch method
with sqrt and cuberoot batch sizes. Divergence tracking is added to
seven gradient-based C++ samplers (HMC, AHMC, HMCDA, MALA, NUTS, THMC,
and Refractive) using a threshold of \|ΔH\| \> 10 nats.
[`summary.demonoid()`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md)
is rewritten to report divergence counts, split R-hat with quality
assessment (excellent, adequate, or not converged), bulk and tail ESS,
alongside standard ESS.

------------------------------------------------------------------------

## lucifer 1.0.1

### Bug fixes

Sixteen algorithmic correctness bugs are corrected across
`specs_validation.R`, `postprocess.R`, and C++ samplers. B1: Gibbs
algorithm assignment used `==` for comparison instead of `<-` for
assignment. B2: shadow variable scoping instead of Specs mutation at
four locations in CHARM/HARM/HMC validation. B3: orphaned if-guard in
AHMC Specs validation with no body. B4: `Posterior2` dimension drop
without `drop = FALSE` in `postprocess.R`. B5: floating-point equality
in HMC momentum reflection replaced with `std::abs() < 1e-14`. B7: 17
unchecked `inv_sympd` calls in `distributions_multivariate.cpp` gain
two-argument form with bool check and eigendecomposition fallback. B8:
division by zero in `col_vars_cpp`/`row_vars_cpp` when n = 1. B9:
`cloglog` NaN at p = 1; `logadd` log(0) when inputs are nearly equal,
fixed with log1p. B10: `as_inverse_cpp` eigenvalue regularisation
improved for near-singular matrices. B11: random index bias from
`ceil(runif)` replaced with `floor(runif)`. B12: slice sampler upper
stepping asymmetry corrected. B14: parallel chains over-broad
environment export scoped to required objects. B15: Gelman diagnostic
silent failure replaced with a warning. B16: `Initial.Values` matrix
column count now validated. Gibbs-to-MWG fallback returns the full name
“Metropolis-within-Gibbs”. MALA defaults gain the missing gamma
parameter.

`specs_validation.R` is refactored from a 954-line cascading if-else
chain to a 735-line table-driven architecture with five declarative
registries and 34 per-algorithm validators.

------------------------------------------------------------------------

## lucifer 1.0.0

### C++ backend

Computationally intensive operations have been rewritten in C++ via Rcpp
and RcppArmadillo with C++17 and OpenMP parallelisation for multicore
acceleration.

`src/matrices.cpp` provides 13 matrix utility functions including
`is_positive_definite_cpp`, `as_positive_definite_cpp`,
`as_symmetric_matrix_cpp`, `tr_cpp`, `as_parm_matrix_cpp`,
`GaussHermiteQuadRule`, `Cov2Prec`, and `Prec2Cov`. `src/math.cpp`
provides 8 mathematical functions: `logit`, `invlogit`, `cloglog`,
`invcloglog`, `interval`, `logadd`, `partial`, and `softmax`.
`src/distributions_multivariate.cpp` (~580 lines, 44 exported functions)
and `src/distributions_univariate.cpp` (~540 lines, ~40 exported
functions) cover density, distribution, quantile, and random generation
for all distribution families. OpenMP parallelises multivariate
distribution evaluations over observations. `src/diagnostics.cpp` (~330
lines) implements ESS via Levinson-Durbin AR fitting with AIC order
selection, MCSE (IMPS and batch methods), IAT, and BMK Hellinger
distance via Gaussian KDE, all OpenMP-parallelised over posterior
columns.

MCMC sampler inner loops are ported to C++ across six source files
organised by family. `src/sampler_metropolis.cpp`: RWM, HARM, CHARM,
OHSS, RSS, RDMH, SAMWG, SMWG. `src/sampler_adaptive.cpp`: AM, AMM, AMWG,
AAMM, ASAMWG, ASMWG. `src/sampler_gradient.cpp`: HMC, AHMC, HMCDA, MALA,
NUTS, THMC, Refractive. `src/sampler_slice.cpp`: Slice, UESS.
`src/sampler_ensemble.cpp`: AIES, DEMC, twalk. `src/sampler_other.cpp`:
AFSS, ESS, IM, MWG, ADMG, RAM. Nine algorithms remain in pure R due to
their inherent R-level dependencies: Gibbs, GG variants, AGG, SGLD,
INCA, MCMCMC, RJ, USAMWG, USMWG. Shared infrastructure is in
`src/sampler_common.h` and `src/sampler_common.cpp`. The package
requires C++17 and detects OpenMP availability automatically at compile
time.

### Monolith decomposition

The original `R/LaplacesDemon.R` (8,619 lines) is decomposed into
`R/lucifer.R` (dispatcher, 487 lines), `R/specs_validation.R` (algorithm
specification validation), `R/postprocess.R` (post-MCMC processing), and
42 `R/sampler_*.R` files, one per MCMC algorithm. `missing(Specs)` is
replaced with `is.null(Specs)` throughout extracted functions.

### Multi-chain support

[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
accepts `Chains`, `CPUs`, and `Type` arguments directly. When
`Chains > 1`, chains run in parallel when `CPUs > 1`, are automatically
combined via
[`Combine()`](https://robustecologies.github.io/lucifer/reference/Combine.md),
and the Gelman-Rubin diagnostic is computed. The `demonoid` object gains
`$Chains`, `$Chains.n`, `$Gelman`, and `$CPUs` fields. Single-chain
output is unchanged for backward compatibility. `lucifer.hpc()` is
deprecated in favour of `lucifer(..., Chains = n, CPUs = n)`.

### ggplot2 visualisation

All 21 plot method files have been rewritten from base R to ggplot2 with
a consistent RElab theme
([`theme_relab()`](https://robustecologies.github.io/lucifer/reference/theme_relab.md)
in `R/theme_relab.R`). The theme provides a 9-colour palette, custom
scales
([`scale_color_relab_d()`](https://robustecologies.github.io/lucifer/reference/scale_color_relab_d.md),
[`scale_fill_relab_c()`](https://robustecologies.github.io/lucifer/reference/scale_fill_relab_c.md),
etc.), and publication-quality defaults. Five PPC plot files (~680 lines
each, highly duplicated) are consolidated into a single shared
implementation (`R/plot_ppc_internal.R`, ~1,100 lines, 26 plot styles)
with thin S3 wrappers for `plot.demonoid.ppc`, `plot.iterquad.ppc`,
`plot.laplace.ppc`, `plot.pmc.ppc`, and `plot.vb.ppc`.

### Documentation and utility infrastructure

All R files are converted from hand-written `.Rd` to roxygen2. The
NAMESPACE is regenerated from roxygen2 with 395 exports, 35 S3 methods,
and 81 imports. The new `R/utils.R` provides `get_colored_symbol()` for
Unicode symbols with ANSI colour output,
`require_package()`/`require_packages()` for dependency checking, and
the `ProgressBar` R6 class for iterative computation progress.
[`summary.demonoid()`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md)
is added as a new S3 method reporting per-chain acceptance rates,
Gelman-Rubin convergence, and combined posterior summaries.
[`ld_cpp_available()`](https://robustecologies.github.io/lucifer/reference/ld_cpp_available.md)
checks whether the C++ backend is operational.
