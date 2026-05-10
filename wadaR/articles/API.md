# Technical architecture of the lucifer inference framework

  

## Design philosophy

lucifer is a unified Bayesian inference framework that consolidates more
than 130 algorithms, spanning MCMC, variational inference, quadrature,
importance sampling, sequential Monte Carlo, approximate Bayesian
computation, simulation-based inference, and several domain-specific
engines, behind a single model interface. The user writes one model
function and gains access to every inference strategy without rewriting
anything; switching from NUTS to ADVI to Laplace approximation requires
changing a single argument, not restructuring the model. This section
lays out the three architectural principles that make this possible.

### Unified model contract

Every inference engine in lucifer consumes the same function signature:
`Model(parm, Data)`, where `parm` is a numeric vector of parameters and
`Data` is a named list containing observations, hyperparameters, and
metadata. The function returns a five-element list containing the log of
the unnormalized joint posterior (`LP`), the deviance (`Dev`),
user-tracked quantities (`Monitor`), fitted values (`yhat`), and the
(potentially constrained) parameter vector (`parm`). This contract is
the single integration point between model definition and inference
machinery; the model function is entirely agnostic to which engine will
process it. The consequence is that a model written for a quick
variational approximation can be handed unchanged to a full MCMC
sampler, a quadrature engine, or the automated Crucible pipeline without
any modification.

### Progressive refinement

lucifer is designed for iterative workflow: start with a fast
approximate method, diagnose the result, and escalate to a more
expensive sampler only when necessary. The Consort diagnostic advisor
evaluates any fit object and returns structured suggestions for
improvement, including specific algorithm switches and tuning parameter
adjustments. The Crucible orchestrator automates this entire loop,
profiling the model via Prescribe, fitting multiple candidate methods,
diagnosing each with Consort, refining those that fail convergence
criteria, and ranking the survivors via Arena. Every output class feeds
into this pipeline because they all share a common diagnostic interface.

### Zero-cost abstraction

The R-level API remains simple: a handful of user-facing functions with
clear signatures and sensible defaults. Computational intensity is
pushed to C++ via Rcpp, where roughly 35,000 lines of optimized code
handle distribution evaluations, MCMC kernels, diagnostics, PSIS
importance weights, particle filters, Kalman filters, and more, all with
OpenMP parallelization where appropriate. The
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
mechanism takes this further by translating a declarative `model_spec`
object into native C++ function pointers, eliminating R interpreter
callbacks entirely during sampling. This means that gradient-based
samplers like NUTS evaluate both the log-posterior and its gradient in
pure C++, with no per-iteration R overhead.

Layered architecture of the lucifer inference framework.

  

## Architecture overview

The package is organized around four subsystems. The **inference
engines** (lucifer, VariationalBayes, LaplaceApproximation,
IterativeQuadrature, PMC, SMC, ABC, SBI, and the domain-specific NODE,
SDE, SSM) are the primary user-facing functions that accept a model and
produce posterior samples or approximations. The **model specification**
system (model_spec) provides a declarative DSL that compiles
probabilistic notation into executable Model functions, optionally with
native C++ backends. The **orchestration layer** (Prescribe, Consort,
Arena, Crucible) automates algorithm selection, convergence diagnosis,
iterative refinement, and cross-method comparison. The
**infrastructure** (algo_registry, C++ backend, parallel_config)
provides the shared machinery that all other subsystems depend on.

### The model contract

The five-element return list from `Model(parm, Data)` is the universal
interface that every inference engine expects:

- **LP:** the log of the unnormalized joint posterior, \\\log p(\theta
  \| y) + C\\. This is the quantity that MCMC samplers target and
  variational methods optimize.
- **Dev:** the deviance, \\-2 \log p(y \| \theta)\\. Used for DIC, WAIC,
  and LOO-PSIS computations.
- **Monitor:** a named numeric vector of user-tracked quantities
  computed at each iteration (e.g., derived parameters, predictions,
  likelihood components).
- **yhat:** fitted or predicted values at the current parameter values.
  These feed posterior predictive checks (PPC) and residual diagnostics.
- **parm:** the parameter vector, potentially modified by constraint
  transformations via
  [`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md).
  Returning the constrained version ensures the sampler operates on the
  correct support.

The
[`is.model()`](https://robustecologies.github.io/lucifer/reference/is.model.md)
function validates this contract by calling the Model function once and
checking the return structure. Any function that satisfies this contract
works with every engine in the package.

``` r

# Minimal model function satisfying the contract
MyModel <- function(parm, Data) {
    # Extract and constrain parameters
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    beta <- parm[Data$pos.beta]

    # Log-likelihood
    mu <- Data$X %*% beta
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))

    # Log-prior
    LP <- LL + sum(dnorm(beta, 0, 100, log = TRUE)) +
        dhalfcauchy(sigma, 25, log = TRUE)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(length(mu), mu, sigma), parm = parm)
}
```

Component relationship map showing data flow between major subsystems.

### Inference engine summary

| Engine                 | Methods | Gradient? | Discrete? | C++ backend? | Output class |
|:-----------------------|:-------:|:---------:|:---------:|:------------:|:-------------|
| lucifer()              |   82    |  varies   |  varies   |     yes      | demonoid     |
| VariationalBayes()     |    8    |   yes\*   |    no     |   partial    | vb           |
| LaplaceApproximation() |   18    |  varies   |    no     |      no      | laplace      |
| IterativeQuadrature()  |    3    |    no     |    no     |      no      | iterquad     |
| PMC()                  |    1    |    no     |    no     |      no      | pmc          |
| SMC()                  |    1    |    no     |    no     |     yes      | smc          |
| ABC()                  |    4    |    no     |    no     |   partial    | abc          |
| SBI()                  |    6    |    no     |    no     |     yes      | sbi          |
| NODE()                 |    1    |    no     |    no     |     yes      | node_fit     |
| SDE()                  |   11    |    no     |    no     |     yes      | sde_fit      |
| SSM()                  |    5    |    no     |    no     |     yes      | ssm_fit      |

Inference engine summary. Methods = number of algorithms or model
families. \*VB gradient requirement depends on method. {.table
style="width:100%;"}

  

## Model specification

lucifer supports two approaches to model definition: hand-written Model
functions and the declarative `model_spec` DSL. Both produce the same
five-element contract described above; they differ in how the user
expresses the model.

### Hand-written model functions

The traditional approach is to write a plain R function that manually
extracts parameters from the `parm` vector by position (using
`Data$pos.*` indices), applies constraints via
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md),
computes the log-likelihood and log-prior, and assembles the return
list. This approach offers maximum flexibility: the user controls every
computation, can embed arbitrary R logic, and can call external
functions. The cost is verbosity and the possibility of subtle errors in
parameter indexing.

``` r

Model <- function(parm, Data) {
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma

    mu <- Data$X %*% beta
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + sum(dnorm(beta, 0, 1000, log = TRUE)) +
        dhalfcauchy(sigma, 25, log = TRUE)

    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(length(mu), mu, sigma), parm = parm)
}
```

### The model_spec DSL

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
function accepts a string of probabilistic notation where stochastic
relationships use `~` and deterministic relationships use `=`. The
parser handles LaTeX-like syntax (Greek letters, subscripts, `\sim`),
resolves distribution names to canonical forms from lucifer’s library of
86+ distributions, and classifies every variable as data, parameter,
hyperparameter, or deterministic node. The result is a `model_spec` S3
object containing a compiled `Model` function, a `data_template()`
builder for constructing the Data list, an `initial_values()` generator,
the intermediate representation (IR), and the generated source code.

``` r

spec <- model_spec("
    y[i] ~ Normal(mu[i], sigma)
    mu[i] = X[i,] %*% beta
    beta[j] ~ Normal(0, 100)
    sigma ~ HalfCauchy(25)
")

# The spec object contains everything needed for inference
fit <- lucifer(spec$Model, spec$data_template(y = y, X = X),
               spec$initial_values(), Algorithm = "NUTS",
               Iterations = 5000)
```

### Compilation pipeline

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
function processes the input through three stages. First, the **parser**
(`.parse_model_notation()`) tokenizes the notation string, strips LaTeX
formatting, identifies index ranges, and resolves distribution names
against the internal registry. Second, the **IR builder**
(`.build_ir()`) classifies variables by their role: variables that
appear on the left-hand side of the first stochastic statement and are
not referenced by other priors are classified as data (likelihood);
other stochastic left-hand sides become parameters; right-hand-side-only
variables become data requirements. The IR also performs topological
sorting of deterministic dependencies and infers constraint types from
distribution supports (e.g., HalfCauchy implies a positive constraint).
Third, the **code generator** (`.compile_ir()`) emits a complete Model
function in R source code, including parameter extraction, constraint
application, log-likelihood accumulation, log-prior accumulation, and
the five-element return list.

### Native compilation

When `backend = "cpp"` is specified (or when
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
is called on an existing model_spec object), a parallel code generation
path activates. The C++ code generator (`.compile_ir_cpp()`) translates
the IR into C++ source that uses lucifer’s header-only gradient library
(`inst/include/lucifer/dist_gradients.h`) for analytical log-density
gradients across 24 distributions, and the transform library
(`inst/include/lucifer/transforms.h`) for constraint transformations
with correct Jacobian adjustments. The generated source is compiled at
runtime via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) with
caching, producing a `compiled_model` S3 object that exposes native
function pointers. When a compiled model is passed to a sampler, the C++
dispatch layer (`call_model()` in `sampler_common.cpp`) invokes the
function pointer directly, bypassing R entirely. The dual-number AD
library (`inst/include/lucifer/dual.h`) provides forward-mode automatic
differentiation as an alternative to analytical gradients for
distributions not yet in the gradient library.

Model specification pipeline with R and C++ compilation backends.

  

## Inference engines

This section documents each of the eight major inference engine families
and the five domain-specific engines. Each engine accepts a Model
function satisfying the five-element contract, applies its inference
strategy, and returns an S3 object with print, summary, and plot
methods.

### MCMC: `lucifer()`

The
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function is the main MCMC dispatcher. It accepts an `Algorithm`
parameter (a string abbreviation), validates algorithm-specific tuning
parameters via `Specs`, routes to the appropriate sampler implementation
(one of ~55 R wrapper files that typically call C++ backends), and
postprocesses the raw chain into a `demonoid` S3 object. The 82 MCMC
algorithms span 14 subcategories covering gradient-based dynamics,
adaptive random walks, componentwise updates, slice sampling, ensemble
methods, quasi-Monte Carlo, Riemannian geometry, piecewise deterministic
processes, tempering, constraint handling, normalizing flows, surrogate
models, and data augmentation.

``` r

fit <- lucifer(
    Model, Data, Initial.Values,
    Algorithm = "NUTS",
    Specs = list(epsilon = 0.1, L = 20),
    Iterations = 10000,
    Thinning = 10,
    Chains = 4,
    CPUs = 4
)
```

When `Chains > 1` and `CPUs > 1`, each chain runs in a separate R
subprocess via
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html) (falling
back to a PSOCK cluster if callr is unavailable), and the results are
merged via
[`Combine()`](https://robustecologies.github.io/lucifer/reference/Combine.md)
with automatic Rhat computation. The `Specs` parameter accepts
algorithm-specific tuning; when omitted, sensible defaults are applied.
For models with user-supplied gradients (`Data$user.grad` or
`Data$gradient`), gradient-based algorithms use the analytical gradient
directly; otherwise they fall back to finite differences via
`partial_cpp()` in C++.

The auto-FC Gibbs mechanism deserves special mention: when
`Algorithm = "Gibbs"` and `Specs$FC = "auto"`, lucifer constructs
full-conditional samplers automatically using stepping-out slice
sampling for continuous parameters and exhaustive enumeration for
discrete parameters, enabling Gibbs sampling without manually specifying
full conditionals.

Decision flowchart for selecting an inference engine. Diamond nodes are
decision points; rectangular nodes are recommended engines.

#### MCMC algorithm families

The 82 MCMC algorithms are organized into 14 subcategories. The
following table summarizes each family; the complete algorithm list
appears in the quick reference at the end of this document.

| Family | Count | Gradient? | Discrete? | Key property | Representative |
|:---|:--:|:--:|:--:|:---|:---|
| Gradient-based | 16 | yes | no | Hamiltonian dynamics, Langevin diffusion | NUTS, GIST, autoMALA |
| Adaptive random-walk | 12 | no | no | Learned proposal covariance | RAM, DRAM |
| Componentwise | 13 | no | yes\* | Parameter-by-parameter updates | AMWG, Gibbs, Zanella |
| Slice | 10 | varies | no | Rejection-free, tuning-free | AFSS, ESS, QSS |
| Ensemble | 6 | varies | no | Affine invariance, population-based | AIES, Zeus, DREAM |
| QMC | 8 | no | no | Quadratic/higher-order proposals | QMC, DQMC |
| Geometric | 3 | yes | no | Riemannian metric tensor adaptation | RMHMC |
| PDMP | 4 | yes | no | Continuous-time, event-driven | BPS, ZigZag |
| Multimodal | 4 | no | no | Tempering, flat histogram | NRPT, NRST, SimTemp |
| Constraint | 2 | yes | no | Projection, proximal operators | ProjLang, ProxMCMC |
| Flow-enhanced | 2 | varies | no | Normalizing flow transport | NeuTra, flowMC |
| Surrogate | 1 | no | no | Two-stage acceptance | DA |
| Other | 2 | no | yes\* | Reversible-jump, pseudo-marginal | RJ, PMCMC |
| Augmentation | 1 | no | no | Polya-Gamma data augmentation | PG |

MCMC algorithm families. \*Discrete support limited to specific
algorithms within the family (Gibbs, AMWG, Zanella, RJ). {.table}

### Variational Bayes: `VariationalBayes()`

[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
provides eight deterministic variational inference methods spanning
fixed-form (Salimans1, Salimans2), automatic differentiation variational
inference (ADVI mean-field and full-rank), black-box variational
inference (BBVI), Stein variational gradient descent (SVGD), natural
gradient descent (NGD), and Pathfinder (quasi-Newton trajectory-based).
All methods optimize an approximation to the posterior rather than
sampling from it, trading exactness for speed. The output is a `vb` S3
object containing the approximate posterior mean, covariance, and
optionally importance-resampled posterior draws.

Pathfinder stands apart from the other VB methods: it traces the L-BFGS
optimization trajectory, fits a multivariate normal approximation at
each point along the path, selects the best by ELBO, and optionally runs
multiple independent paths in parallel (via PSOCK workers) with PSIS
resampling to combine draws across paths. Its C++ core
(`src/pathfinder.cpp`) handles the L-BFGS trajectory and L-BFGS inverse
Hessian approximation. For most problems, Pathfinder is the fastest
route to a reasonable posterior approximation.

### Laplace approximation: `LaplaceApproximation()`

[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
finds the posterior mode and approximates the posterior as a
multivariate normal centered at the mode with covariance derived from
the inverse Hessian of the negative log-posterior. It offers 18
optimization methods across five subcategories: quasi-Newton (BFGS,
L-BFGS, DFP, BHHH, SR1), gradient-based (AGA, CG, Rprop, SGD, SPG),
second-order (NR, LM, TR), derivative-free (NM, HJ, HAR), and
metaheuristic (PSO, SOMA). The output is a `laplace` S3 object. This
engine is the fastest path to a point estimate with uncertainty
quantification when the posterior is approximately Gaussian.

### Iterative quadrature: `IterativeQuadrature()`

[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)
provides Gaussian quadrature for low-dimensional problems (typically \\d
\< 10\\). Three methods are available: adaptive Gauss-Hermite (AGH),
adaptive Gauss-Hermite with sparse grids (AGHSG), and componentwise
adaptive Gauss-Hermite (CAGH). The engine iteratively refines the
quadrature grid location and scale to center on the posterior mode,
providing numerically exact posterior moments for smooth,
low-dimensional posteriors. The output is an `iterquad` S3 object. The
\\d \< 10\\ limitation arises from the exponential growth of quadrature
points with dimension; for higher-dimensional problems, MCMC or VB
should be used instead.

### Population Monte Carlo: `PMC()`

[`PMC()`](https://robustecologies.github.io/lucifer/reference/PMC.md)
implements adaptive importance sampling with iterative mixture fitting.
At each iteration, importance weights are computed for the current
particle set, the mixture proposal is updated to match the reweighted
particle distribution, and new particles are drawn from the updated
proposal. The algorithm converges when the effective sample size
stabilizes. The output is a `pmc` S3 object. PMC is particularly
effective for moderately multimodal posteriors where the mixture
proposal can cover multiple modes.

### Sequential Monte Carlo: `SMC()`

[`SMC()`](https://robustecologies.github.io/lucifer/reference/SMC.md)
implements sequential Monte Carlo with adaptive tempering, progressively
bridging from the prior to the posterior through a sequence of tempered
distributions \\p(\theta)^{1-\beta_t} \cdot p(y\|\theta)^{\beta_t}\\
with \\\beta_t\\ increasing from 0 to 1. The tempering schedule is
determined adaptively by targeting a minimum ESS ratio at each stage.
Particle rejuvenation uses random-walk Metropolis moves; for compiled
models, the inner particle loop is parallelized via conditional OpenMP.
The C++ backend resides in `src/smc.cpp`. The output is an `smc` S3
object.

### Approximate Bayesian computation: `ABC()`

[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
provides likelihood-free inference via four methods: rejection sampling,
MCMC-ABC, SMC-ABC, and simulated annealing ABC (SA-ABC). These methods
require a simulator function rather than a likelihood, comparing
simulated data to observed data through user-specified summary
statistics and a distance metric. When `Method = "rejection"` and
`CPUs > 1`, the target sample count is split across PSOCK workers for
parallel execution. MCMC-ABC, SMC-ABC, and SA-ABC remain sequential
because their state is not decomposable. The output is an `abc` S3
object.

### Simulation-based inference: `SBI()`

[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
implements neural density, likelihood, and ratio estimation for
amortized or sequential inference. Six methods are available: NPE
(neural posterior estimation), NLE (neural likelihood estimation), NRE
(neural ratio estimation), and their sequential variants SNPE, SNLE, and
SNRE, plus TSNPE (truncated sequential NPE). The C++ backends handle
mixture density networks (`src/mdn.cpp`), neural ratio classifiers
(`src/nre_classifier.cpp`), and parallel forward simulation
(`src/sbi_simulate.cpp`). SBI requires torch and is designed for
problems where the likelihood is intractable but a fast simulator is
available. The output is an `sbi` S3 object. Diagnostic tools include
C2ST (classifier two-sample test), TARP (tests of accuracy with random
points), and SBI-specific posterior predictive checks.

  

## Algorithm registry

The algorithm registry (`R/algo_registry.R`) is the centralized metadata
dictionary that powers the orchestration layer. Every algorithm in the
package, across all eight inference categories plus SSM, has an entry
built by `.make_algo_entry()` with 16 standardized fields. The combined
registry `.algo_registry` is a named list of 122 entries constructed by
concatenating the category-specific sub-registries (`.mcmc_registry`,
`.vb_registry`, `.laplace_registry`, `.iq_registry`, `.pmc_registry`,
`.smc_registry`, `.abc_registry`, `.sbi_registry`,
`.ssm_algo_registry`).

Lookup is performed by abbreviation via `.algo_lookup()` or through the
exported
[`algo_info()`](https://robustecologies.github.io/lucifer/reference/algo_info.md)
function, which returns the full metadata for any algorithm by name or
abbreviation. Category-level queries use `.algo_by_category()`. The
registry drives automatic decision-making throughout the orchestration
layer: Prescribe scores algorithms based on gradient availability,
dimensionality, and problem characteristics drawn from registry fields;
Consort uses `acceptance_range` to diagnose acceptance rate failures and
`subcategory` to determine escalation paths; Crucible uses
`quality_tier` to prioritize high-quality algorithms when selecting
candidates.

| Field | Type | Description |
|:---|:---|:---|
| abbrev | character | Short algorithm name used in Algorithm parameter |
| category | character | MCMC, VB, Laplace, IQ, PMC, SMC, ABC, SBI, or SSM |
| subcategory | character | Fine-grained family (gradient, adaptive_rw, slice, …) |
| requires_gradient | logical | Whether the algorithm requires gradient evaluations |
| supports_discrete | logical | Whether the algorithm handles discrete parameters |
| componentwise | logical | Whether the algorithm updates parameters one at a time |
| acceptance_range | numeric\[2\] | Ideal acceptance rate window \[low, high\], or NULL |
| has_adaptive_phase | logical | Whether the algorithm has a warmup/adaptation phase |
| sde_penalty | numeric | Cost modifier for SDE model fitting (0 = neutral) |
| full_name | character | Full human-readable algorithm name |
| multimodal_affinity | numeric | Affinity score 0-1 for multimodal posteriors |
| constraint_affinity | numeric | Affinity score 0-1 for constrained parameter spaces |
| requires_torch | logical | Whether the algorithm requires the torch package |
| eval_cost_multiplier | numeric | Relative cost per iteration vs. baseline |
| dim_range | integer\[2\] | Min/max dimensionality \[low, high\], or NULL |
| quality_tier | integer | Quality tier: 1 (recommended), 2 (standard), 3 (legacy) |

Algorithm registry fields. Each algorithm entry contains all 16 fields,
enabling automated scoring and recommendation. {.table}

  

## Orchestration

The orchestration layer provides four composable tools that automate the
inference workflow. Prescribe profiles a model without fitting it and
ranks algorithms by expected performance. Consort evaluates any fit
object against convergence criteria and generates structured improvement
suggestions. Arena compares multiple fits from different methods using
standardized metrics. Crucible composes all three into a fully automated
pipeline. Each tool is designed to be useful independently; Crucible
simply chains them together.

### Prescribe: pre-fit profiling

[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
characterizes a model’s computational properties without performing full
inference. It calls the Model function approximately \\n + 2d\\ times
(where \\n\\ is the profile sample count and \\d\\ is the
dimensionality) to measure evaluation speed, check gradient
availability, assess constraint structure, estimate posterior
conditioning via local Hessian analysis, and detect potential
multimodality through kurtosis and range-to-SD ratios of log-posterior
samples. Each algorithm in the registry then receives a composite score
computed as the product of multiplicative factors; a zero factor
disqualifies the algorithm entirely (e.g., a gradient-based algorithm
receives a zero factor when no gradient is available). The output is a
`prescription` S3 object containing the model profile, a ranked score
table, and a primary recommendation with ready-to-paste code.

``` r

rx <- Prescribe(Model, Data, Initial.Values)
print(rx)   # top 5 recommendations with code snippets
plot(rx)     # score comparison across categories
```

### Consort: post-fit diagnostics

[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
dispatches to family-specific diagnostic engines based on the S3 class
of the fit object. For MCMC fits (class `demonoid`), it evaluates eight
conditions: the chain must not be in its adaptive phase; the acceptance
rate must fall within the algorithm-specific range stored in the
registry; the maximum MCSE-to-SD ratio must be below 6.27%; bulk ESS
must exceed 400 and tail ESS must exceed 200 across all parameters;
split Rhat must be below 1.01; there must be zero divergent transitions;
and adaptation diagnostics must show diminishing adaptation. For VB
fits, Consort checks convergence, ELBO stability (coefficient of
variation below 0.05 over the tail 10% of iterations), and SIR Pareto k
below 0.7. For Laplace fits, it checks optimizer convergence, Hessian
positive-definiteness, and gradient norm below \\10^{-3}\\.

When conditions fail, Consort generates a structured suggestion
containing the recommended action (algorithm switch, iteration increase,
or tuning adjustment), the specific `Specs` to use in the next call,
ready-to-paste code, and a rationale explaining the escalation logic.
The suggestion priority order is: divergences first (switch to NUTS or
AFSS), then diminishing adaptation failures (escalation map with 31
algorithm-to-algorithm transitions), then acceptance rate issues, and
finally ESS/MCSE/Rhat shortfalls (which typically indicate insufficient
iterations).

A fit is “appeased” when all conditions pass. The
[`is.appeased()`](https://robustecologies.github.io/lucifer/reference/is.appeased.md)
function provides a quick Boolean check.

Crucible orchestration flow. The Consort feedback loop refines
non-appeased fits up to max_rounds before passing all results to Arena.

Consort MCMC diagnostic decision tree. Conditions are checked in
priority order; the first failure determines the suggestion.

### Arena: cross-method comparison

[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
accepts a named list of fit objects from different methods and computes
standardized metrics for comparison. Each fit is processed through a
class-specific extraction adapter (`.arena_extract_demonoid()`,
`.arena_extract_vb()`, etc.) that normalizes the output into a common
structure: posterior mean, posterior SD, posterior samples, minimum ESS,
ESS per second, wall-clock time, convergence status, and log marginal
likelihood estimate. Metrics include efficiency (ESS per second,
bottlenecked by the minimum across parameters), accuracy (marginal KLD
and Wasserstein-1 distance to a reference distribution), and reliability
(convergence diagnostics specific to each method). A consensus reference
is constructed by ESS-weighted averaging across all posteriors when no
external reference is provided. Arena identifies the Pareto frontier,
methods not dominated in the time-accuracy tradeoff, and computes a
composite ranking.

### Crucible: automated pipeline

[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
composes Prescribe, the inference engines, Consort, and Arena into a
six-stage automated pipeline. Stage 1 profiles the model via Prescribe.
Stage 2 selects the top N methods from the ranked list, enforcing
diversity across families when `diverse = TRUE`. Stages 3-5 form a
per-method loop: fit the method, evaluate with Consort, and if not
appeased, extract the structured suggestion and re-fit with the
recommended adjustments (warm-starting from the last posterior draw for
MCMC), repeating up to `max_rounds`. Stage 6 passes all final fits to
Arena for cross-method comparison. The output is a `crucible` S3 object
containing the Prescribe profile, per-method round histories, the Arena
comparison, and the overall best fit.

``` r

result <- Crucible(
    Model, Data, Initial.Values,
    n_methods = 5, max_rounds = 3,
    Chains = 4, CPUs = 4,
    families = c("MCMC", "VB", "Laplace")
)

print(result)              # summary table with rankings
plot(result)               # convergence traces
plot(result, type = "arena")  # efficiency vs accuracy
```

  

## S3 class system

Every inference engine returns an S3 object, and every S3 class
implements the triad of [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. This
uniformity means that the user interacts with the results of any engine
through the same interface:
[`print()`](https://rdrr.io/r/base/print.html) produces a concise
one-screen synopsis showing the algorithm used, key parameters, and
primary diagnostics; [`summary()`](https://rdrr.io/r/base/summary.html)
provides extended output with convergence metrics, parameter tables, and
interpretation guidance;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) generates the
most informative single visualization by default and supports a `type`
argument for alternative views (e.g., `type = "convergence"`,
`type = "diagnostics"`, `type = "autocorrelation"`).

Posterior predictive checks are available for all fitting classes via
[`as.ppc()`](https://robustecologies.github.io/lucifer/reference/as.ppc.md),
which extracts fitted values and observed data into a PPC object with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods. The PPC
framework is implemented once in `R/plot_ppc_internal.R` and shared
across all families; the class-specific `plot.*.ppc` and `summary.*.ppc`
methods dispatch to this shared implementation. The
[`predict()`](https://rdrr.io/r/stats/predict.html) generic is
implemented for demonoid, vb, laplace, iterquad, pmc, and sbi objects,
generating posterior predictive samples.

The
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function dispatches on the S3 class of the fit object to apply
family-specific diagnostic engines: demonoid fits use MCMC-specific
checks (acceptance rate, ESS, Rhat, divergences); vb fits use
convergence and ELBO stability checks; laplace fits use Hessian
conditioning checks; and so on. This class-based dispatch is what allows
Crucible to work transparently across all inference strategies.

| Class        | print | summary | plot | predict | Consort | PPC | LOO/WAIC |
|:-------------|:-----:|:-------:|:----:|:-------:|:-------:|:---:|:--------:|
| demonoid     |  yes  |   yes   | yes  |   yes   |   yes   | yes |   yes    |
| vb           |  yes  |   yes   | yes  |   yes   |   yes   | yes |   yes    |
| laplace      |  yes  |   yes   | yes  |   yes   |   yes   | yes |   yes    |
| iterquad     |  yes  |   yes   | yes  |   yes   |   yes   | yes |   yes    |
| pmc          |  yes  |   yes   | yes  |   yes   |   yes   | yes |   yes    |
| smc          |  yes  |   yes   | yes  |   no    |   yes   | yes |    no    |
| abc          |  yes  |   yes   | yes  |   no    |   yes   | no  |    no    |
| sbi          |  yes  |   yes   | yes  |   yes   |   yes   | yes |    no    |
| node_fit     |  yes  |   yes   | yes  |   yes   |   no    | no  |    no    |
| sde_fit      |  yes  |   yes   | yes  |   yes   |   yes   | no  |    no    |
| ssm_fit      |  yes  |   yes   | yes  |   no    |   no    | no  |    no    |
| crucible     |  yes  |   yes   | yes  |   no    |   no    | no  |    no    |
| prescription |  yes  |   yes   | yes  |   no    |   no    | no  |    no    |
| consort      |  yes  |   yes   | yes  |   no    |   no    | no  |    no    |
| arena        |  yes  |   yes   | yes  |   no    |   no    | no  |    no    |

S3 class method availability matrix. {.table style="width:100%;"}

  

### S3 family index

lucifer defines 63 S3 classes with 173 registered methods. The index
below groups classes by semantic category and lists every generic
currently dispatched for each class. Every row is also reachable from
the constructor help page through the `@seealso` cross-reference graph,
which is enforced bidirectionally in the roxygen sources.

  

#### Core Bayesian engines

| S3 class | Generics dispatched |
|:---|:---|
| abc.ppc | plot |
| demonoid | log_lik, LOO, plot, predict, print, summary, to_draws_array, to_draws_df, to_draws_matrix, to_mcmc_list, WAIC |
| demonoid.ppc | plot, summary |
| importance | plot |
| iterquad | plot, predict, print |
| iterquad.ppc | plot, summary |
| laplace | plot, predict, print |
| laplace.ppc | plot, summary |
| miss | plot, print, summary |
| pmc | log_lik, LOO, plot, predict, print, to_mcmc_list, WAIC |
| pmc.ppc | plot, summary |
| smc.ppc | plot |
| vb | log_lik, LOO, plot, predict, print, summary, to_mcmc_list, WAIC |
| vb.ppc | plot, summary |

Core Bayesian engines (14 classes). {.table}

  

#### Advanced inference

| S3 class  | Generics dispatched                |
|:----------|:-----------------------------------|
| abc       | plot, print, summary               |
| bayesquad | plot, predict, print, summary      |
| node_fit  | plot, predict, print, summary      |
| sbi       | plot, print, summary               |
| smc       | plot, print, summary, to_mcmc_list |

Advanced inference (5 classes). {.table}

  

#### Domain-specific engines

| S3 class         | Generics dispatched                               |
|:-----------------|:--------------------------------------------------|
| bayes_mode       | plot, print, summary                              |
| bayes_mode_multi | plot, print, summary                              |
| compiled_model   | code, print                                       |
| mode_estimate    | as.double, plot, print, summary                   |
| model_spec       | code, compile, plot, print, summary               |
| sde_fit          | log_lik, LOO, plot, predict, print, summary, WAIC |
| sde_model        | compile, plot, print, simulate, summary           |
| sde_prediction   | plot, print                                       |
| ssm_fit          | as.demonoid, plot, print, summary                 |
| ssm_model        | print                                             |

Domain-specific engines (10 classes). {.table}

  

#### Diagnostics

| S3 class            | Generics dispatched |
|:--------------------|:--------------------|
| bmk                 | plot                |
| coupling_diagnostic | plot, print         |
| heidelberger        | print               |
| nuts_check          | plot, print         |
| raftery             | print               |

Diagnostics (5 classes). {.table}

  

#### Orchestration and diagnostics

| S3 class     | Generics dispatched  |
|:-------------|:---------------------|
| arena        | plot, print, summary |
| consort      | plot, print, summary |
| coupled_mcmc | plot, print, summary |
| coverage_sim | plot, print, summary |
| crucible     | plot, print, summary |
| prescription | plot, print, summary |

Orchestration and diagnostics (6 classes). {.table}

  

#### Frequentist bridge

| S3 class           | Generics dispatched  |
|:-------------------|:---------------------|
| freq_residuals     | plot, print          |
| freq_summary       | plot, print          |
| lr_test            | print                |
| profile_likelihood | plot, print, summary |
| score_test         | print                |
| wald_test          | print                |

Frequentist bridge (6 classes). {.table}

  

#### Model comparison and robustness

| S3 class         | Generics dispatched  |
|:-----------------|:---------------------|
| data_cloning     | plot, print, summary |
| kfold            | plot, print, summary |
| lfo              | plot, print, summary |
| loo              | plot, print, summary |
| loo_comparison   | plot, print          |
| product_space    | plot, print, summary |
| robust_bayes     | plot, print, summary |
| ssvs_summary     | plot, print, summary |
| stacking_weights | print                |
| waic             | plot, print          |

Model comparison and robustness (10 classes). {.table}

  

#### Interop

| S3 class     | Generics dispatched    |
|:-------------|:-----------------------|
| default      | as.demonoid, LOO, WAIC |
| draws_array  | as.demonoid            |
| draws_df     | as.demonoid            |
| draws_matrix | as.demonoid            |
| mcmc         | as.demonoid            |
| mcmc.list    | as.demonoid            |

Interop (6 classes). {.table}

  

#### Other

| S3 class      | Generics dispatched |
|:--------------|:--------------------|
| node_forecast | plot, print         |

Other (1 class). {.table}

  

## Burn-in and posterior splitting

Every call to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
produces two posterior matrices. `Posterior1` contains all thinned
samples, including the transient (burn-in) phase where the chain has not
yet reached its stationary distribution. `Posterior2` contains only the
stationary samples, with burn-in removed. The fields `Summary2`, `DIC2`,
`LML`, and [`predict()`](https://rdrr.io/r/stats/predict.html) all
operate on `Posterior2`; `Summary1` and `DIC1` operate on `Posterior1`.

### Automatic burn-in detection

By default (`BurnIn = NULL`), the burn-in boundary is estimated using
the BMK stationarity diagnostic, which computes the Hellinger distance
between successive 10% batches of the chain. The first batch at which
all subsequent batches are stationary (HD \\\leq\\ 0.5) defines the end
of burn-in. If no stationary region is found, the first 50% of samples
are discarded as a conservative fallback, and a warning is issued.

### User-specified burn-in

The `BurnIn` argument accepts a non-negative integer indicating how many
thinned samples to discard from the start. When provided, the automatic
BMK diagnostic is bypassed entirely:

``` r

# Discard first 500 thinned samples as burn-in
fit <- lucifer(Model, Data, IV, Iterations = 10000, Thinning = 10,
               Algorithm = "NUTS", BurnIn = 500)
```

The `BurnIn.Method` field in the returned object records whether burn-in
was `"user"` (user-specified) or `"BMK"` (automatic). This is shown in
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html).

### Post-hoc adjustment

After fitting,
[`deburn()`](https://robustecologies.github.io/lucifer/reference/deburn.md)
allows further burn-in removal without re-running the sampler. It strips
the specified number of thinned samples from `Posterior1` and sets
`Posterior2 = Posterior1` (no further splitting). The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method also
accepts a `BurnIn` argument for visualization-only filtering that does
not modify the object.

  

## Diagnostics

lucifer provides three tiers of convergence diagnostics: modern
diagnostics that form the basis of Consort’s automated checks, classical
diagnostics retained for compatibility and exploratory use, and
algorithm-specific diagnostics for specialized sampler types.

### Modern convergence diagnostics

The modern diagnostic suite implements the recommendations of Vehtari et
al. (2021) and is the default tier used by Consort.
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
computes the split \\\hat{R}\\ with rank normalization, which handles
multimodal and heavy-tailed distributions more reliably than the
classical Gelman-Rubin statistic.
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md)
returns the bulk effective sample size (measuring mixing in the center
of the distribution) and
[`ESS.tail()`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md)
returns the tail ESS (measuring mixing in the extremes).
[`MCSE()`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
computes the Monte Carlo standard error;
[`MCSE.quantile()`](https://robustecologies.github.io/lucifer/reference/MCSE.quantile.md)
extends this to arbitrary quantiles. All modern diagnostics have C++
backends in `src/diagnostics.cpp` with OpenMP parallelization across
parameters.

Consort uses the following thresholds: split \\\hat{R} \< 1.01\\
(indicating convergence), bulk ESS \\\geq 400\\ (sufficient mixing for
reliable mean estimation), tail ESS \\\geq 200\\ (sufficient mixing for
reliable quantile estimation), and MCSE/SD \\\< 6.27\\\\ (Monte Carlo
error is a small fraction of posterior uncertainty). These thresholds
are applied per-parameter, and the worst-case parameter determines the
overall diagnostic status.

### Classical diagnostics

Six classical diagnostics remain available:
[`Geweke.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/Geweke.Diagnostic.md)
compares the first 10% and last 50% of the chain via a Z-test for
difference in means;
[`Heidelberger.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/Heidelberger.Diagnostic.md)
tests for stationarity and estimates the half-width of a confidence
interval for the mean;
[`Raftery.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/Raftery.Diagnostic.md)
estimates the dependence factor and required chain length for a target
quantile accuracy;
[`Gelman.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md)
computes the classical multi-chain \\\hat{R}\\ without rank
normalization;
[`BMK.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md)
applies the batch means kernel test for single-chain stationarity; and
[`KS.Diagnostic()`](https://robustecologies.github.io/lucifer/reference/KS.Diagnostic.md)
performs a Kolmogorov-Smirnov test against a Student-t reference
distribution. These diagnostics are primarily useful for educational
purposes and for comparison with legacy analyses; Consort prefers the
modern suite.

### Specialized diagnostics

[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
provides NUTS-specific per-iteration diagnostics: energy statistics
(identifying poor energy conservation), tree depth distribution
(detecting max-treedepth saturation), divergent transition counts and
locations, and leapfrog step counts.
[`coupling_diagnostic()`](https://robustecologies.github.io/lucifer/reference/coupling_diagnostic.md)
computes L-lag total variation upper bounds via coupled Markov chains
(Biswas et al. 2019), providing a rigorous bound on the distance from
stationarity. PSIS diagnostics in
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
use the Pareto shape parameter \\k\\: \\k \< 0.5\\ indicates reliable
importance sampling, \\0.5 \leq k \< 0.7\\ indicates moderate
reliability, \\0.7 \leq k \< 1\\ indicates poor reliability requiring
moment-matching correction, and \\k \geq 1\\ indicates the importance
weights have infinite variance.

| Diagnostic | Measures | Multi-chain? | C++ backend? | Consort default? |
|:---|:---|:---|:--:|:--:|
| Rhat() | Split R-hat (rank-normalized) | yes\* | yes | yes |
| ESS() / ESS.tail() | Bulk/tail effective sample size | no | yes | yes |
| MCSE() / MCSE.quantile() | Monte Carlo standard error | no | yes | yes |
| Geweke.Diagnostic() | First vs. last chain segment means | no | no | no |
| Heidelberger.Diagnostic() | Stationarity + half-width CI | no | no | no |
| Raftery.Diagnostic() | Dependence factor, required length | no | no | no |
| Gelman.Diagnostic() | Classical multi-chain R-hat | yes | no | fallback |
| BMK.Diagnostic() | Batch means stationarity test | no | no | fallback |
| KS.Diagnostic() | K-S test vs. Student-t | no | no | no |
| check_nuts() | Energy, tree depth, divergences | no | no | yes |
| coupling_diagnostic() | L-lag TV upper bounds | yes | no | no |

Diagnostic method comparison. \*Rhat uses split chains internally;
single-chain input is split in half. Fallback = used when multi-chain
Rhat is unavailable. {.table style="width:100%;"}

  

## Cross-validation and model comparison

lucifer provides a complete model comparison toolkit for evaluating
predictive performance. The primary method is leave-one-out
cross-validation via Pareto smoothed importance sampling (LOO-PSIS),
which estimates out-of-sample predictive accuracy without refitting the
model. When PSIS diagnostics indicate that importance sampling is
unreliable (high Pareto \\k\\ values), K-fold cross-validation provides
a more robust but computationally expensive alternative. For time-series
models, leave-future-out cross-validation preserves temporal structure.
Model comparison is performed via
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
which computes pairwise ELPD differences with standard errors, and
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md),
which determines optimal prediction-averaging weights via constrained
optimization.

[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
requires a pointwise log-likelihood matrix (observations \\\times\\
posterior draws), which can be extracted from any fitting class via
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md).
The PSIS computation is implemented in C++ (`src/psis.cpp` with
`src/psis_internal.h`, ported from the loo R package v2.9) with OpenMP
parallelization across observations.
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md)
provides the widely applicable information criterion as a faster but
less diagnostic-rich alternative.

[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
refits the model \\K\\ times, each time leaving out one fold, and
computes the out-of-sample log-predictive density for the held-out
observations. When `CPUs > 1`, folds are evaluated in parallel via PSOCK
workers.
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
implements leave-future-out cross-validation for time series, refitting
the model at time points where the Pareto \\k\\ diagnostic exceeds a
threshold; refits are inherently sequential because each depends on the
posterior from the previous time window.

| Method | Function | Refitting? | Parallelizable? | Key diagnostic |
|:---|:---|:---|:---|:--:|
| LOO-PSIS | LOO() | no | n/a | Pareto k |
| WAIC | WAIC() | no | n/a | p_waic |
| K-fold | Kfold() | yes (K times) | yes (PSOCK) | ELPD |
| LFO | LFO() | yes (at trigger points) | no | Pareto k |
| loo_compare() | loo_compare() | no | n/a | ELPD difference SE |
| stacking_weights() | stacking_weights() | no | n/a | weight vector |

Cross-validation and model comparison methods. {.table}

  

## C++ backend architecture

Approximately 35,000 lines of C++ code provide the computational core of
lucifer. The source files in `src/` are organized by function: sampler
implementations, distribution evaluations, mathematical utilities,
diagnostic engines, and specialized inference backends. Header-only
libraries in `inst/include/lucifer/` provide reusable components for
analytical gradients, parameter transforms, thread-safe random number
generation, and automatic differentiation.

### Source file organization

The sampler files implement the MCMC kernels that
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
dispatches to. `sampler_gradient.cpp` is the largest single file
(approximately 6,600 lines), containing all gradient-based samplers
(NUTS, HMC, MALA, SGLD, autoMALA, BPS, ZigZag, Relativistic, and many
others) plus the NUTS tree-building machinery. `sampler_metropolis.cpp`,
`sampler_adaptive.cpp`, `sampler_slice.cpp`, `sampler_ensemble.cpp`,
`sampler_qmc.cpp`, and `sampler_other.cpp` contain the remaining MCMC
families. `sampler_geometric.cpp` implements Riemannian samplers (RMHMC,
LMC, MHMC) using the metric tensor engine from `metric_tensor.h` and
`metric_tensor.cpp`.

Distribution evaluations span three files:
`distributions_univariate.cpp` (40+ scalar distributions),
`distributions_multivariate.cpp` (MVN, MVT, Inverse-Wishart, Matrix
Normal), and `distributions_extended.cpp`. Mathematical utilities in
`math.cpp` and `matrices.cpp` provide safe log-sum-exp, log-determinant,
column/row variances, and symmetric matrix operations.

The specialized engines include `pathfinder.cpp` (L-BFGS trajectory
variational inference), `psis.cpp` (Pareto smoothed importance
sampling), `smc.cpp` (sequential Monte Carlo tempering),
`sde_engines.cpp` (SDE likelihood evaluation), `kalman.cpp` / `ukf.cpp`
/ `enkf.cpp` (Kalman filter family), `particle_filter.cpp` / `pgas.cpp`
/ `rbpf.cpp` (particle filter family), `node_bngm.cpp` (neural ODE
gradient matching), and `mdn.cpp` / `nre_classifier.cpp` /
`sbi_simulate.cpp` (SBI neural backends).

### The call dispatch

The central interface between C++ samplers and R model functions is
defined in `sampler_common.h` and `sampler_common.cpp`. Two functions
handle all model evaluation:

`call_model(Model, parm, Data)` invokes the R Model function via
`Rcpp::Function`, unpacks the five-element return list, validates the
output (checking for non-finite LP, correct parm length, try-error
results), and returns a `ModelOutput` struct. Before every R callback,
`PutRNGstate()` is called to synchronize the C-level RNG state with R’s
`.Random.seed`, preventing the state-loss bug that previously caused
deterministic (non-stochastic) sampler output.

`compute_gradient(parm, Data, grad_fn)` dispatches gradient computation
through three paths in priority order: (1) a user-supplied gradient
function passed via `Data$user.grad` or `Data$gradient`; (2) the
compiled model’s native gradient function pointer if a `compiled_model`
is being used; (3) finite-difference approximation via `partial_cpp()`,
which uses central differences with Richardson extrapolation. For
compiled models, both `call_model()` and `compute_gradient()` are
bypassed entirely in favor of direct function pointer invocation,
eliminating all R interpreter overhead.

### Header-only libraries

The `inst/include/lucifer/` directory contains four header-only C++
libraries designed for reuse across sampler files and in compiled
models:

`dist_gradients.h` provides analytical gradients of log-density
functions for 24 distributions (Normal, Student-t, Cauchy, Laplace,
Gamma, Beta, LogNormal, and others), including chain rule support for
composed transformations. These are used by the C++ code generator to
avoid finite-difference approximations in compiled models.

`transforms.h` implements five constraint transform types (log for
positive parameters, logit for unit interval, generalized logit for
bounded intervals, stick-breaking for simplexes, and Cholesky for
correlation matrices), each with the log-Jacobian-determinant adjustment
required for correct posterior sampling on constrained supports.

`thread_rng.h` provides the xoshiro256+ pseudorandom number generator
(Blackman & Vigna, 2018) with per-thread state for safe use in OpenMP
parallel regions. R’s built-in RNG (`R::rnorm()`, `R::runif()`) cannot
be called from parallel threads because it is guarded by a process-wide
lock; xoshiro256+ provides a fast, high-quality alternative with
deterministic seeding from R’s master-thread RNG.

`dual.h` implements forward-mode automatic differentiation via dual
numbers (approximately 170 lines, supporting 20+ elementary functions).
This provides an alternative to both analytical gradients and finite
differences for computing gradients of compiled models.

C++ call flow showing the dispatch path from the R-level dispatcher
through the C++ sampler loop to model evaluation.

| File | Lines | Purpose |
|:---|---:|:---|
| sampler_gradient.cpp | ~6600 | All gradient-based, PDMP, tempering, constraint, flow samplers |
| sampler_common.cpp | ~800 | call_model, compute_gradient, interrupt, progress |
| sampler_metropolis.cpp | ~1200 | RWM, IM, MTM, DRM, RDMH |
| sampler_adaptive.cpp | ~1400 | AM, AMM, AMWG, ADMG, DRAM, RAM, INCA |
| sampler_slice.cpp | ~900 | Slice, AFSS, ESS, OHSS, RSS, UESS |
| sampler_ensemble.cpp | ~1100 | AIES, DEMC, MCMCMC, t-walk |
| sampler_qmc.cpp | ~1000 | QMC family (8 methods) |
| sampler_other.cpp | ~600 | pCN, RJ, PMCMC, misc |
| sampler_geometric.cpp | ~400 | RMHMC, LMC, MHMC + metric tensor |
| distributions_univariate.cpp | ~2500 | 40+ scalar distribution densities |
| distributions_multivariate.cpp | ~1500 | MVN, MVT, IW, MatrixNormal |
| distributions_extended.cpp | ~800 | Extended distribution set |
| math.cpp | ~500 | LogAdd, LogDet, matrix utilities |
| matrices.cpp | ~400 | Column/row variances, symmetric ops |
| diagnostics.cpp | ~600 | Rhat, ESS, MCSE (modern, OpenMP) |
| pathfinder.cpp | ~500 | L-BFGS trajectory VI |
| psis.cpp | ~800 | PSIS-LOO (ported from loo v2.9) |
| smc.cpp | ~700 | SMC adaptive tempering + rejuvenation |
| sde_engines.cpp | ~1000 | EKF, Kalman, SDE likelihood |
| kalman.cpp | ~600 | Kalman filter/smoother/FFBS |
| ukf.cpp | ~300 | Unscented Kalman filter |
| enkf.cpp | ~350 | Ensemble Kalman filter |
| particle_filter.cpp | ~500 | Conditional PF, auxiliary PF, bridge PF |
| pgas.cpp | ~400 | Particle Gibbs with ancestor sampling |
| ssm_engines.cpp | ~1400 | KSC (SV), SMC-squared |
| node_bngm.cpp | ~620 | Neural ODE Bayesian gradient matching |
| mdn.cpp | ~800 | Mixture density network (SBI) |
| nre_classifier.cpp | ~500 | Neural ratio estimator (SBI) |
| dist_gradients.h\* | ~545 | Analytical gradients, 24 distributions |
| transforms.h\* | ~213 | Log, logit, bounded, simplex, Cholesky transforms |
| thread_rng.h\* | ~150 | xoshiro256+ per-thread PRNG |
| dual.h\* | ~170 | Dual-number forward-mode AD |

C++ source file summary. \*Header-only libraries in
inst/include/lucifer/. Line counts are approximate. {.table}

  

## Parallelism

### Architecture overview

lucifer exploits parallelism at three levels, each targeting a different
granularity of work and controlled by a centralized thread budget
manager.

**C++ OpenMP (intra-process).** Distribution density evaluations, MCMC
diagnostics (Rhat, ESS), PSIS importance weight smoothing, neural
network forward passes, NODE ensemble optimization, and
particle/ensemble Kalman filter propagation are parallelized via OpenMP
pragmas in C++. These operate on shared memory within a single R process
and require no user intervention; the number of threads is controlled by
the thread budget. All parallel regions are guarded with
`#ifdef _OPENMP` and degrade gracefully to serial execution when
compiled without OpenMP.

**Multi-chain MCMC (inter-process).** When
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
is called with `Chains > 1` and `CPUs > 1`, each chain runs in a
separate R subprocess via
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html) (or falls
back to a PSOCK cluster if callr is unavailable). Subprocesses have
isolated memory spaces, eliminating all thread-safety concerns. The
thread budget controls how many OpenMP threads each subprocess gets,
allowing C++ parallelism to operate within each chain.

**Task-level parallelism (cross-platform PSOCK).** Embarrassingly
parallel workloads such as K-fold cross-validation (`Kfold`), SBI
forward simulation (`SBI`), VB Pathfinder multi-path optimization
(`VariationalBayes` with `Method = "Pathfinder"`), and ABC rejection
sampling (`ABC` with `CPUs > 1`) use a cross-platform PSOCK dispatch
helper that works on Linux, macOS, and Windows.

### The thread budget API

The
[`lucifer_threads()`](https://robustecologies.github.io/lucifer/reference/lucifer_threads.md)
function is the single entry point for controlling how computational
resources are distributed.

``` r

library(lucifer)

# Query current settings
lucifer_threads()

# Set total CPUs and let the auto strategy distribute them
lucifer_threads(total = 8)

# Force all threads to chain parallelism, no OpenMP (safest)
lucifer_threads(total = 8, strategy = "chains_only")

# Force all threads to OpenMP, chains run sequentially
lucifer_threads(total = 8, strategy = "omp_only")

# Manual override
lucifer_threads(omp = 2, chains = 4)
```

#### Strategies

The `strategy` parameter controls how CPUs are split between inter-chain
workers and intra-chain OpenMP threads:

- **`"auto"` (default):** distributes threads evenly. On 8 CPUs with 4
  chains, each chain gets 2 OpenMP threads. This gives the best overall
  throughput for models that benefit from vectorized distribution
  evaluations.

- **`"chains_only"`:** all CPUs to chain workers, each chain runs
  single-threaded. This is the safest option and avoids any potential
  for nested parallelism conflicts with threaded BLAS libraries
  (OpenBLAS, MKL). Use this if you observe deadlocks or erratic
  behavior.

- **`"omp_only"`:** all CPUs to OpenMP within a single process. Chains
  run sequentially. Best for a single chain with large data and many
  vectorized distribution evaluations.

- **`"balanced"`:** synonym for `"auto"`.

#### Diagnostics

``` r

# Full diagnostic report
lucifer_parallel_info()
```

This prints detected cores, OpenMP availability, current budget, and
platform-specific guidance (e.g., how to install OpenMP on macOS).

### Multi-chain MCMC

The `Chains` and `CPUs` parameters in
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
control multi-chain execution:

``` r

fit <- lucifer(
    Model = MyModel, Data = MyData,
    Initial.Values = iv,
    Algorithm = "NUTS",
    Iterations = 5000,
    Chains = 4,
    CPUs = 4
)
```

When `Chains > 1` and `CPUs > 1`, the thread budget determines how many
OpenMP threads each chain subprocess receives. With the default `"auto"`
strategy, 4 CPUs and 4 chains gives each chain 1 OpenMP thread; 8 CPUs
and 4 chains gives each chain 2 OpenMP threads.

The dispatch mechanism is:

1.  If `callr` is available (recommended): each chain runs as a
    [`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html)
    subprocess with per-iteration progress bars.
2.  If `callr` is unavailable: a PSOCK cluster is used (no per-iteration
    progress output).
3.  FORK clusters are never used because forking is unsafe with OpenMP
    (forked children inherit the parent’s OpenMP runtime but only one
    thread survives, causing deadlocks).

BLAS threading is always set to single-threaded in child processes to
avoid contention when multiple processes share CPU cores. This is
separate from the OpenMP thread count, which is set per the budget.

#### Wrappers

[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
and
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
accept `Chains` and `CPUs` and pass them to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
The same thread budget applies.

### Cross-validation parallelism

#### K-fold cross-validation

``` r

loo_k <- Kfold(fit, K = 10, CPUs = 4)
```

The `CPUs` parameter in
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
controls how many folds are evaluated in parallel. The dispatch uses
PSOCK clusters, so it works on all platforms including Windows.

#### Leave-future-out (LFO)

[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
refits the model at time points where the Pareto k diagnostic exceeds a
threshold. These refits are inherently sequential because each depends
on the posterior from the previous one, so LFO does not support parallel
execution.

### SBI parallel simulation

``` r

result <- SBI(
    sim_fn = my_simulator,
    prior = my_prior,
    obs = observed_data,
    Method = "NPE",
    n_sims = 10000,
    n_cores = 4
)
```

The `n_cores` parameter in
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
controls how many parallel workers run the simulator function. Each
worker draws from the prior, simulates data, and returns the result. The
dispatch uses PSOCK clusters and works on all platforms.

### VB Pathfinder multi-path

Pathfinder runs multiple L-BFGS optimization paths from jittered initial
values, each producing a normal approximation. The paths are independent
and run in parallel when multiple CPUs are available:

``` r

lucifer_threads(total = 4)
vb_fit <- VariationalBayes(
    Model = MyModel,
    Data = MyData,
    parm = iv,
    Method = "Pathfinder",
    Covar = 4,        # number of paths
    Iterations = 1000
)
```

When `Covar > 1` (number of paths) and the thread budget has
`total_cpus > 1`, paths run in parallel via PSOCK workers. The best path
is selected by ELBO, and multi-path PSIS resampling combines draws from
all paths.

### ABC rejection parallelism

``` r

abc_fit <- ABC(
    Model = my_simulator,
    Data = my_data,
    Prior = my_prior,
    Summary.Stats = my_stats,
    Observed.Stats = obs_stats,
    Method = "rejection",
    N = 5000,
    CPUs = 4
)
```

When `CPUs > 1` and `Method = "rejection"`, the target sample count is
split across workers, each running an independent rejection loop.
MCMC-ABC, SMC-ABC, and SA-ABC remain sequential because their state is
not decomposable.

### Platform-specific notes

#### Linux

Full support. OpenMP is available via GCC’s `-fopenmp`. Both callr-based
and PSOCK-based dispatch work. `mclapply` has been replaced with PSOCK
clusters throughout, so behavior is identical to other platforms.

#### Windows

OpenMP is available via the Rtools mingw-w64 toolchain. PSOCK clusters
work. There is no FORK support on Windows, but lucifer never uses FORK
(it is unsafe with OpenMP).

#### macOS

Apple Clang does not bundle OpenMP. Without it, all C++ code runs
single-threaded (the package still compiles and works correctly). To
enable OpenMP:

1.  Install libomp: `brew install libomp`
2.  Add to `~/.R/Makevars`:

&nbsp;

    CPPFLAGS += -Xclang -fopenmp
    LDFLAGS += -lomp
    SHLIB_OPENMP_CFLAGS = -fopenmp
    SHLIB_OPENMP_CXXFLAGS = -fopenmp

3.  Reinstall: `install.packages("lucifer", type = "source")`

Run
[`lucifer_parallel_info()`](https://robustecologies.github.io/lucifer/reference/lucifer_parallel_info.md)
to verify OpenMP is available after reinstallation.

#### CRAN compliance

lucifer respects the `_R_CHECK_LIMIT_CORES_` environment variable. When
set to `"true"` or `"warn"` (as on CRAN build machines), the total CPU
budget is capped at 2. All parallel tests are wrapped in
`skip_on_cran()`.

### Performance guidelines

The optimal strategy depends on the model and data size:

**Small models (few parameters, fast likelihood).** Chain parallelism
dominates. Use `strategy = "chains_only"` or the default `"auto"`. The
overhead of OpenMP thread creation is not worth it when each likelihood
evaluation takes microseconds.

**Large-n models with vectorized distributions.** When the model
evaluates vectorized C++ distribution functions on thousands of
observations, OpenMP within each chain provides meaningful speedup. Use
`strategy = "auto"` or `"balanced"` so each chain gets multiple OpenMP
threads.

**Single-chain exploratory runs.** Use `strategy = "omp_only"` to give
all threads to OpenMP within the single chain.

**Compiled models.** When using
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
on a `model_spec`, the model evaluation runs entirely in C++ with no R
callbacks. This enables additional parallelism in SMC particle
rejuvenation (conditional OpenMP in the inner particle loop).

### Thread-safe RNG

OpenMP parallel regions cannot call R’s random number generator
(`R::rnorm()`, `R::runif()`) because R’s RNG is guarded by a
process-wide lock. lucifer provides the xoshiro256+ PRNG (Blackman &
Vigna, 2018), a fast, high-quality generator with per-thread state. Each
thread gets its own RNG seeded deterministically from R’s RNG in the
master thread before entering the parallel region. This is used in
particle filter propagation, EnKF ensemble propagation, NODE ensemble
initialization, and SMC rejuvenation with compiled models.

### Internals for developers

#### `.lucifer_parallel_apply()`

Package contributors should use this internal helper for embarrassingly
parallel tasks:

``` r

results <- .lucifer_parallel_apply(
    X = seq_len(100),
    FUN = function(i) expensive_computation(i),
    n_workers = 4,
    export_vars = c("expensive_computation"),
    export_env = environment()
)
```

It handles CRAN limits, always uses PSOCK (never FORK), loads lucifer in
workers, and falls back to `lapply` when `n_workers <= 1`.

#### `.make_cluster()`

For code that needs a persistent cluster across iterations (e.g., GG,
AGG, MCMCMC samplers), use `.make_cluster()` to create the cluster with
proper setup, then call `parLapply(cl, ...)` directly:

``` r

cl <- .make_cluster(CPUs, export_vars = c("Model", "Data"),
                     export_env = environment())
on.exit(stopCluster(cl))
# ... use cl in loop ...
```

#### R API safety in OpenMP regions

Never call R functions (including `Rcpp::Function`, `R::rnorm`,
`Rcpp::wrap`) inside `#pragma omp parallel for` regions. R’s internal
memory manager and GC are not re-entrant. Use the xoshiro256+ RNG from
`inst/include/lucifer/thread_rng.h` for random number generation in
parallel regions.

  

## Graceful interruption

Long-running inference is inherently unpredictable: a model that was
expected to converge in minutes may require hours, or a quick
exploratory run may reveal a problem that makes the current fit
unnecessary. lucifer implements a cooperative interruption system that
lets the user press Esc (or Ctrl+C) at any point during inference and
receive a usable partial result rather than losing all progress. The
system operates across three execution contexts, single-process C++
samplers, R-level iterative methods, and multi-process parallel chains,
each with its own signaling mechanism but sharing a common
partial-result caching strategy.

### Single-chain C++ MCMC

Every C++ sampler loop calls `maybe_check_interrupt()` once per
iteration. This function checks R’s interrupt-pending flag via
`Rcpp::checkUserInterrupt()`, but only every 200 milliseconds (tracked
via `std::chrono::steady_clock`) to avoid excessive overhead; the cost
between actual checks is approximately 20-50 nanoseconds per iteration.
When the user presses Esc, R sets its internal interrupt flag; the next
`Rcpp::checkUserInterrupt()` call throws an `InterruptedException`,
which is caught inside `call_model()`. The catch block sets the global
volatile flag `g_mcmc_interrupted = true` and returns an invalid
`ModelOutput`. The sampler loop checks this flag at the top of its next
iteration and breaks cleanly.

Before the interrupt is detected, the sampler has been periodically
caching its state via `.cache_partial()` from `R/interrupt_utils.R`.
This function stores the current posterior matrix, LP trace, Monitor
values, and acceptance statistics in a package-private environment
(`.interrupt_env$partial`) at every thinning point. When the C++ sampler
returns (with incomplete output or NULL), the
`tryCatch(interrupt = ...)` wrapper in
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
retrieves the cached state via `.get_cached_partial()`, reconstructs a
valid `demonoid` S3 object from it, and marks it with
`$Interrupted = TRUE`. The cache is cleared immediately after retrieval
to prevent stale data from leaking into future runs.

### R-level iterative methods

R-level samplers (Gibbs, GG, AGG, LSS, QSS, GPSS, NeuTra, flowMC,
MCMCMC) and non-MCMC engines (VariationalBayes, ABC, PMC,
LaplaceApproximation, IterativeQuadrature) follow the same pattern:
every 50 iterations, they call `.cache_partial(current_state)` followed
by `check_interrupt_r()`, an exported C++ function that wraps
`Rcpp::checkUserInterrupt()`. If the interrupt is detected, R throws a
condition of class `interrupt`, which is caught by the outer `tryCatch`.
The handler retrieves the cached state and returns it as a partial
result.

For variational inference, the cached state includes the current
approximate posterior mean, covariance estimate, ELBO trajectory, and
step size, so the returned `vb` object is a genuine intermediate
approximation rather than a stub. For PMC and iterative quadrature, the
cached state includes the current particle set or quadrature grid,
respectively. For Laplace approximation, the optimizer’s current best
point and Hessian estimate are cached.

### Multi-chain parallel execution

Multi-chain interruption requires cross-process coordination. lucifer
supports two parallel backends, each with its own signaling mechanism:

**callr backend.** Each chain runs in a separate R subprocess spawned by
[`callr::r_bg()`](https://callr.r-lib.org/reference/r_bg.html). When the
user presses Esc during the progress-polling loop, the master process
catches the interrupt and executes a three-stage shutdown: first, it
sends `SIGINT` to each alive child via `handle$interrupt()`; second, it
waits up to 3 seconds for children to finish gracefully (checking
`handle$is_alive()` every 200ms); third, it forcibly kills any survivors
via `handle$kill()`. Each child process that receives the signal runs
the same single-chain interrupt path described above, caching partial
results and returning them. The master collects whatever results are
available, marks all chains with `$Interrupted = TRUE`, and assembles
the combined object, computing Rhat and ESS only from chains with
sufficient samples.

**PSOCK backend (fallback).** PSOCK cluster workers cannot receive
`SIGINT` directly because they are spawned through intermediate wrapper
processes. Instead, lucifer uses a cooperative file-based protocol:
before launching workers, the master creates a stop-signal file path (a
tempfile that does not yet exist) and registers it in each worker’s C++
global via `set_stop_signal_cpp(path)`. When the user presses Esc, the
master writes “1” to the file via `.signal_stop(path)`. Worker processes
detect the file within 200ms via `maybe_check_interrupt()`, which checks
`file.exists(g_stop_signal_file)` alongside the standard R interrupt
flag. Once detected, `g_mcmc_interrupted` is set to true and the sampler
loop breaks. The stop-signal file is cleaned up via
`.clear_stop_signal(path)` in an
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) handler.

Interrupt handling flow across single-chain, callr multi-chain, and
PSOCK multi-chain execution contexts.

### Output recovery guarantees

Every engine in lucifer returns a usable partial result on interrupt.
The `$Interrupted` field is set to `TRUE` on the output object, and
[`print()`](https://rdrr.io/r/base/print.html) methods display a yellow
warning banner: “\\\mathtt{\unicode{26A0}}\\ Results from interrupted
run (partial).” For multi-chain fits, chains killed before returning any
samples are replaced with minimal stubs (`$Iterations = 0`); these stubs
are excluded from diagnostic computations like Rhat and ESS, so the
surviving chains’ results are not corrupted by the inclusion of empty
data.

The caching strategy is designed for minimal overhead:
`.cache_partial()` performs a single environment assignment (~1
microsecond), and `maybe_check_interrupt()` performs a clock comparison
(~20-50 nanoseconds per call) with the full R interrupt check happening
only every 200 milliseconds. For a sampler running at 1000 iterations
per second, the combined overhead of caching and interrupt polling is
well below 0.01% of total runtime.

  

## Interoperability

lucifer provides bidirectional bridges to Stan, brms, JAGS, coda, and
the posterior R package, enabling users to import external fits into the
lucifer ecosystem for post-processing and to export lucifer fits into
formats expected by other tools.

### Importing external fits

The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
generic converts foreign fit objects into lucifer’s `demonoid` class.
Supported input types include `stanfit` (rstan), `brmsfit` (brms),
`CmdStanMCMC` (cmdstanr), `mcmc.list` (coda/rjags), `runjags` objects,
and `draws_array` / `draws_matrix` / `draws_df` from the posterior
package. The converted object gains access to the full lucifer
post-processing ecosystem: Rhat, ESS, MCSE, plot methods, Consort
diagnostics, LOO, and WAIC. Generated quantities like
`log_lik[1]..log_lik[N]` are automatically filtered from the posterior
matrix during Stan imports.

### Backend fitting

[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
and
[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)
enable fitting models using Stan or JAGS as the sampling backend while
returning lucifer `demonoid` objects. This is useful for users who have
existing Stan or JAGS model code and want to use lucifer’s diagnostic,
comparison, and orchestration tools without rewriting their models. The
returned objects support
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
and all other demonoid methods.

### Exporting lucifer fits

[`to_mcmc_list()`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md)
converts a demonoid object to a coda `mcmc.list`, enabling compatibility
with packages that expect coda objects.
[`to_draws_array()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md),
[`to_draws_matrix()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md),
and
[`to_draws_df()`](https://robustecologies.github.io/lucifer/reference/to_draws_array.md)
convert to posterior package formats, enabling use with the bayesplot,
loo, and priorsense ecosystems. Experimental translation functions
[`stan_to_spec()`](https://robustecologies.github.io/lucifer/reference/stan_to_spec.md)
and
[`jags_to_spec()`](https://robustecologies.github.io/lucifer/reference/jags_to_spec.md)
attempt to convert Stan or JAGS model code into lucifer’s `model_spec`
notation.

### Namespace collision handling

Several functions exported by brms and loo (`log_lik`, `LOO`,
`loo_compare`, `Rhat`) share names with lucifer exports. The `.onLoad`
and `.onAttach` hooks in `R/lucifer-package.R` register
`log_lik.demonoid` on brms/loo generics via
[`registerS3method()`](https://rdrr.io/r/base/ns-internal.html), and
`setHook(packageEvent())` handlers ensure correct method registration
when brms or loo loads after lucifer. In vignettes and scripts where
both packages are loaded, using fully qualified names
([`lucifer::log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`lucifer::LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md))
is the safest approach.

| External format | Import function | Export function | log_lik? | Full diagnostics? |
|:---|:---|:---|:--:|:--:|
| stanfit (rstan) | as.demonoid() | – | yes | full |
| brmsfit (brms) | as.demonoid() | – | yes | full |
| CmdStanMCMC (cmdstanr) | as.demonoid() | – | yes | full |
| mcmc.list (coda) | as.demonoid() | to_mcmc_list() | manual | full |
| runjags | as.demonoid() | – | manual | full |
| draws\_\* (posterior) | as.demonoid() | to_draws\_\*() | manual | full |

Interoperability matrix. ‘manual’ = user must supply pointwise
log-likelihood separately. ‘full’ = Rhat, ESS, MCSE, plot, Consort,
LOO/WAIC all available after conversion. {.table}

  

## Domain-specific inference engines

Three specialized engines extend lucifer’s general-purpose inference
into specific application domains. Each accepts domain-specific model
structures and produces an S3 fit object with tailored print, summary,
and plot methods. These engines internally dispatch to lucifer’s MCMC
samplers for posterior computation while adding domain-specific
pre-processing (model construction, likelihood specification) and
post-processing (eigenanalysis, forecasting, trajectory visualization).

### Neural ODEs: `NODE()`

[`NODE()`](https://robustecologies.github.io/lucifer/reference/NODE.md)
implements Bayesian gradient matching for system identification using
the method of Bonnaffé (2023). Given multivariate time-series data, it
fits a neural ODE model with dual-pathway single-layer perceptrons and
sinusoidal interpolation. The C++ backend (`src/node_bngm.cpp`,
approximately 620 lines) implements the sinusoidal network architecture,
BFGS optimization, and ensemble learning. Ensemble members are optimized
in parallel when multiple CPUs are available. The output `node_fit`
class provides 8 plot types including network diagrams and RK4
forecasting trajectories.

### Stochastic differential equations: `SDE()`

[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
provides Bayesian inference for stochastic differential equations with
11 model families registered in the SDE family system. Likelihood
evaluation engines include exact transition densities (when available),
Euler-Maruyama and Milstein discretization, the extended Kalman filter
(EKF), and particle filters (bootstrap, auxiliary, and bridge). The C++
backends span `src/sde_engines.cpp`, `src/particle_filter.cpp`, and the
Kalman filter family. The `sde_fit` class supports predict methods for
trajectory forecasting and is compatible with Consort’s SDE-specific
diagnostic engine.

### State-space models: `SSM()`

[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
constructs and fits state-space models with five model types: local
level, local linear trend, seasonal, basic structural model (BSM), and
VAR(p). Filtering and smoothing are performed via the Kalman filter
(`src/kalman.cpp`), unscented Kalman filter (`src/ukf.cpp`), ensemble
Kalman filter (`src/enkf.cpp`), Rao-Blackwellized particle filter
(`src/rbpf.cpp`), or particle Gibbs with ancestor sampling
(`src/pgas.cpp`), depending on the model’s linearity and Gaussianity.
Six observation families are supported: Gaussian, Poisson, negative
binomial, Student-t, zero-inflated Poisson, and binomial.
Markov-switching regime detection uses the Hamilton filter with
forward-filtering backward-sampling (MS-FFBS).

| Engine | Model families |    C++ backend    | Output class | Vignette |
|:-------|:---------------|:-----------------:|:-------------|:---------|
| NODE() | 1 (BNGM)       |   node_bngm.cpp   | node_fit     | node.Rmd |
| SDE()  | 11             | sde_engines.cpp + | sde_fit      | sde.Rmd  |
| SSM()  | 5              |   kalman.cpp +    | ssm_fit      | ssm.Rmd  |

Domain-specific inference engines. ‘+’ indicates multiple C++ source
files. {.table}

  

## Extending lucifer

lucifer’s registry-driven architecture is designed so that adding a new
algorithm automatically integrates it with the entire orchestration
ecosystem. The critical insight is that Prescribe, Consort, Arena, and
Crucible all consume algorithm metadata from the centralized registry
rather than maintaining their own lists; once an algorithm is
registered, these components discover it without modification.

### Adding a new MCMC algorithm

The complete workflow involves seven files across three layers: the
sampler implementation, the registry and dispatcher, and the
orchestration integration. The following diagram shows the dependencies
between these steps.

Workflow for adding a new MCMC algorithm to lucifer. Files to modify are
shown with their roles; green nodes indicate automatic integration via
the registry.

**Step 1: sampler implementation.** Create `R/sampler_myalgo.R`
containing a function
`.mcmcmyalgo(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, ...)`.
The function must return a list with elements `Posterior` (thinned
samples matrix), `Summary` (parameter summary statistics), `LML` (log
marginal likelihood estimate), `Monitor` (monitored quantities matrix),
`Acceptance` (acceptance rate), and `Minutes` (elapsed wall-clock time).
Examine any existing sampler file (e.g., `R/sampler_nuts.R`) for the
exact return structure. For interrupt support, call `.cache_partial()`
at each thinning point and wrap the main loop in
`tryCatch(interrupt = ...)`.

**Step 1b (optional): C++ core.** For performance-critical algorithms,
implement the iteration loop in `src/sampler_*.cpp`. Use `call_model()`
from `sampler_common.h` to evaluate the Model function, and
`compute_gradient()` for gradient computation. Call
`maybe_check_interrupt()` at the top of each iteration for interrupt
support. Export the C++ function via `// [[Rcpp::export]]` and call it
from the R wrapper.

**Step 2: registry entry.** Add an entry to `.mcmc_registry` in
`R/algo_registry.R` via `.make_algo_entry()`. The 16 metadata fields
drive all downstream behavior: `requires_gradient` determines whether
Prescribe scores the algorithm for gradient-available models;
`acceptance_range` defines the window that Consort uses to diagnose
acceptance rate problems; `quality_tier` influences Crucible’s method
selection priority; `dim_range` allows Prescribe to disqualify the
algorithm for problems outside its effective dimensionality range. Set
these fields carefully, as they directly control the orchestration
layer’s behavior.

**Step 3: dispatcher case.** Add a case to the `switch(Algorithm)` block
in `R/lucifer.R` that calls `.mcmcmyalgo()` with the appropriate
arguments. Follow the pattern of existing cases.

**Step 4: default Specs.** Add a default Specs entry in `.algo_defaults`
within `R/specs_validation.R`. This ensures that
`lucifer(Algorithm = "MyAlgo")` works without the user specifying Specs
explicitly, and that Crucible can instantiate the algorithm with
sensible defaults.

**Step 5: escalation path.** Add the algorithm to the escalation map in
`R/consort_suggest.R`. This map defines what Consort suggests when the
algorithm fails convergence: for example, if your algorithm is a
gradient-based sampler, the escalation path might be `MyAlgo -> NUTS`.
If Consort does not find the algorithm in its map, it falls back to
generic suggestions (increase iterations), but explicit escalation paths
produce better refinement in Crucible.

**Step 6: documentation.** Add the algorithm to
`vignettes/mcmc-algorithms.Rmd` with a brief description, the `Specs`
parameters, and a reference citation. Add roxygen documentation to the
sampler file if it contains any internal helper functions worth
documenting.

**Step 7: build.** Run
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
(if C++ was added), then
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
then
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html).

After these steps, the algorithm works transparently with
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
(scored based on registry metadata),
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
(diagnosed using `acceptance_range` and standard convergence criteria),
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
(extracted via the existing `demonoid` adapter, since all MCMC
algorithms produce demonoid objects), and
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
(selected, fitted, and refined via the orchestration pipeline).

### Adding a new inference engine

Adding an entirely new engine family requires integration at more points
than adding an MCMC algorithm, because the output class is new and the
orchestration tools need class-specific adapters to process it.

Workflow for adding a new inference engine to lucifer. Unlike MCMC
algorithms, a new engine requires explicit adapters for Arena, Consort,
and optionally PPC.

**Step 1: engine function.** Create the main R function with the
standard `Model(parm, Data)` interface. It must accept the same Model
contract (five-element return list) used by all other engines. The
engine implements its inference strategy and packages the results into a
new S3 class.

**Step 2: S3 output class.** Define a new class (e.g.,
`"my_engine_fit"`) for the engine’s output. The object should contain at
minimum: posterior samples or approximation parameters, the algorithm
name, wall-clock time, convergence status, and the original call.

**Step 3: method triad.** Implement
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
[`print()`](https://rdrr.io/r/base/print.html) should produce a concise
one-screen synopsis; [`summary()`](https://rdrr.io/r/base/summary.html)
should provide extended diagnostics;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) should default
to the most informative visualization and support a `type` argument for
alternatives.

**Step 4: registry sub-list.** Create a new registry in
`R/algo_registry.R` (e.g., `.my_engine_registry`) using
`.make_algo_entry()` for each algorithm variant, then append it to
`.algo_registry`. This makes the engine visible to
[`algo_info()`](https://robustecologies.github.io/lucifer/reference/algo_info.md)
and to Prescribe’s scoring system.

**Step 5: Arena adapter.** Add `.arena_extract_my_engine_fit()` in
`R/arena_extract.R`. This function must normalize the engine’s output
into Arena’s standardized extraction structure: posterior mean,
posterior SD, posterior samples, minimum ESS, ESS per second, wall-clock
time, convergence status, and log marginal likelihood. Without this
adapter, Arena cannot compare the engine’s results against other
methods.

**Step 6: Consort engine.** Add `.consort_my_engine()` in
`R/consort_families.R`. This function evaluates the fit against
convergence criteria appropriate to the method and returns a structured
list with `conditions` (a data frame of condition/status/detail rows),
`diagnostics`, and optionally a `suggestion` for refinement.

**Step 7: Consort dispatch.** Add class detection to `R/Consort.R` so
that `Consort(my_fit)` routes to the new diagnostic engine. This is
typically a single [`inherits()`](https://rdrr.io/r/base/class.html)
check added before the existing class switch.

**Steps 8-10 (optional).** Implement
[`predict()`](https://rdrr.io/r/stats/predict.html) for posterior
predictive sampling, PPC support via
[`as.ppc()`](https://robustecologies.github.io/lucifer/reference/as.ppc.md)
and `plot.my_engine_fit.ppc()`, and add the engine family to Crucible’s
`families` parameter list so that Crucible can include it in automated
pipeline runs.

  

## Quick reference

| Function | Purpose | Key arguments | Output |
|:---|:---|:---|:---|
| lucifer() | MCMC sampling (82 algorithms) | Algorithm, Specs, Chains, CPUs | demonoid |
| VariationalBayes() | Variational inference (8 methods) | Method, Iterations, Covar | vb |
| LaplaceApproximation() | Posterior mode finding (18 optimizers) | Method, CovEst, Iterations | laplace |
| IterativeQuadrature() | Gaussian quadrature (3 rules) | Method, Iterations | iterquad |
| PMC() | Population Monte Carlo | N, Iterations, Thinning | pmc |
| SMC() | Sequential Monte Carlo | N, Iterations | smc |
| ABC() | Approximate Bayesian computation | Method, N, CPUs | abc |
| SBI() | Simulation-based inference (neural) | Method, n_sims, n_cores | sbi |
| model_spec() | Declarative model specification DSL | formula, backend | model_spec |
| compile() | Compile model_spec to native C++ | x (model_spec object) | compiled_model |
| Prescribe() | Pre-fit algorithm profiling and ranking | Model, Data, Initial.Values | prescription |
| Consort() | Post-fit convergence diagnostics | x (any fit object) | consort |
| Arena() | Cross-method benchmarking | x (list of fits), reference | arena |
| Crucible() | Automated full inference pipeline | Model, Data, n_methods, max_rounds | crucible |
| LOO() | Leave-one-out cross-validation (PSIS) | x (fit), log_lik | loo |
| WAIC() | Widely applicable information criterion | x (fit), log_lik | waic |
| Kfold() | K-fold cross-validation | x (fit), K, CPUs | kfold |
| LFO() | Leave-future-out cross-validation | x (fit), L, threshold | lfo |
| loo_compare() | Compare models by ELPD | … (loo objects) | data.frame |
| stacking_weights() | Bayesian model averaging weights | … (loo objects) | numeric |
| Rhat() | Split R-hat (rank-normalized) | x (matrix or demonoid) | numeric |
| ESS() | Effective sample size (bulk/tail) | x (matrix or demonoid) | numeric |
| MCSE() | Monte Carlo standard error | x (matrix or demonoid) | numeric |
| as.demonoid() | Import external fit to demonoid | x (stanfit, brmsfit, etc.) | demonoid |
| to_mcmc_list() | Export demonoid to coda mcmc.list | x (demonoid) | mcmc.list |
| lucifer_threads() | Configure thread budget | total, strategy, omp, chains | invisible |
| lucifer_parallel_info() | Parallel system diagnostics | (none) | invisible |

Primary API function reference. {.table}

### MCMC algorithm quick reference

| Algorithm    | Family        | Gradient? | Discrete? | Tier |
|:-------------|:--------------|:---------:|:---------:|:----:|
| AHMC         | gradient      |    yes    |    no     |  2   |
| HMC          | gradient      |    yes    |    no     |  2   |
| HMCDA        | gradient      |    yes    |    no     |  1   |
| NUTS         | gradient      |    yes    |    no     |  1   |
| MALA         | gradient      |    yes    |    no     |  2   |
| THMC         | gradient      |    yes    |    no     |  2   |
| SGLD         | gradient      |    yes    |    no     |  3   |
| SGHMC        | gradient      |    yes    |    no     |  3   |
| MCHMC        | gradient      |    yes    |    no     |  2   |
| Barker       | gradient      |    yes    |    no     |  1   |
| autoMALA     | gradient      |    yes    |    no     |  1   |
| MALT         | gradient      |    yes    |    no     |  1   |
| AAPS         | gradient      |    yes    |    no     |  1   |
| Relativistic | gradient      |    yes    |    no     |  2   |
| NROLangevin  | gradient      |    yes    |    no     |  2   |
| GIST         | gradient      |    yes    |    no     |  1   |
| AM           | adaptive_rw   |    no     |    no     |  2   |
| AMM          | adaptive_rw   |    no     |    no     |  2   |
| RAM          | adaptive_rw   |    no     |    no     |  1   |
| DRAM         | adaptive_rw   |    no     |    no     |  1   |
| RWM          | adaptive_rw   |    no     |    no     |  3   |
| DRM          | adaptive_rw   |    no     |    no     |  2   |
| pCN          | adaptive_rw   |    no     |    no     |  2   |
| MTM          | adaptive_rw   |    no     |    no     |  2   |
| RDMH         | adaptive_rw   |    no     |    no     |  2   |
| MCMCMC       | adaptive_rw   |    no     |    no     |  2   |
| IM           | adaptive_rw   |    no     |    no     |  3   |
| INCA         | adaptive_rw   |    no     |    no     |  2   |
| AMWG         | componentwise |    no     |    yes    |  1   |
| MWG          | componentwise |    no     |    yes    |  2   |
| Gibbs        | componentwise |    no     |    yes    |  1   |
| ADMG         | componentwise |    no     |    yes    |  2   |
| GG           | componentwise |    no     |    yes    |  3   |
| AGG          | componentwise |    no     |    yes    |  3   |
| SAMWG        | componentwise |    no     |    yes    |  2   |
| SMWG         | componentwise |    no     |    yes    |  2   |
| USAMWG       | componentwise |    no     |    yes    |  2   |
| USMWG        | componentwise |    no     |    yes    |  2   |
| CHARM        | componentwise |    no     |    no     |  2   |
| HARM         | componentwise |    no     |    no     |  2   |
| Zanella      | componentwise |    no     |    yes    |  1   |
| AFSS         | slice         |    no     |    no     |  1   |
| Slice        | slice         |    no     |    no     |  2   |
| ESS          | slice         |    no     |    no     |  1   |
| OHSS         | slice         |    no     |    no     |  2   |
| RSS          | slice         |    no     |    no     |  2   |
| UESS         | slice         |    no     |    no     |  2   |
| Refractive   | slice         |    yes    |    no     |  2   |
| QSS          | slice         |    no     |    no     |  1   |
| LSS          | slice         |    no     |    no     |  2   |
| GPSS         | slice         |    no     |    no     |  2   |
| AIES         | ensemble      |    no     |    no     |  1   |
| DEMC         | ensemble      |    no     |    no     |  2   |
| t-walk       | ensemble      |    no     |    no     |  2   |
| MEADS        | ensemble      |    yes    |    no     |  2   |
| Zeus         | ensemble      |    no     |    no     |  1   |
| DREAM        | ensemble      |    no     |    no     |  1   |
| QMC          | qmc           |    no     |    no     |  2   |
| QMC3         | qmc           |    no     |    no     |  2   |
| QMCN         | qmc           |    no     |    no     |  2   |
| DQMC         | qmc           |    no     |    no     |  2   |
| SQMC         | qmc           |    no     |    no     |  2   |
| MAMC         | qmc           |    no     |    no     |  2   |
| SAMC         | qmc           |    no     |    no     |  2   |
| WMC          | qmc           |    no     |    no     |  2   |
| RMHMC        | geometric     |    yes    |    no     |  1   |
| LMC          | geometric     |    yes    |    no     |  2   |
| MHMC         | geometric     |    yes    |    no     |  2   |
| BPS          | pdmp          |    yes    |    no     |  2   |
| Boomerang    | pdmp          |    yes    |    no     |  2   |
| ZigZag       | pdmp          |    yes    |    no     |  2   |
| RHMC         | pdmp          |    yes    |    no     |  2   |
| SimTemp      | multimodal    |    no     |    no     |  1   |
| NRST         | multimodal    |    no     |    no     |  1   |
| WL           | multimodal    |    no     |    no     |  2   |
| NRPT         | multimodal    |    no     |    no     |  1   |
| ProjLang     | constraint    |    yes    |    no     |  1   |
| ProxMCMC     | constraint    |    yes    |    no     |  1   |
| NeuTra       | flow          |    yes    |    no     |  2   |
| flowMC       | flow          |    no     |    no     |  2   |
| DA           | surrogate     |    no     |    no     |  2   |
| RJ           | other         |    no     |    yes    |  2   |
| PMCMC        | other         |    no     |    no     |  2   |
| PG           | augmentation  |    no     |    no     |  2   |

Complete MCMC algorithm reference. Tier: 1 = recommended, 2 = standard,
3 = legacy/specialized. {.table}
