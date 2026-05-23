# Interoperability: bridging Stan, JAGS, and the R posterior ecosystem in lucifer

## The landscape of Bayesian computation in R

The R ecosystem for Bayesian computation has evolved through three
generations of probabilistic programming. The first generation began
with the BUGS language (Bayesian inference Using Gibbs Sampling),
created by the MRC Biostatistics Unit at Cambridge in 1989
[\[4\]](#ref4). WinBUGS (1997) made this accessible via a Windows GUI
[\[6\]](#ref6), and its open-source successor OpenBUGS (2004) extended
the reach further. JAGS (Just Another Gibbs Sampler), written by Martyn
Plummer in 2007 [\[7\]](#ref7), broke free of the Windows dependency and
introduced a clean, portable C++ implementation that remains widely used
today via the `rjags` and `runjags` R packages.

The second generation arrived with Stan in 2012, developed by Andrew
Gelman’s group at Columbia University [\[2\]](#ref2). Stan replaced
Gibbs sampling with Hamiltonian Monte Carlo (HMC) and its adaptive
variant, the No-U-Turn Sampler (NUTS) [\[5\]](#ref5). By compiling
models to C++ and exploiting automatic differentiation, Stan achieved
dramatic improvements in sampling efficiency for continuous parameter
spaces. The `rstan` package brought Stan to R, and `brms` (2017)
[\[1\]](#ref1) layered an lme4-style formula interface on top, making
Stan accessible to applied researchers who had never written a
probabilistic program.

lucifer occupies a distinct position in this landscape. It is the
next-generation refactor of LaplacesDemon, a library originally released
in 2010 by Statisticat LLC, removed from CRAN in 2015, and recently
rescued in a different
[branch](https://github.com/LaplacesDemonR/LaplacesDemon) with little
changes from the original. The refactor preserved the original design
principle (models are R functions, data are R lists, the entire
inference pipeline runs inside the R session) but rebuilt the
computational core in modern C++ via Rcpp and RcppArmadillo with OpenMP
parallelization, expanded the algorithm catalogue from a handful to
**over 130 inference engines** (82 MCMC samplers, 8 variational Bayes
methods including Pathfinder, 18 Laplace optimizers, 3 iterative
quadrature rules, population and sequential Monte Carlo, 4 ABC variants,
6 simulation-based inference methods, and 5 state-space model engines),
and added a declarative
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL that compiles models to native C++ for zero R-callback overhead at
sampling time. On top of these inference engines, lucifer ships an
automated Prescribe / Consort / Arena / Crucible pipeline that selects
algorithms and diagnoses convergence, the RobustBayes
sensitivity-analysis suite, LOO-PSIS / K-fold / leave-future-out
cross-validation with Bayesian stacking, a frequentist-Bayesian bridge,
and domain-specific modules for stochastic differential equations,
neural ODEs, and integrated population models. The package absorbed
everything LaplacesDemon ever did and added a decade of subsequent
methodological progress on top.

The interoperability bridge described in this vignette connects these
three ecosystems. Rather than competing with Stan on HMC-NUTS or JAGS on
Gibbs, lucifer absorbs their output: users can fit models with their
preferred engine and then apply lucifer’s post-processing toolkit,
including
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
for sensitivity analysis,
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
for cross-validation,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
for diagnostics, and the full suite of visualization and summary
functions.

## Architecture: how Stan, JAGS, and lucifer differ

The three systems make fundamentally different design choices about
model representation, compilation strategy, and sampler architecture.

**JAGS** uses an interpreted, graph-based approach. The user writes a
model in the BUGS language, specifying stochastic and deterministic
nodes. JAGS builds a directed acyclic graph (DAG) from this
specification and exploits conjugacy relationships to select efficient
samplers for each node. Normal-normal, gamma-Poisson, and beta-binomial
conjugate pairs are handled by Gibbs sampling with exact conditional
distributions; non-conjugate nodes receive Metropolis-Hastings updates
or slice sampling. This graph-aware strategy means that JAGS can be
remarkably efficient for hierarchical models with many conjugate
relationships, but it struggles with highly correlated posteriors or
models that lack conjugacy.

**Stan** takes the opposite approach. Models are compiled to C++ code
that evaluates the log-posterior density and its gradient via
reverse-mode automatic differentiation. The NUTS sampler then explores
the posterior using Hamiltonian dynamics, adapting the step size and
mass matrix during warmup. Stan’s strength is its generality: any model
with a differentiable log-density can be sampled efficiently, regardless
of conjugacy. The price is compilation time (30-60 seconds for a new
model) and the restriction to continuous parameters (discrete parameters
must be marginalized out).

**lucifer** occupies a middle ground. Models are R functions that return
the log-posterior, deviance, monitored quantities, posterior predictive
samples, and the (possibly transformed) parameter vector. This
five-element return contract is both the key constraint and the key
strength: it enables lucifer to wrap any model evaluation, including
calls to external libraries, C++ code, or even other probabilistic
programming languages. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL provides a declarative interface similar to BUGS/Stan for standard
models, while the raw R function interface handles everything else.

The interoperability bridge exploits a simple observation: once sampling
is complete, the output from all three systems reduces to the same
structure, namely a matrix of posterior samples with parameter names as
columns and iterations as rows. Converting between systems is therefore
a matter of extracting this matrix, computing standard diagnostics, and
assembling the result into the target format.

## Installing the backends

Before running the case study below, the JAGS and Stan backends must be
installed on the host system. Each ecosystem has its own toolchain
requirements, and the installation path depends on the operating system.
This section summarises the minimum setup for every backend referenced
in the vignette; readers who already have their preferred backend
installed can skip ahead.

### JAGS (system binary + rjags, runjags)

JAGS is a standalone C++ program that must be installed separately from
the R packages that interface with it. On Debian/Ubuntu use
`sudo apt install jags`; on Fedora/RHEL, `sudo dnf install jags`; on
macOS, either `brew install jags` or the sourceforge installer; on
Windows, the [SourceForge
installer](https://sourceforge.net/projects/mcmc-jags/) is the canonical
route. Once the `jags` binary is on the PATH, the two R wrappers are
ordinary CRAN packages:

``` r

install.packages(c("rjags", "runjags"))
```

`rjags` exposes the low-level
[`jags.model()`](https://rdrr.io/pkg/rjags/man/jags.model.html) and
[`coda.samples()`](https://rdrr.io/pkg/rjags/man/coda.samples.html)
pipeline, while `runjags` layers a higher-level wrapper
([`run.jags()`](https://rdrr.io/pkg/runjags/man/run.jags.html)) that
supports parallel chain execution and richer output summaries. Both
share the same underlying binary, so installing one JAGS is enough for
both.

### Stan (C++ toolchain + rstan / cmdstanr / brms / rstanarm)

Stan compiles each model to C++ and requires a working C++ toolchain on
the host. The prerequisites per platform are: **Rtools** on Windows
(match the version to your R installation, available from CRAN at
<https://cran.r-project.org/bin/windows/Rtools/>); **Xcode Command Line
Tools** on macOS (install with `xcode-select --install`, then configure
via the `rstan` quickstart guide); **build-essential** plus `libv8-dev`
on Debian/Ubuntu (`sudo apt install build-essential libv8-dev`) or
`@development-tools` on Fedora/RHEL.

Once the toolchain is in place, the four R packages that wrap Stan are
installed independently:

``` r

# Low-level Stan interface; writing Stan code by hand
install.packages("rstan")

# High-level formula interface; compiles to Stan code under the hood
install.packages("brms")

# Pre-compiled GLM/GLMM collection; no compilation at runtime
install.packages("rstanarm")

# CmdStanR: the modern, lightweight alternative to rstan
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           getOption("repos")))
# After installing the R package, fetch and compile the CmdStan binary
# (one-time, ~5 min). This is the only step that requires the full
# C++ toolchain; subsequent model fits reuse the compiled CmdStan.
cmdstanr::install_cmdstan()
```

Of the four Stan wrappers, **rstanarm** is the only one that works
without a user-level C++ toolchain, because its models are pre-compiled
at package build time. This makes it the easiest entry point for applied
users who want to fit generalised linear mixed models without writing
Stan code and without setting up Rtools or Xcode. **brms** provides a
comparably high-level interface but compiles each new model at runtime,
so it inherits the full toolchain requirement. **rstan** is the
low-level hand-written Stan route, and **cmdstanr** is a modern,
lightweight alternative to rstan with faster compilation and smaller
memory footprint, at the cost of a one-time
[`install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)
step that fetches and builds the CmdStan binary.

The `rstantools` package occasionally comes up in this context; it is a
*developer-side* scaffolding tool used to create new R packages that
embed Stan code. End users of brms, rstanarm, or rstan never need to
interact with it directly.

### Summary of backend requirements

| Backend  | System binary | Model language | Compiles at runtime |
|:---------|:--------------|:---------------|:--------------------|
| rjags    | JAGS          | BUGS           | No (interpreted)    |
| runjags  | JAGS          | BUGS           | No (interpreted)    |
| rstan    | C++ toolchain | Stan           | Yes                 |
| cmdstanr | C++ + CmdStan | Stan           | Yes                 |
| brms     | C++ toolchain | formula        | Yes                 |
| rstanarm | none          | formula        | No (pre-compiled)   |

Installation requirements per backend. rstanarm is the only
zero-toolchain option; all other Stan wrappers require
Rtools/Xcode/build-essential for model compilation. {.table}

## A case study: 8-schools in every backend and every lucifer family

This section walks through a single model fit with **every** JAGS
wrapper (`rjags`, `runjags`), **every** Stan wrapper (`rstan`, `brms`,
`rstanarm`; `cmdstanr` shown with code only because it requires the
separately installed CmdStan binary), and **every** applicable lucifer
inference family: Laplace approximation, iterative quadrature,
variational Bayes (Pathfinder), sequential Monte Carlo, population Monte
Carlo, and four distinct MCMC algorithms (adaptive random-walk,
affine-invariant ensemble, ensemble slice, differential-evolution Markov
chain). The same posterior, every engine, one vignette. The model is
Rubin’s **8-schools** [\[9\]](#ref9), the canonical example of Bayesian
hierarchical inference popularised by Gelman and colleagues in
[\[10\]](#ref10) and used as the running example in the original NUTS
paper [\[5\]](#ref5). Its funnel-shaped posterior under the centered
parameterisation is the single best pedagogical vehicle for exposing how
different samplers react to problematic geometry, and its real data
(eight educational coaching experiments with known standard errors) make
the comparison concrete.

### The model

Each of eight schools ran a coaching program; each reports an effect
estimate \\y_j\\ with its own known standard error \\\sigma_j\\. The
question is how to borrow strength across schools while respecting the
substantial between-school variation:

\\y_j \sim \mathcal{N}(\theta_j, \sigma_j), \quad j = 1, \ldots, 8
\qquad (\sigma_j \text{ known})\\ \\\theta_j \sim \mathcal{N}(\mu,
\tau)\\ \\\mu \sim \mathcal{N}(0, 10), \qquad \tau \sim
\mathcal{HC}(5)\\

The ten parameters are \\\theta_1, \ldots, \theta_8\\, \\\mu\\, and
\\\tau\\.

### The data (Rubin 1981)

``` r

library(lucifer)

J           <- 8
school_name <- c("A", "B", "C", "D", "E", "F", "G", "H")
y_obs       <- c(28, 8, -3,  7, -1,  1, 18, 12)
sigma_obs   <- c(15, 10, 16, 11,  9, 11, 10, 18)

stan_data <- list(J = J, y = y_obs, sigma = sigma_obs)
jags_data <- stan_data
df_long   <- data.frame(school = school_name,
                        y = y_obs, sigma = sigma_obs,
                        stringsAsFactors = FALSE)
```

### Centered and non-centered parameterisations

The 8-schools likelihood pulls \\\tau\\ toward zero because eight
observations carry very little information about between-school
variance. Under the **centered** parameterisation \\\theta_j \sim
\mathcal{N}(\mu, \tau)\\, this creates a funnel: when \\\tau\\ is small,
all the \\\theta_j\\ must cluster tightly around \\\mu\\, so the region
of high posterior density becomes a thin cone in 10-dimensional
parameter space. Gradient-based samplers (NUTS, HMC) report divergent
transitions on centered 8-schools; mode-finders such as Laplace
approximation and iterative quadrature converge to the pathological
point \\\tau = 0\\ where the likelihood is ill-defined.

The **non-centered** parameterisation rewrites the prior as \\\eta_j
\sim \mathcal{N}(0, 1)\\ with \\\theta_j = \mu + \tau\\ \eta_j\\. The
funnel unfolds: the joint posterior of \\(\eta_1, \ldots, \eta_8, \mu,
\tau)\\ has no singularity at \\\tau = 0\\, and gradient-based samplers
move freely. This reparameterisation is the single most important trick
in hierarchical Bayesian modelling with modern HMC; every Stan tutorial
eventually teaches it.

We fit the **centered** version with the JAGS and Stan wrappers (because
JAGS Gibbs handles it fine by exploiting conjugacy, and we want to show
what Stan’s NUTS does on it) and the **non-centered** version with
lucifer (to show that every lucifer inference family converges cleanly
once the parameterisation is correct). The final posterior-comparison
plot overlays both sets on a common scale, demonstrating that the
disagreement visible in rstan’s high \\\hat{R}\\ on centered 8-schools
is a parameterisation artefact, not a sampler failure.

### JAGS via rjags

``` r

library(rjags)

jags_code <- "
model {
  for (j in 1:J) {
    y[j] ~ dnorm(theta[j], 1 / (sigma[j] * sigma[j]))
    theta[j] ~ dnorm(mu, tau_prec)
  }
  mu ~ dnorm(0, 0.01)                   # sd = 10
  tau_prec <- pow(tau, -2)
  tau ~ dt(0, 1 / (5 * 5), 1) T(0,)     # half-Cauchy(5)
}
"

t0 <- Sys.time()
jm <- jags.model(textConnection(jags_code), data = jags_data,
                 n.chains = 3, n.adapt = 1000, quiet = TRUE)
update(jm, 2000)
rj_samples <- coda.samples(jm,
    variable.names = c("theta", "mu", "tau"),
    n.iter = 3000, thin = 2)
time_rjags <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### JAGS via runjags

The `runjags` package provides a higher-level wrapper around JAGS that
supports parallel chain execution and richer output summaries. The BUGS
code is identical to the rjags example above; only the driver changes.

``` r

library(runjags)

t0 <- Sys.time()
rj2 <- run.jags(
    model    = jags_code,
    data     = jags_data,
    monitor  = c("theta", "mu", "tau"),
    n.chains = 3,
    adapt    = 1000,
    burnin   = 2000,
    sample   = 3000,
    thin     = 2,
    method   = "rjags",
    silent.jags = TRUE,
    summarise = FALSE)
time_runjags <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### Stan via rstan (hand-written Stan code)

``` r

library(rstan)

stan_code <- "
data {
  int<lower=1> J;
  vector[J] y;
  vector<lower=0>[J] sigma;
}
parameters {
  vector[J] theta;
  real mu;
  real<lower=0> tau;
}
model {
  mu  ~ normal(0, 10);
  tau ~ cauchy(0, 5);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
generated quantities {
  vector[J] log_lik;
  for (j in 1:J)
    log_lik[j] = normal_lpdf(y[j] | theta[j], sigma[j]);
}
"

t0 <- Sys.time()
stan_fit <- stan(model_code = stan_code, data = stan_data,
                 chains = 3, iter = 2000, warmup = 1000,
                 seed = 42, refresh = 0)
time_rstan <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### Stan via brms (formula interface)

brms accepts a formula with the special `se(sigma)` term to encode the
known per-observation standard deviation, then compiles the model to
Stan under the hood. This is the highest-level way to fit 8-schools in
the Stan ecosystem.

``` r

library(brms)

brms_priors <- c(
    prior(normal(0, 10),  class = "Intercept"),
    prior(cauchy(0, 5),   class = "sd", group = "school")
)

t0 <- Sys.time()
brms_fit <- brm(
    y | se(sigma) ~ 1 + (1 | school),
    data   = df_long,
    prior  = brms_priors,
    chains = 3, iter = 2000, warmup = 1000,
    seed = 42, refresh = 0, silent = 2
)
time_brms <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### Stan via rstanarm (pre-compiled GLMM interface)

rstanarm specialises in generalised linear mixed models with a
homoscedastic residual assumption. The 8-schools model has a different
known \\\sigma_j\\ for every school, which rstanarm’s `stan_glmer`
interface does not directly support. The best rstanarm approximation is
a weighted Gaussian GLMM with `weights = 1 / sigma^2`, which is **not**
the exact 8-schools likelihood; the resulting posterior is biased toward
the inverse-variance-weighted estimate and underestimates \\\tau\\.
Including it here is deliberate: it illustrates a genuine limitation of
rstanarm’s formula interface for this model class, and the
posterior-comparison plot below will show the discrepancy visually.

``` r

library(rstanarm)

t0 <- Sys.time()
arm_fit <- stan_glmer(
    y ~ 1 + (1 | school),
    data = df_long,
    family = gaussian(),
    weights = 1 / sigma^2,
    prior_intercept = normal(0, 10),
    prior_covariance = decov(),
    chains = 3, iter = 2000, warmup = 1000,
    seed = 42, refresh = 0
)
time_rstanarm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### Stan via cmdstanr (reference only)

`cmdstanr` is the modern, lightweight Stan wrapper recommended by the
Stan development team. It reuses the hand-written Stan code from the
`rstan` block above, but requires a separately installed CmdStan binary
via
[`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)
which is not available in the vignette build environment. The code below
is shown for reference; see `cmdstanr`’s own documentation at
<https://mc-stan.org/cmdstanr/> for a full walkthrough.

``` r

library(cmdstanr)

mod <- cmdstan_model(write_stan_file(stan_code))
cmdstan_fit <- mod$sample(
    data = stan_data,
    chains = 3,
    iter_sampling = 1000,
    iter_warmup = 1000,
    seed = 42,
    refresh = 0)

fit_cmdstan <- as.demonoid(cmdstan_fit)
```

### The lucifer Model function (non-centered)

lucifer accepts either a
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL specification or a raw `Model` function that returns the
five-element contract `list(LP, Dev, Monitor, yhat, parm)`. For
8-schools we use the raw function because it lets us express the
non-centered parameterisation in three lines and monitor the derived
\\\theta_j\\ alongside the sampled \\(\eta_j, \mu, \tau)\\. The
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
wrapper enforces \\\tau \> 0\\ transparently, so samplers can propose
unconstrained moves without explicit log-transform bookkeeping.

``` r

parm_names <- c(paste0("eta[", 1:J, "]"), "mu", "tau")
mon_names  <- c("LP", "mu", "tau", paste0("theta[", 1:J, "]"))

Model_8s <- function(parm, Data) {
    eta   <- parm[1:Data$J]
    mu    <- parm[Data$J + 1]
    tau   <- parm[Data$J + 2]
    tau   <- interval(tau, 1e-10, Inf)
    parm[Data$J + 2] <- tau

    # Deterministic: theta_j = mu + tau * eta_j (non-centered)
    theta <- mu + tau * eta

    # log-priors on the unconstrained (eta, mu, tau)
    lp <- sum(dnorm(eta, 0, 1, log = TRUE)) +
          dnorm(mu, 0, 10, log = TRUE) +
          dhalfcauchy(tau, 5, log = TRUE)

    # log-likelihood on the derived theta
    LL <- sum(dnorm(Data$y, theta, Data$sigma, log = TRUE))

    LP <- LL + lp
    list(LP      = LP,
         Dev     = -2 * LL,
         Monitor = c(LP = LP, mu = mu, tau = tau, theta),
         yhat    = rnorm(length(Data$y), theta, Data$sigma),
         parm    = parm)
}

Data_luc <- list(J = J, y = y_obs, sigma = sigma_obs,
                 parm.names = parm_names, mon.names = mon_names)
IV <- c(rep(0, J), mean(y_obs), 5)
names(IV) <- parm_names
```

### lucifer LaplaceApproximation (BFGS)

[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
finds the posterior mode via numerical optimisation and then draws
samples from a Gaussian approximation centered at the mode with
covariance \\(-H)^{-1}\\ where \\H\\ is the Hessian. On 8-schools
non-centered, BFGS converges slowly because the log-posterior is nearly
flat along the \\\mu\\ axis; a loose `Stop.Tolerance` and the `Identity`
covariance estimator are sufficient to produce a usable posterior
approximation.

``` r

t0 <- Sys.time()
fit_la <- LaplaceApproximation(Model_8s, parm = IV, Data = Data_luc,
                                Iterations = 2000, Method = "BFGS",
                                Stop.Tolerance = 1e-3,
                                CovEst = "Identity", sir = TRUE)
time_la <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### lucifer IterativeQuadrature (CAGH)

[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)
integrates the posterior via adaptive Gauss-Hermite quadrature. The
componentwise variant (`CAGH`) is dimension-safe and can handle the
10-parameter 8-schools model in principle, but the hierarchical funnel
remains a challenge even after reparameterisation: the quadrature nodes
must be placed carefully around the posterior mode, and CAGH’s default
tolerance schedule is too strict for this problem. On 8-schools the call
below does not converge to a formal posterior-sample draw, but the
diagnostic output is still useful; it demonstrates that iterative
quadrature is best reserved for lower-dimensional posteriors with
well-behaved marginal densities.

``` r

IV_iq <- if (!is.null(fit_la$Summary1)) {
    fit_la$Summary1[seq_along(IV), "Mode"]
} else IV
names(IV_iq) <- parm_names

t0 <- Sys.time()
fit_iq <- tryCatch(
    IterativeQuadrature(Model_8s, parm = IV_iq, Data = Data_luc,
                        Iterations = 100, Algorithm = "CAGH",
                        Stop.Tolerance = c(1e-3, 1e-10),
                        Specs = list(N = 3, Nmax = 10,
                                     Packages = NULL, Dyn.libs = NULL)),
    error = function(e) NULL)
time_iq <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### lucifer VariationalBayes (Pathfinder)

Pathfinder [\[11\]](#ref11) runs a single L-BFGS optimisation and
selects the best Gaussian approximation along the L-BFGS trajectory by
maximising the evidence lower bound (ELBO). The method is dramatically
faster than ADVI or BBVI for well-behaved posteriors and often recovers
a usable approximation in seconds. We force `Covar = 1` (single path) to
avoid the parallel-dispatch overhead on tiny models; for production use,
higher `Covar` values run multiple paths concurrently and stitch them
together via PSIS resampling.

``` r

t0 <- Sys.time()
fit_vb <- VariationalBayes(Model_8s, parm = IV, Data = Data_luc,
                           Iterations = 500, Method = "Pathfinder",
                           Covar = 1L, sir = TRUE)
time_vb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### lucifer SMC (sequential Monte Carlo)

Sequential Monte Carlo tempers the posterior from the prior (\\\beta =
0\\) to the full likelihood (\\\beta = 1\\) in a sequence of
intermediate distributions, reweighting and rejuvenating particles at
each stage. SMC is uniquely valuable because it produces a
log-marginal-likelihood estimate as a by-product, which JAGS and Stan do
not. On 8-schools the adaptive schedule skips directly to \\\beta = 1\\
because the prior is already close to the posterior (weak likelihood
with only eight observations), so we force a linear 15-stage schedule to
make the rejuvenation step meaningful.

``` r

t0 <- Sys.time()
set.seed(42)
fit_smc <- SMC(Model_8s, Data = Data_luc, Initial.Values = IV,
               Iterations = 1000, N.particles = 1500,
               ESS.threshold = 0.5, Rejuvenation = "RWM",
               Rejuvenation.steps = 15,
               Schedule = seq(0.05, 1, length.out = 15),
               Status = 9999)
time_smc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### lucifer PMC (population Monte Carlo)

Population Monte Carlo maintains a population of adaptively-tuned
importance-sampling proposals and refines them iteratively. Unlike SMC
it does not temper the target, but it shares SMC’s parallel-sample
structure and also yields a marginal-likelihood estimate. On 8-schools
PMC converges quickly with modest sample sizes.

``` r

t0 <- Sys.time()
set.seed(42)
fit_pmc <- PMC(Model_8s, Data = Data_luc, Initial.Values = IV,
               Iterations = 10, Thinning = 1, N = 1000, M = 1,
               nu = 10)
time_pmc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
```

### lucifer MCMC: four gradient-free algorithms

The
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
dispatcher exposes 82 MCMC algorithms behind a single interface. For the
case study we fit four that span the non-gradient landscape: **RAM**
(Robust Adaptive Metropolis, single-chain adaptive random walk),
**AIES** (Affine-Invariant Ensemble Sampler, the population method
popularised by emcee), **Zeus** (Ensemble Slice Sampler, gradient-free
and auto-tuning), and **DREAM** (Differential Evolution Adaptive
Metropolis, which crosses chains to generate proposals). Each is
well-suited to the 8-schools posterior once reparameterised.

``` r

run_mcmc <- function(algo, iters = 10000, chains = 3) {
    t0 <- Sys.time()
    set.seed(42)
    f <- lucifer(Model_8s, Data = Data_luc, Initial.Values = IV,
                 Iterations = iters, Status = 9999, Thinning = 5,
                 Algorithm = algo, Specs = NULL,
                 Chains = chains)
    list(fit = f,
         time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

res_ram   <- run_mcmc("RAM")
res_aies  <- run_mcmc("AIES")
res_zeus  <- run_mcmc("Zeus")
res_dream <- run_mcmc("DREAM")
```

### ABC: why not here

Approximate Bayesian computation replaces likelihood evaluation with
simulation and summary-statistic matching. It is designed for problems
where the likelihood is *intractable*: population-genetics coalescent
models, stochastic agent-based simulators, implicit ODE solutions, and
so on. For 8-schools the likelihood is an analytical Gaussian product,
so ABC is the wrong tool for this job. The code pattern is shown below
for completeness, but running it on 8-schools would produce a rough
approximation at much higher cost than any likelihood-based method.

``` r

# ABC-SMC with hand-crafted summary statistics.
# See vignette("examples") for cases where ABC shines (Wright-Fisher
# diffusion, SIR epidemics, individual-based population models).
abc_summary <- function(y) c(mean = mean(y), sd = sd(y),
                              min = min(y), max = max(y))

fit_abc <- ABC(Model = NULL,            # likelihood-free
                Summary = abc_summary,
                Data = Data_luc,
                Method = "SMC",
                Tolerance = c(10, 5, 2, 1, 0.5),
                N.particles = 1000)
```

### Unified conversion via as.demonoid()

Every external fit is funnelled into a native `demonoid` object through
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
after which the rest of lucifer’s post-processing pipeline applies
uniformly. The lucifer-native fits (`laplace`, `iterquad`, `vb`, `smc`,
`pmc`, `demonoid`) already live in lucifer’s own class hierarchy, so no
conversion is needed.

``` r

fit_rjags    <- as.demonoid(rj_samples)
fit_runjags  <- as.demonoid(rj2)
fit_rstan    <- as.demonoid(stan_fit)
fit_brms     <- as.demonoid(brms_fit)
fit_rstanarm <- as.demonoid(arm_fit)

fits <- list(
    rjags    = fit_rjags,
    runjags  = fit_runjags,
    rstan    = fit_rstan,
    brms     = fit_brms,
    rstanarm = fit_rstanarm,
    luc_la   = fit_la,
    luc_iq   = fit_iq,
    luc_vb   = fit_vb,
    luc_smc  = fit_smc,
    luc_pmc  = fit_pmc,
    luc_ram   = res_ram$fit,
    luc_aies  = res_aies$fit,
    luc_zeus  = res_zeus$fit,
    luc_dream = res_dream$fit
)
fits <- Filter(Negate(is.null), fits)

times_s <- c(
    rjags    = time_rjags, runjags  = time_runjags,
    rstan    = time_rstan, brms     = time_brms,
    rstanarm = time_rstanarm,
    luc_la   = time_la,    luc_iq   = time_iq,
    luc_vb   = time_vb,    luc_smc  = time_smc,
    luc_pmc  = time_pmc,
    luc_ram   = res_ram$time, luc_aies  = res_aies$time,
    luc_zeus  = res_zeus$time, luc_dream = res_dream$time
)
```

### Diagnostic comparison via Consort()

[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
runs modern convergence diagnostics (split \\\hat{R}\\, bulk and tail
ESS, MCSE, divergent transitions) on any `demonoid` object, regardless
of whether it came from lucifer’s own samplers or an imported external
fit. We call it on all 14 engines and assemble a compact comparison
table.

``` r

suppressMessages({
    cx_list <- lapply(fits, function(f)
        tryCatch(Consort(f, verbose = FALSE), error = function(e) NULL))
})

get_fld <- function(cx, path, default = NA_real_) {
    tryCatch({
        v <- cx[[path[1]]][[path[2]]]
        if (is.null(v) || length(v) == 0) default else v
    }, error = function(e) default)
}

family_map <- c(
    rjags = "JAGS", runjags = "JAGS",
    rstan = "Stan", brms = "Stan", rstanarm = "Stan",
    luc_la = "lucifer", luc_iq = "lucifer", luc_vb = "lucifer",
    luc_smc = "lucifer", luc_pmc = "lucifer",
    luc_ram = "lucifer", luc_aies = "lucifer",
    luc_zeus = "lucifer", luc_dream = "lucifer"
)

consort_tbl <- data.frame(
    engine   = names(fits),
    family   = family_map[names(fits)],
    time_s   = round(unname(times_s[names(fits)]), 1),
    score    = sapply(cx_list, function(c)
                       round(get_fld(c, c("diagnostics",
                                           "convergence_score")), 1)),
    rhat_max = sapply(cx_list, function(c)
                       round(get_fld(c,
                             c("diagnostics", "rhat_max")), 3)),
    ess_min  = sapply(cx_list, function(c)
                       round(get_fld(c,
                             c("diagnostics", "ess_bulk_min")), 0)),
    n_div    = sapply(cx_list, function(c)
                       get_fld(c, c("diagnostics", "n_divergent"),
                                NA_integer_)),
    appeased = sapply(cx_list, function(c) isTRUE(c$appeased)),
    row.names = NULL
)
knitr::kable(consort_tbl, align = c("l", "l", rep("r", 4), "l", "l"),
             caption = "Convergence diagnostics across 14 engines. 'score' is Consort's 0-100 composite; 'appeased' means every criterion passed. rstan and rstanarm operate on centered / weighted formulations that stress-test the samplers, while every lucifer engine fits the non-centered reparameterisation.")
```

### NUTS-specific diagnostics on the brms fit

For fits produced by NUTS (rstan, brms, rstanarm, cmdstanr), lucifer’s
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
extracts the per-iteration divergent-transition flag, tree depth,
energy, and leapfrog-step count from `$NUTS.Diagnostics` and produces a
4-panel multipanel figure. On the centered 8-schools fit via brms, the
tree-depth panel reveals how often the NUTS integrator saturates its
trajectory-length budget, which is the characteristic funnel signature.

``` r

nd_brms <- check_nuts(fit_brms)
print(nd_brms)
plot(nd_brms)
```

### Posterior comparison plot

We extract posterior means and 95% credible intervals for \\\mu\\,
\\\tau\\, and \\\theta_1, \ldots, \theta_8\\ from every fit, harmonise
the parameter names across engines (brms uses `b_Intercept` for \\\mu\\
and `r_school[j,Intercept]` for the centered offset; rstanarm uses
`b[(Intercept) school:X]`; lucifer stores \\\theta_j\\ in the Monitor
matrix), and plot the 14 engines side by side faceted by engine family.

``` r

library(ggplot2)

extract_post <- function(name, f) {
    get_luc_theta_mutau <- function(eta_mat, mu_vec, tau_vec) {
        theta <- sapply(seq_len(J), function(j)
            mu_vec + tau_vec * eta_mat[, j])
        colnames(theta) <- paste0("theta[", 1:J, "]")
        cbind(theta, mu = mu_vec, tau = tau_vec)
    }
    if (name %in% c("rjags", "runjags", "rstan")) {
        P <- f$Posterior1
        return(P[, c(paste0("theta[", 1:J, "]"), "mu", "tau"),
                 drop = FALSE])
    }
    if (name == "brms") {
        P <- f$Posterior1
        mu <- P[, "b_Intercept"]
        tau <- P[, "sd_school__Intercept"]
        theta <- sapply(1:J, function(j)
            P[, paste0("r_school[", LETTERS[j], ",Intercept]")] + mu)
        colnames(theta) <- paste0("theta[", 1:J, "]")
        return(cbind(theta, mu = mu, tau = tau))
    }
    if (name == "rstanarm") {
        P <- f$Posterior1; cn <- colnames(P)
        intercept_col <- cn[cn == "(Intercept)"][1]
        tau_col <- grep("^Sigma\\[school:", cn, value = TRUE)[1]
        b_cols <- grep("^b\\[\\(Intercept\\) school:", cn, value = TRUE)
        b_cols <- grep("_NEW_", b_cols, invert = TRUE, value = TRUE)
        b_cols <- b_cols[order(sub(".*school:(.).*", "\\1", b_cols))]
        mu <- P[, intercept_col]
        tau <- sqrt(P[, tau_col])
        theta <- sapply(seq_along(b_cols),
                        function(i) P[, b_cols[i]] + mu)
        colnames(theta) <- paste0("theta[", 1:J, "]")
        return(cbind(theta, mu = mu, tau = tau))
    }
    # lucifer fits (non-centered): reconstruct theta from eta, mu, tau
    P <- if (!is.null(f$Posterior2) && is.matrix(f$Posterior2) &&
             nrow(f$Posterior2) > 0) f$Posterior2 else
         if (!is.null(f$Posterior1)) f$Posterior1 else
         if (!is.null(f$Posterior) && is.matrix(f$Posterior) &&
             nrow(f$Posterior) > 1) f$Posterior else NULL
    if (is.null(P)) return(NULL)
    eta_mat <- P[, paste0("eta[", 1:J, "]"), drop = FALSE]
    return(get_luc_theta_mutau(eta_mat, P[, "mu"], P[, "tau"]))
}

posts <- lapply(names(fits), function(nm) extract_post(nm, fits[[nm]]))
names(posts) <- names(fits)
posts <- Filter(Negate(is.null), posts)

pretty_params <- c("mu", "tau", paste0("theta[", 1:J, "]"))
comp_rows <- list()
for (nm in names(posts)) {
    for (p in pretty_params) {
        v <- posts[[nm]][, p]
        comp_rows[[length(comp_rows) + 1L]] <- data.frame(
            engine    = nm,
            family    = family_map[[nm]],
            parameter = p,
            mean      = mean(v),
            lo        = unname(quantile(v, 0.025)),
            hi        = unname(quantile(v, 0.975)),
            stringsAsFactors = FALSE)
    }
}
comp <- do.call(rbind, comp_rows)
comp$parameter <- factor(comp$parameter, levels = rev(pretty_params))
comp$family    <- factor(comp$family, levels = c("JAGS", "Stan", "lucifer"))

ggplot(comp, aes(x = mean, y = parameter, colour = engine)) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0, position = position_dodge(width = 0.7),
                  linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.7), size = 1.8) +
    facet_wrap(~ family, ncol = 1, scales = "free_y") +
    labs(x = "Posterior mean with 95% credible interval",
         y = NULL, colour = "Engine",
         title = "8-schools posterior comparison across 14 engines",
         subtitle = "JAGS and Stan wrappers fit the centered parameterisation; lucifer families fit the non-centered version") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          panel.grid.major.y = element_line(linetype = "dotted"),
          strip.text = element_text(face = "bold"))
```

### LOO-PSIS cross-engine comparison

[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
applies Pareto-smoothed importance sampling to leave-one-out
cross-validation estimates. It requires pointwise log-likelihoods, which
are available natively for rstan/brms/rstanarm fits (via the
`generated quantities` block or the
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
extractor built into each package) and for lucifer MCMC fits (via the
`Monitor` matrix when the model function stores pointwise LL values). We
extract log-likelihoods from every fit that supports it and compare
expected log predictive density differences across engines.

``` r

has_ll <- sapply(fits, function(f) {
    ll <- tryCatch(log_lik(f), error = function(e) NULL)
    !is.null(ll) && is.matrix(ll) && nrow(ll) == J
})
cat("Engines with usable log_lik:",
    paste(names(fits)[has_ll], collapse = ", "), "\n")

loo_list <- lapply(fits[has_ll], function(f) {
    ll <- log_lik(f)
    tryCatch(LOO(ll), error = function(e) NULL)
})
loo_list <- Filter(Negate(is.null), loo_list)

if (length(loo_list) >= 2) {
    cmp <- lucifer::loo_compare(loo_list)
    print(cmp)
}
```

### Takeaways

Seven lessons emerge from fitting the same model fourteen ways. First,
on a well-identified hierarchical posterior the external backends and
lucifer’s own samplers converge to the same point estimates up to Monte
Carlo error, confirming that lucifer’s interoperability layer is
faithful to the source engines. Second, rstan and brms on the centered
parameterisation produce high \\\hat{R}\\ and low effective sample size
because the funnel breaks NUTS’s adaptation, which is exactly the
problem the non-centered reparameterisation solves; the lucifer fits,
all of which use the non-centered form, reach convergence scores that
are uniformly higher for comparable wall-clock. Third, rstanarm’s
`stan_glmer` cannot express the known per-observation \\\sigma_j\\
directly, and the weighted-Gaussian approximation it substitutes gives a
biased posterior; if a user’s model has observation-specific known
variances, brms (via
[`se()`](https://mc-stan.org/rstanarm/reference/se.html)) or
hand-written rstan is the correct tool. Fourth, iterative quadrature
remains challenging for hierarchical models even after
reparameterisation, because the posterior is not well-approximated by a
product of 1D marginals around the mode; IQ is better reserved for
lower-dimensional smooth posteriors where it is competitive with
Laplace. Fifth, Pathfinder VB is an excellent speed-to-quality
trade-off: in a few seconds it produces an approximation close to the
full MCMC posterior for \\\mu\\ and \\\tau\\, making it a useful
initialiser for downstream samplers. Sixth, SMC provides a
log-marginal-likelihood estimate that neither JAGS nor Stan expose,
which is critical for Bayes-factor model comparison. Seventh, the four
lucifer MCMC algorithms (RAM, AIES, Zeus, DREAM) reach comparable
effective sample size per second; on problems where one of them stalls,
switching the `Algorithm` argument is a one-line change, and
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
will point the user toward the alternative that fits best.

## Importing Stan fits

The case study above already showed the full workflow for rstan, brms,
rstanarm, and cmdstanr. This section is an API reference: for each
wrapper, the one-line
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
call that turns its native object into a lucifer `demonoid`. All code
below assumes you already have a fit in hand from the case study.

### From rstan

[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
converts a `stanfit` object (the output of
[`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html) or
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html))
into a `demonoid`. Pointwise log-likelihoods are extracted automatically
if the Stan model includes a `generated quantities` block computing
`log_lik`, as in the case-study Stan code. NUTS-specific diagnostics
(divergent transitions, tree depth, step size, energy) are pulled from
[`rstan::get_sampler_params()`](https://mc-stan.org/rstan/reference/stanfit-class.html)
and stored in `$NUTS.Diagnostics`, making
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
work out of the box.

``` r

fit <- as.demonoid(stan_fit)       # stan_fit from the case study
print(fit); summary(fit)
```

### From brms

For `brmsfit` objects, the conversion automatically extracts pointwise
log-likelihoods via
[`brms::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
which works without writing a `generated quantities` block because brms
computes log-likelihoods at post-processing time.

``` r

fit <- as.demonoid(brms_fit)       # brms_fit from the case study
lucifer::LOO(lucifer::log_lik(fit))
```

When brms and lucifer are both loaded in the same session,
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
and
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
are masked by brms’s namespace. Use explicit
[`lucifer::log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
/
[`lucifer::LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
qualifiers or attach lucifer after brms to get lucifer’s implementations
on the search path.

### From rstanarm

rstanarm fits are `stanreg` objects that internally store a `stanfit` in
`$stanfit`. The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
dispatch extracts this slot, delegates to the `stanfit` path, and
additionally pulls pointwise log-likelihoods via
[`rstanarm::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html).
The converted fit advertises itself as
`Algorithm = "Stan:NUTS (rstanarm)"`.

``` r

fit <- as.demonoid(arm_fit)        # arm_fit from the case study
summary(fit)
```

### From cmdstanr

The `CmdStanMCMC` R6 object from cmdstanr is also supported.

``` r

fit <- as.demonoid(cmdstan_fit)
summary(fit)
```

## Importing JAGS fits

As with the Stan wrappers above, the case study walked through both
`rjags` and `runjags`. This subsection is again an API reference.

### From rjags

[`rjags::coda.samples()`](https://rdrr.io/pkg/rjags/man/coda.samples.html)
returns a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html), and
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
handles the `mcmc.list` class natively.

``` r

fit <- as.demonoid(rj_samples)     # rj_samples from the case study
print(fit); summary(fit)
```

Note that JAGS uses precision parameterisation (`tau = 1/sigma^2`) for
the normal distribution. The posterior samples are in whichever
parameterisation the model uses; the bridge does not reparameterise. If
you monitored `sigma <- 1/sqrt(tau)` as a deterministic node inside your
JAGS model, the posterior for `sigma` is in standard-deviation units.

### From runjags

The `runjags` object is handled transparently via the `runjags` class
dispatch.

``` r

fit <- as.demonoid(rj2)            # rj2 from the case study
summary(fit)
```

## Importing posterior::draws objects

The `posterior` package provides a standardized format for posterior
samples used across the Stan ecosystem. All three draws formats are
supported.

``` r

library(posterior)

# draws_array: 3D array [iterations, chains, variables]
arr <- array(rnorm(1000 * 4 * 5), dim = c(1000, 4, 5),
             dimnames = list(NULL, paste0("chain:", 1:4),
                             c("alpha", "beta[1]", "beta[2]",
                               "sigma", "lp__")))
draws <- as_draws_array(arr)
fit <- as.demonoid(draws)

# draws_matrix: 2D matrix (chains merged)
mat_draws <- as_draws_matrix(draws)
fit2 <- as.demonoid(mat_draws)

# draws_df: data frame with .chain, .iteration, .draw columns
df_draws <- as_draws_df(draws)
fit3 <- as.demonoid(df_draws)
```

## Using lucifer diagnostics on imported fits

Every imported fit supports the full diagnostic toolkit. The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
conversion computes BMK stationarity assessment, effective sample sizes,
Monte Carlo standard errors, and Gelman-Rubin diagnostics (for
multi-chain fits) automatically. In the case study above,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
was called on all 14 engines to build the comparison table; the same
call works on any single imported fit.

``` r

# Modern diagnostics on the brms fit from the case study
Rhat(fit_brms$Posterior1)
ESS(fit_brms$Posterior1)

# Caterpillar plot
caterpillar.plot(fit_brms, Parms = "r_school")

# Joint density
joint.density.plot(fit_brms$Posterior1[, "b_Intercept"],
                   fit_brms$Posterior1[, "sd_school__Intercept"],
                   Title = "Joint posterior: mu vs tau")
```

For fits produced by NUTS (rstan, brms, rstanarm, cmdstanr), the 4-panel
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
plot shown in the case study (`plot(check_nuts(fit_brms))`) exposes
divergences, tree-depth saturation, leapfrog counts, and energy traces
and is usually the single most informative diagnostic.

## LOO-PSIS and WAIC with imported fits

Leave-one-out cross-validation [\[8\]](#ref8) and WAIC require pointwise
log-likelihoods. The
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
extractor provides a uniform interface across all converted fits. The
case study already ran
[`lucifer::LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`lucifer::loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
across the 14 engines; the general pattern is:

``` r

# Extract pointwise log-likelihoods from any converted fit
ll <- lucifer::log_lik(fit_brms)      # N x S matrix

# LOO-PSIS
loo_result <- lucifer::LOO(ll)
print(loo_result)

# WAIC
waic_result <- lucifer::WAIC(ll)
print(waic_result)

# Compare multiple models
loos <- list(
    rstan    = lucifer::LOO(lucifer::log_lik(fit_rstan)),
    brms     = lucifer::LOO(lucifer::log_lik(fit_brms)),
    rstanarm = lucifer::LOO(lucifer::log_lik(fit_rstanarm))
)
lucifer::loo_compare(loos)
```

When brms or the loo package is also loaded, several function names
overlap with lucifer (`log_lik`, `LOO`, `loo_compare`, `Rhat`). Use
explicit namespaces
([`lucifer::log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`lucifer::LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md))
to ensure the correct function is called. When pointwise log-likelihoods
are not available (because the Stan model lacks `generated quantities`
or JAGS does not compute them),
[`log_lik()`](https://robustecologies.github.io/lucifer/reference/log_lik.md)
will raise an informative error; in this case, pass the log-likelihood
matrix manually via the `log_lik` argument to
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md).

## RobustBayes with imported fits

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
performs five sensitivity analyses on posterior distributions. The
`compare` module works with any imported fit (it only needs posterior
samples). The `influence` module works when pointwise log-likelihoods
are available. The `power` and `conflict` modules require a callable
lucifer Model function and Data list, which are not available for
externally fitted models unless you provide them explicitly.

With the case-study `fits` list in scope, a cross-model comparison is
one call:

``` r

rb <- RobustBayes(
    list(Stan = fits$rstan, brms = fits$brms, JAGS = fits$rjags),
    modules = "compare",
    verbose = FALSE
)
summary(rb)
plot(rb)

# Observation influence on any fit with a log_lik matrix
rb_inf <- RobustBayes(fits$brms, modules = "influence", verbose = FALSE)
plot(rb_inf, type = "influence")
```

## Backend fitting: lucifer_stan()

The
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
function provides a one-call interface: write Stan code, get a demonoid
back. The compiled Stan model is stored internally, enabling
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
and
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)
cross-validation.

``` r

fit <- lucifer_stan(
    stan_code = stan_code,        # the Stan code from the case study
    data = stan_data,             # the case-study data list
    chains = 3, iter = 2000, warmup = 1000,
    log_lik = "log_lik"
)
print(fit); summary(fit)
lucifer::LOO(lucifer::log_lik(fit))
```

The stored Stan model enables refitting for cross-validation without
recompilation.

## Backend fitting: lucifer_jags()

The
[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)
function provides the same interface for JAGS models.

``` r

fit <- lucifer_jags(
    jags_code = jags_code,        # the JAGS code from the case study
    data = jags_data,
    parameters.to.save = c("theta", "mu", "tau"),
    n.chains = 3, n.iter = 20000, n.burnin = 5000, n.thin = 5
)
print(fit); summary(fit)
```

## Exporting lucifer fits

lucifer fits can be exported to coda and posterior formats for use with
external diagnostic packages.

### To coda::mcmc.list

``` r

library(coda)

# Export any demonoid to mcmc.list (e.g. the brms fit from the case study)
ml <- to_mcmc_list(fit_brms)

# Use coda diagnostics
gelman.diag(ml)
effectiveSize(ml)
autocorr.plot(ml)
```

### To posterior::draws formats

``` r

library(posterior)

# draws_array preserves chain structure
draws_arr <- to_draws_array(fit_brms)
summarise_draws(draws_arr)

# draws_matrix for flat operations
draws_mat <- to_draws_matrix(fit_brms)

# draws_df for tidyverse-style analysis
draws_df <- to_draws_df(fit_brms)
```

## Model translation (experimental)

The
[`stan_to_spec()`](https://robustecologies.github.io/lucifer/reference/stan_to_spec.md)
and
[`jags_to_spec()`](https://robustecologies.github.io/lucifer/reference/jags_to_spec.md)
functions attempt to translate model code into lucifer’s
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
notation. These translators handle simple models with standard
distributions but cannot handle custom functions, complex indexing,
matrix operations, or ODE solvers. They are intended as a starting
point, not a complete solution.

``` r

stan_code_simple <- "
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  y ~ normal(mu, sigma);
}
"

spec_code <- stan_to_spec(stan_code_simple)
cat(spec_code)
# y ~ Normal(mu, sigma)
# mu ~ Normal(0, 10)
# sigma ~ HalfCauchy(5)

# Review and use
spec <- model_spec(spec_code)
```

``` r

jags_code_simple <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu, tau)
  }
  mu ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
}
"

spec_code <- jags_to_spec(jags_code_simple)
cat(spec_code)
```

The translators warn about constructs they cannot handle and mark
untranslatable lines as comments. Always review the output before
passing it to
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

## Feature compatibility matrix

The table below shows which lucifer functions work with imported fits at
each integration tier. The tiers are cumulative: “samples only” is the
base level, “+ log_lik” adds pointwise log-likelihoods, and “+ refit”
adds the ability to re-fit the model with modified data (available via
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
and
[`lucifer_jags()`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)).

| Function               | Samples only | \+ log_lik | \+ refit_fn |
|:-----------------------|:------------:|:----------:|:-----------:|
| print()                |     Yes      |    Yes     |     Yes     |
| summary()              |     Yes      |    Yes     |     Yes     |
| plot()                 |     Yes      |    Yes     |     Yes     |
| Rhat(), ESS(), MCSE()  |     Yes      |    Yes     |     Yes     |
| caterpillar.plot()     |     Yes      |    Yes     |     Yes     |
| joint.density.plot()   |     Yes      |    Yes     |     Yes     |
| Consort()              |   Partial    |  Partial   |     Yes     |
| LOO() via log_lik()    |      No      |    Yes     |     Yes     |
| WAIC() via log_lik()   |      No      |    Yes     |     Yes     |
| loo_compare()          |      No      |    Yes     |     Yes     |
| RobustBayes(compare)   |     Yes      |    Yes     |     Yes     |
| RobustBayes(influence) |      No      |    Yes     |     Yes     |
| RobustBayes(power)     |      No      |     No     |     No      |
| RobustBayes(conflict)  |      No      |     No     |     No      |
| predict()              |      No      |     No     |     No      |
| Kfold()                |      No      |     No     |     Yes     |
| LFO()                  |      No      |     No     |     Yes     |

Feature compatibility by integration tier {.table}

And the second table shows what integration tier each backend reaches
out of the box, so you can pick the right wrapper for the features you
need.

| Backend | log_lik (auto) | NUTS diag. | Refit path | Notes |
|:---|:--:|:--:|:---|:---|
| rstan | if generated quantities | yes | via lucifer_stan() | Full control; write Stan code by hand |
| brms | yes | yes | via recompile | Highest-level formula interface; compiles at runtime |
| rstanarm | yes | yes | no (pre-compiled only) | Pre-compiled; zero toolchain; limited model class |
| cmdstanr | if generated quantities | yes | via cmdstanr::cmdstan_model | Lightweight modern Stan wrapper; needs CmdStan binary |
| rjags | manual | n/a | via lucifer_jags() | Low-level BUGS interface |
| runjags | manual | n/a | via lucifer_jags() | Parallel JAGS runs; richer output than rjags |
| posterior::draws | manual | manual | no | Generic draws format; no model or log_lik |

Per-backend integration summary. ‘log_lik (auto)’ is whether pointwise
log-likelihoods are extracted automatically by as.demonoid(); ‘NUTS
diag.’ is whether divergent transitions, tree depth, and energy are
extracted into \$NUTS.Diagnostics. {.table}

## Advantages and limitations

The interoperability bridge makes lucifer a natural companion to Stan
and JAGS rather than a competitor. Users gain access to 55 MCMC
algorithms,
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
for comprehensive sensitivity analysis,
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
for automated algorithm comparison,
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
for declarative model specification, and the full diagnostic and
visualization toolkit, all without abandoning their existing modeling
workflows.

The bridge does not replicate what Stan and JAGS do well. Stan’s static
typing, GPU support, ODE solvers, Gaussian process primitives, and
algebraic equation solvers are features that lucifer does not attempt to
match. JAGS’ automatic conjugacy detection and the resulting efficiency
for large hierarchical models with natural conjugate structure is
similarly beyond the bridge’s scope. The bridge is a tool for
post-processing and cross-pollination, not a replacement for either
system.

The model translation functions
([`stan_to_spec()`](https://robustecologies.github.io/lucifer/reference/stan_to_spec.md)
and
[`jags_to_spec()`](https://robustecologies.github.io/lucifer/reference/jags_to_spec.md))
are explicitly experimental and limited to simple models. Complex Stan
programs with custom functions, matrix calculus, or ODE integration
should be used via
[`lucifer_stan()`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)
rather than translated to
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

## References

**\[1\]** Bürkner, P.-C. (2017). *brms: an R package for Bayesian
multilevel models using Stan.* Journal of Statistical Software, 80(1),
1-28. [doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

**\[2\]** Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich,
B., Betancourt, M., Brubaker, M., Guo, J., Li, P., and Riddell, A.
(2017). *Stan: a probabilistic programming language.* Journal of
Statistical Software, 76(1), 1-32.
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

**\[3\]** Gelman, A. and Rubin, D. B. (1992). *Inference from iterative
simulation using multiple sequences.* Statistical Science, 7(4),
457-472.
[doi:10.1214/ss/1177011136](https://doi.org/10.1214/ss/1177011136)

**\[4\]** Gilks, W. R., Thomas, A., and Spiegelhalter, D. J. (1994). *A
language and program for complex Bayesian modelling.* The Statistician,
43(1), 169-177. [doi:10.2307/2348941](https://doi.org/10.2307/2348941)

**\[5\]** Hoffman, M. D. and Gelman, A. (2014). *The No-U-Turn Sampler:
adaptively setting path lengths in Hamiltonian Monte Carlo.* Journal of
Machine Learning Research, 15(1), 1593-1623.
<https://jmlr.org/papers/v15/hoffman14a.html>

**\[6\]** Lunn, D. J., Thomas, A., Best, N., and Spiegelhalter, D.
(2000). *WinBUGS – a Bayesian modelling framework: concepts, structure,
and extensibility.* Statistics and Computing, 10(4), 325-337.
[doi:10.1023/A:1008929526011](https://doi.org/10.1023/A:1008929526011)

**\[7\]** Plummer, M. (2003). *JAGS: a program for analysis of Bayesian
graphical models using Gibbs sampling.* Proceedings of the 3rd
International Workshop on Distributed Statistical Computing, Vienna,
Austria. <https://www.r-project.org/conferences/DSC-2003/>

**\[8\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). *Practical
Bayesian model evaluation using leave-one-out cross-validation and
WAIC.* Statistics and Computing, 27(5), 1413-1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[9\]** Rubin, D. B. (1981). *Estimation in parallel randomized
experiments.* Journal of Educational Statistics, 6(4), 377-401.
[doi:10.3102/10769986006004377](https://doi.org/10.3102/10769986006004377)

**\[10\]** Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
Vehtari, A., and Rubin, D. B. (2013). *Bayesian Data Analysis*, Third
Edition. Chapman & Hall/CRC, Boca Raton. ISBN: 978-1-4398-4095-5

**\[11\]** Zhang, L., Carpenter, B., Gelman, A., and Vehtari, A. (2022).
*Pathfinder: parallel quasi-Newton variational inference.* Journal of
Machine Learning Research, 23(306), 1-49.
<https://jmlr.org/papers/v23/21-0889.html>
