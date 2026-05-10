# Introduction to lucifer: a complete Bayesian inference workflow

## What is lucifer?

The **lucifer** package is the modern continuation of LaplacesDemon, one
of the earliest and most ambitious R packages for general-purpose
Bayesian inference. Where most Bayesian software constrains the user to
a fixed model grammar (as in BUGS, JAGS, or Stan), lucifer takes a
fundamentally different approach: the user writes a plain R function
that evaluates the log-posterior for any parameter vector, and the
package provides the computational machinery to explore that posterior.
This design places no restrictions on the model structure. Any
probability model that can be expressed as an R function, including
hierarchical specifications, mixture models, state-space systems,
stochastic differential equations, and models with custom likelihoods
that have no closed-form density, can be fitted without modification to
the interface.

The package has grown substantially beyond its origins. The current
release implements over 130 inference algorithms spanning eight
families: 82 MCMC samplers organized across 13 subcategories
(gradient-based, adaptive random-walk, componentwise, slice, ensemble,
quasi-Monte Carlo, geometric, piecewise-deterministic, multimodal
tempering, constraint-handling, flow-enhanced, surrogate, and
augmentation methods), 8 variational Bayes algorithms, 18 Laplace
optimizers, 3 iterative quadrature methods, population Monte Carlo,
sequential Monte Carlo, 4 approximate Bayesian computation methods, and
6 simulation-based inference methods. A C++ backend with OpenMP
parallelization accelerates core computations, and a model specification
DSL
([`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md))
allows declarative model building that compiles to native C++ for
additional speed.

The unifying contract across all of these methods is a single Model
function with signature `Model(parm, Data)` that returns a list
containing five elements: `LP` (the log-posterior), `Dev` (the deviance,
equal to \\-2 \times\\ the log-likelihood), `Monitor` (quantities to
track), `yhat` (fitted values for posterior predictive checks), and
`parm` (the parameter vector after any constraint enforcement). This
contract means that every diagnostic, comparison, and visualization tool
in the package works identically regardless of which inference engine
produced the fit.

The package name continues a historical thread. Pierre-Simon Laplace
independently developed Bayes’ theorem in 1774, and his “Essai
Philosophique sur les Probabilites” (1814) envisioned a vast intellect
that could, given complete knowledge of the present state of the
universe, predict all future states. Later biographers called this
hypothetical intellect Lucifer, the light-bearer. The package aspires to
that same goal: bringing computational light to high-dimensional
posterior distributions that would otherwise remain opaque.

This vignette is a guided tour through a complete Bayesian workflow.
Starting from raw data with known ground truth, it demonstrates how to
specify models both declaratively and manually, select algorithms, fit
multiple competing models, diagnose convergence, compare models via
information criteria and cross-validation, assess sensitivity to prior
and likelihood assumptions, test identifiability through data cloning,
extend to hierarchical structures, and bridge to frequentist inference.
Each section introduces one component of the workflow; together they
form a template that applies to virtually any applied problem.

### The model contract

Every lucifer model is an R function with two arguments: `parm` (a
numeric vector of parameter values) and `Data` (a named list containing
observations and metadata). The function must return a named list with
five elements.

The `LP` element contains the log-posterior, which is the sum of the
log-likelihood and the log-prior. The `Dev` element contains the
deviance, defined as \\-2\\ times the log-likelihood. This pair allows
automatic decomposition into log-prior and log-likelihood without the
user writing separate functions for each: \\\log \pi(\theta) =
\text{LP} + \text{Dev}/2\\ and \\\ell(\theta) = -\text{Dev}/2\\. The
`Monitor` element tracks user-specified quantities across iterations,
such as derived parameters or predictive summaries. The `yhat` element
provides fitted values (predictions at the observed covariate values)
used for posterior predictive checks. The `parm` element returns the
parameter vector after any constraint enforcement, which is necessary
because some algorithms propose values outside the feasible region and
rely on the Model function to project them back.

The `Data` list must contain two character vectors: `mon.names`, which
names the monitored quantities corresponding to `Monitor`, and
`parm.names`, which names the parameters corresponding to `parm`. These
metadata entries allow the package to label output, diagnostics, and
plots without additional user input.

This contract is the single design decision that makes the entire
ecosystem possible. Because every model, whether a two-parameter linear
regression or a thousand-parameter state-space model, obeys the same
interface, every algorithm, diagnostic, and visualization tool works
without modification across all models.

## Simulating data with known ground truth

A well-designed analysis begins with a problem where the answer is
known, so that estimation accuracy can be assessed before moving to real
data. The scenario here involves linear regression with heavy-tailed
errors, a situation common in environmental and ecological applications
where extreme observations arise naturally but should not dominate
inference. The data-generating process uses a Student-\\t\\ distribution
with \\\nu = 4\\ degrees of freedom, which produces heavier tails than
the Gaussian but retains finite variance.

``` r

library(lucifer)
library(ggplot2)

set.seed(666)
N <- 300
x <- rnorm(N)

# Ground truth
beta0_true <- 2.0
beta1_true <- 1.5
sigma_true <- 1.2
nu_true <- 4

mu_true <- beta0_true + beta1_true * x
y <- mu_true + sigma_true * rt(N, df = nu_true)
```

The true regression line and the simulated observations can be
visualized together. The heavy tails of the Student-\\t\\ distribution
produce a handful of observations that deviate substantially from the
mean function, precisely the kind of data that will distinguish robust
models from those that assume Gaussian errors.

``` r

df_sim <- data.frame(x = x, y = y, mu = mu_true)

ggplot(df_sim, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, size = 1.5, color = "grey40") +
    geom_line(aes(y = mu), color = "#2166AC", linewidth = 1) +
    labs(x = "x", y = "y",
         title = "Simulated data with linear mean and Student-t errors") +
    theme_minimal()
```

The blue line is the true mean function \\\mu = 2 + 1.5x\\. Several
observations lie far from this line, reflecting the heavy tails of the
\\t_4\\ distribution. These outliers will be informative for model
comparison: a Gaussian error model will be pulled toward them, inflating
its variance estimate, while a Student-\\t\\ model will accommodate them
through its tail parameter.

The choice of this data-generating process is deliberate. Linear
regression is simple enough that the true parameters are unambiguously
identifiable, yet the heavy-tailed errors create a genuine modeling
decision: the analyst must choose whether to use a Gaussian (risking
bias from outliers), a Student-\\t\\ (requiring estimation of the tail
parameter \\\nu\\), or some other robust distribution. This tension
between model complexity and robustness is ubiquitous in applied
statistics, and the workflow demonstrated here applies directly to more
complex problems where the same tradeoff arises.

The Student-\\t\\ distribution with \\\nu = 4\\ has finite moments up to
order 3 but infinite kurtosis. The empirical kurtosis of the simulated
data will therefore vary substantially across realizations, unlike the
Gaussian case where the empirical kurtosis concentrates rapidly around
3. This excess variability in the tails is precisely what the model
comparison framework must detect: the Student-\\t\\ model can
accommodate it through the \\\nu\\ parameter, while the Gaussian model
cannot.

## Building models with model_spec()

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
function provides a declarative interface for probabilistic model
specification. Rather than writing the Model function, Data list, and
initial values by hand, the user writes a compact probabilistic program
and the system compiles it into all three components. Three competing
models illustrate the approach: one with Gaussian errors (deliberately
misspecified), one with Student-\\t\\ errors (matching the
data-generating process), and one with Laplace errors (an alternative
heavy-tailed specification).

### Model A: normal errors (misspecified)

The simplest specification assumes Gaussian errors, which will
underestimate the frequency of extreme observations and, as a
consequence, overestimate the residual variance to compensate.

``` r

spec_normal <- model_spec("
    y ~ Normal(mu, sigma)
    mu = beta0 + beta1 * x
    beta0 ~ Normal(0, 10)
    beta1 ~ Normal(0, 10)

    sigma ~ HalfCauchy(2.5)
")
```

Each line is either a stochastic statement (containing `~`) or a
deterministic assignment (containing `=`). The parser distinguishes data
variables from parameters by checking which symbols appear on the
left-hand side of stochastic statements and which appear only on the
right-hand side without being assigned. The priors follow the weakly
informative recommendations of Gelman [\[1\]](#ref1): \\\text{Normal}(0,
10)\\ for regression coefficients covers plausible effect sizes without
allowing the sampler to waste time in implausible regions, and
\\\text{HalfCauchy}(2.5)\\ for the scale parameter provides gentle
regularization.

The prior on \\\sigma\\ deserves comment. The half-Cauchy distribution
places substantial mass near zero while maintaining heavy tails, which
means it regularizes toward small variances without imposing a hard
upper bound. The scale parameter of 2.5 is the standard recommendation
for residual standard deviations in regression settings: it comfortably
covers most plausible values while discouraging extreme scales that
would otherwise slow convergence. If the data strongly support a large
\\\sigma\\, the posterior will follow the data. This is the defining
property of a weakly informative prior: it rules out absurd values
(negative variances) without pulling the posterior toward any particular
plausible value.

### Model B: Student-t errors (correctly specified)

The Student-\\t\\ model adds a degrees-of-freedom parameter \\\nu\\,
which the data will inform. The reparameterization `nu = nu_raw + 2`
ensures \\\nu \> 2\\, guaranteeing a finite variance, and the
Exponential(0.1) prior on `nu_raw` with mean 10 gently regularizes
toward moderate tail behavior while allowing the data to push \\\nu\\
toward heavier tails if warranted.

``` r

spec_t <- model_spec("
    y ~ StudentT(mu, sigma, nu)
    mu = beta0 + beta1 * x
    beta0 ~ Normal(0, 10)
    beta1 ~ Normal(0, 10)

    sigma ~ HalfCauchy(2.5)
    nu_raw ~ Exponential(0.1)
    nu = nu_raw + 2
")
```

### Model C: Laplace errors (alternative heavy-tailed)

The Laplace (double exponential) distribution provides another
heavy-tailed alternative. Its tails decay exponentially rather than
polynomially, making it heavier than the Gaussian but lighter than the
Student-\\t\\ for small \\\nu\\. It serves as a useful comparison point:
if the data truly come from a \\t_4\\ distribution, the Laplace model
should fit better than the Gaussian but worse than the Student-\\t\\.

``` r

spec_laplace <- model_spec("
    y ~ Laplace(mu, sigma)
    mu = beta0 + beta1 * x
    beta0 ~ Normal(0, 10)
    beta1 ~ Normal(0, 10)

    sigma ~ HalfCauchy(2.5)
")
```

### Inspecting and compiling specifications

The `model_spec` object provides several methods for inspection.
Printing shows the parsed structure with parameter roles, prior
distributions, and detected constraints. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
method reveals the generated R source code for the Model function, which
is useful both for learning and for extracting a starting point when the
model needs manual modification.

``` r

print(spec_t)
code(spec_t)
plot(spec_t)
```

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
renders a directed acyclic graph (DAG) showing the dependency structure
among parameters, deterministic transformations, and observed data. This
graphical representation is particularly valuable for hierarchical
models where the conditional independence structure is not immediately
obvious from the specification text. Parameters are shown as open
circles, observed data as shaded squares, and deterministic
transformations as diamonds, following standard plate notation
conventions.

### Preparing data and initial values

Each specification object provides two helper methods:
`$data_template()` binds observed data into the list format required by
the Model function, automatically adding `mon.names`, `parm.names`, and
any indexing variables; `$initial_values()` generates a plausible
starting point by sampling from the prior or using default values that
respect parameter constraints.

``` r

# Student-t model
Data_t <- spec_t$data_template(y = y, x = x)
IV_t <- spec_t$initial_values(Data_t)

# Normal model
Data_normal <- spec_normal$data_template(y = y, x = x)
IV_normal <- spec_normal$initial_values(Data_normal)

# Laplace model
Data_laplace <- spec_laplace$data_template(y = y, x = x)
IV_laplace <- spec_laplace$initial_values(Data_laplace)
```

### Compiling to C++

For computationally demanding models, the specification can be compiled
to C++ for substantially faster log-posterior evaluation. The
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
function generates C++ source code from the intermediate representation,
compiles it via Rcpp, and returns a new specification object whose
`$Model` function calls the native code. The compiled model evaluates
the log-posterior without any R interpreter overhead, which can provide
speedups of 10-100x for models where the likelihood evaluation dominates
the cost. The compilation is cached, so subsequent calls to
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
with the same specification reuse the compiled object.

The compiled model obeys the same contract as the R-level model: it
accepts the same `parm` and `Data` arguments, returns the same
five-element list, and is interchangeable in every function that accepts
a Model. This means that
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
is a pure performance optimization with no behavioral changes.

``` r

cspec_t <- compile(spec_t, Data_t)
```

### Prior predictive simulation

Before fitting the model, it is worth checking that the prior produces
predictions that are at least plausible. Prior predictive simulation
draws parameter values from the prior distributions and generates data,
without conditioning on the observations. If the prior predictive
distribution places substantial mass on impossible or absurd outcomes
(negative concentrations, probabilities outside \\\[0,1\]\\, counts in
the billions), the prior is too vague and should be tightened.

The
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
function automates this process. It takes the Model function, the Data
list, and an `rprior` function that draws samples from the joint prior,
then overlays the prior predictive distribution against the observed
data. The `rprior` function must accept a single integer `n` and return
an \\n \times K\\ matrix where each row is a draw from the joint prior
and each column corresponds to a parameter in `parm.names` order.

``` r

# Define a function that draws from the joint prior
rprior_t <- function(n) {
    cbind(
        beta0   = rnorm(n, 0, 10),
        beta1   = rnorm(n, 0, 10),
        sigma   = abs(rcauchy(n, 0, 2.5)),
        nu_raw  = rexp(n, 0.1)
    )
}

# Density overlay: prior predictive y vs. observed y
prior_predictive_check(spec_t$Model, Data_t, rprior = rprior_t, n = 200)
```

The density overlay shows the distribution of \\\tilde{y}\\ generated
under the prior (grey) against the observed data (dark line). If the
observed density sits comfortably within the prior predictive envelope,
the prior is compatible with the data without being overly concentrated.
With \\\text{Normal}(0, 10)\\ priors on the regression coefficients, the
prior predictive distribution is wide enough to cover the data but does
not extend to absurd values, confirming that the priors are weakly
informative rather than vague.

A complementary diagnostic compares a scalar test statistic computed on
the prior predictive draws against its observed value. If the observed
statistic falls in the tails of the prior predictive distribution, the
prior assigns low probability to data that look like the actual
observations.

``` r

# Test statistic: prior predictive standard deviation vs. observed
prior_predictive_check(spec_t$Model, Data_t, rprior = rprior_t, n = 200,
                       type = "stat", stat_fun = sd)
```

## Building a model without model_spec()

The declarative interface handles the majority of standard models, but
some problems require the full flexibility of a hand-written Model
function: custom likelihoods without standard distributional forms,
models with discrete latent variables that require marginalization, or
performance-critical inner loops that benefit from manual optimization.
Understanding the manual approach also clarifies what
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
automates, making the user better equipped to diagnose problems when
they arise.

The Model function for the Student-\\t\\ regression illustrates the five
required return components.

``` r

Model_manual <- function(parm, Data) {
    # Extract parameters
    beta0 <- parm[1]
    beta1 <- parm[2]

    # Enforce constraints: sigma > 0, nu > 2
    sigma <- interval(parm[3], 1e-100, Inf)
    parm[3] <- sigma
    nu <- interval(parm[4], 2, Inf)
    parm[4] <- nu

    # Deterministic transformation
    mu <- beta0 + beta1 * Data$x

    # Log-likelihood
    LL <- sum(dst(Data$y, mu = mu, sigma = sigma, nu = nu, log = TRUE))

    # Log-prior
    LP <- LL +
        dnormv(beta0, 0, 10^2, log = TRUE) +
        dnormv(beta1, 0, 10^2, log = TRUE) +
        dhalfcauchy(sigma, 2.5, log = TRUE) +
        dexp(nu - 2, rate = 0.1, log = TRUE)

    # Fitted values for posterior predictive checks
    yhat <- mu

    # Monitored quantities
    Monitor <- LP

    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}
```

Several details deserve attention. The
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
function enforces constraints by reflecting values that violate them
back into the feasible region, and the constrained value is written back
into `parm` so that the sampler operates on the constrained parameter.
This write-back is essential: if the sampler proposes \\\sigma = -0.3\\,
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
reflects it to \\0.3\\, and storing the reflected value in `parm[4]`
ensures that the chain records the constrained value rather than the
infeasible proposal.

The
[`dst()`](https://robustecologies.github.io/lucifer/reference/dist.Student.t.md)
function evaluates the Student-\\t\\ density with location and scale
parameters, following the parameterization where \\\sigma\\ is a scale
(not the standard deviation of the distribution, which would be \\\sigma
\sqrt{\nu / (\nu - 2)}\\ for \\\nu \> 2\\). The
[`dnormv()`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Variance.md)
function is the lucifer equivalent of
[`dnorm()`](https://rdrr.io/r/stats/Normal.html) but parameterized by
variance rather than standard deviation, so
`dnormv(beta0, 0, 10^2, log = TRUE)` specifies \\\beta_0 \sim N(0,
10^2)\\, equivalent to `dnorm(beta0, 0, 10, log = TRUE)`.

The log-posterior `LP` is the sum of the log-likelihood and the
log-prior; the deviance `Dev` equals \\-2 \times\\ the log-likelihood,
following the convention used throughout the package. This decomposition
is not merely bookkeeping: it enables the package to automatically
separate the data contribution from the prior contribution in downstream
analyses such as
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md),
and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md).
The `yhat` vector of fitted values enables posterior predictive checks,
and `Monitor` tracks any quantities the user wants to inspect across
iterations.

The Data list for the manual approach must include `mon.names` and
`parm.names`, which are character vectors naming the monitored
quantities and parameters, respectively.

``` r

Data_manual <- list(
    x = x, y = y, N = N,
    mon.names = "LP",
    parm.names = c("beta0", "beta1", "sigma", "nu")
)
IV_manual <- c(0, 0, 1, 5)
```

The tradeoffs between the two approaches are clear. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
route is declarative, handles constraints and transforms automatically,
compiles to C++, and generates initial values from the prior. The manual
route is fully flexible, supports arbitrary R code in the likelihood,
and gives the user complete control over parameterization. The two are
not mutually exclusive: `code(spec)` extracts the generated source,
which can then be modified by hand for models that start as standard
specifications but require custom extensions.

A practical guideline: start with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
for any model that can be expressed as a directed graphical model with
standard distributions. Switch to manual specification when the model
involves custom likelihoods (e.g., pseudo-likelihoods, composite
likelihoods), discrete latent variables that must be marginalized out,
complex data transformations within the likelihood, or non-standard
parameterizations that the DSL does not support. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
output from
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
provides a well-structured starting point even in these cases, saving
the effort of writing the boilerplate from scratch.

## Algorithm selection with Prescribe()

Before committing to a potentially expensive inference run,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
profiles the model by evaluating it a small number of times and
characterizing the problem along seven axes: dimensionality, evaluation
speed, gradient availability, constraint structure, posterior
conditioning, multimodality risk, and correlation structure. Based on
this profile, it scores all available inference methods and returns a
ranked recommendation.

``` r

rx <- Prescribe(spec_t$Model, Data_t, IV_t)
print(rx)
summary(rx)
```

The printed output shows the top-ranked methods across all inference
families, along with the model profile (a radar chart is available via
`plot(rx)`). For this four-parameter regression with a fast model
evaluation and available gradients,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
will typically recommend NUTS or one of the other gradient-based MCMC
methods as the primary choice, with ADVI or Pathfinder as a fast
variational alternative. The profile also flags potential issues: if the
condition number of the local Hessian is large, it warns about
ill-conditioning; if the multimodality score is elevated, it suggests
tempering methods.

The recommendation is advisory, not prescriptive. An analyst who knows
that the posterior is unimodal and well-conditioned might proceed
directly with NUTS; one who suspects multimodality might instead choose
NRST (non-reversible simulated tempering) or parallel tempering. The
value of
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
lies in surfacing properties of the model that would otherwise require
manual experimentation to discover.

``` r

plot(rx)
```

The radar chart displays the model profile across all seven axes,
providing an immediate visual summary of the problem’s computational
character. A model that is high-dimensional, slow to evaluate, and
poorly conditioned will show a very different profile from a
low-dimensional, fast, well-conditioned model, and the recommended
algorithms will differ accordingly. The radar chart is particularly
useful when comparing the profiles of different models applied to the
same data: if two models have similar profiles, the same algorithm will
likely work well for both.

## Fitting the models

With three model specifications prepared and an algorithm selected, the
models can be fitted using
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
The function signature mirrors the manual workflow: it takes a Model
function, Data list, and initial values, along with control parameters
for the number of iterations, thinning interval, status reporting
frequency, and the inference algorithm. Multi-chain parallelism is
controlled by the `Chains` and `CPUs` arguments; running four chains in
parallel provides both faster sampling and the between-chain diagnostics
needed for convergence assessment.

``` r

set.seed(1150)
fit_normal <- lucifer(spec_normal$Model, Data_normal, IV_normal,
                      Iterations = 20000, Status = 5000, Thinning = 10,
                      Algorithm = "NUTS", Chains = 4, CPUs = 4)

set.seed(1150)
fit_t <- lucifer(spec_t$Model, Data_t, IV_t,
                 Iterations = 20000, Status = 5000, Thinning = 10,
                 Algorithm = "NUTS", Chains = 4, CPUs = 4)

set.seed(1150)
fit_laplace <- lucifer(spec_laplace$Model, Data_laplace, IV_laplace,
                       Iterations = 20000, Status = 5000, Thinning = 10,
                       Algorithm = "NUTS", Chains = 4, CPUs = 4)
```

Note that `fit_t` uses the compiled C++ model (`spec_t$Model`) while the
other two use the R-level models. In practice, all three could be
compiled for speed, but mixing compiled and interpreted models in the
same analysis is perfectly valid since the interface contract is
identical.

Each model is fitted independently, so the three calls could be run in
parallel if the system has sufficient resources. In practice, since each
call already uses four cores for parallel chains, running them
sequentially avoids CPU oversubscription.

The `Iterations` argument specifies the total number of MCMC iterations
including the warmup (adaptation) phase. The `Status` argument controls
how frequently progress is printed to the console. The `Thinning`
argument retains every \\k\\th sample to reduce autocorrelation in the
stored output and memory usage; `Thinning = 10` stores 1000 samples per
chain from the 10,000 iterations. With `Chains = 4` and `CPUs = 4`, each
chain runs in a separate process, and the results are combined into a
single fit object with full multi-chain diagnostics.

After sampling, burn-in is removed automatically: the BMK stationarity
diagnostic partitions the thinned samples into a transient (burn-in)
phase and a stationary phase. The returned object contains two posterior
matrices: `Posterior1` with all thinned samples and `Posterior2` with
only the stationary samples. All inference quantities (`Summary2`,
`DIC2`, [`predict()`](https://rdrr.io/r/stats/predict.html)) operate on
`Posterior2` by default. If the automatic diagnostic is unsatisfactory,
the `BurnIn` argument accepts an integer specifying how many thinned
samples to discard, bypassing BMK entirely; for example, `BurnIn = 200`
discards the first 200 thinned samples from each chain. Post-hoc
adjustment is also available via
[`deburn()`](https://robustecologies.github.io/lucifer/reference/deburn.md).

The NUTS algorithm [\[9\]](#ref9) is generally the recommended default
for continuous parameters with differentiable log-posteriors. It
adaptively selects the number of leapfrog steps by building a binary
tree until a U-turn criterion is met, eliminating the need to tune the
trajectory length. During the warmup phase, it adapts the step size to
achieve a target acceptance rate (default 0.80) and estimates the mass
matrix from the warmup samples. These adaptations make NUTS largely
self-tuning, which is why
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
frequently recommends it as the first choice.

### Interpreting the output

The [`print()`](https://rdrr.io/r/base/print.html) method provides a
concise summary: posterior means, standard deviations, MCSE (Monte Carlo
standard error), ESS (effective sample size), and \\\hat{R}\\ for each
parameter. The [`summary()`](https://rdrr.io/r/base/summary.html) method
adds quantiles, highest posterior density intervals, and diagnostic
flags.

``` r

print(fit_t)
summary(fit_t)
```

For the Student-\\t\\ model, the posterior means should be close to the
true values: \\\hat{\beta}\_0 \approx 2.0\\, \\\hat{\beta}\_1 \approx
1.5\\, \\\hat{\sigma} \approx 1.2\\, and \\\hat{\nu} \approx 4\\. The
precision of recovery depends on the sample size and the realized noise,
but with \\N = 300\\ observations the posterior should concentrate
tightly around the truth.

### Visualizing the posterior

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
fitted objects supports multiple visualization types. The default
displays a panel of trace plots and density estimates for each
parameter.

``` r

# Default: trace + density panels
plot(fit_t)

# Marginal posterior densities
plot(fit_t, type = "density")

# Trace plots for mixing assessment
plot(fit_t, type = "trace")
```

The trace plots should show well-mixed chains that explore the same
region of parameter space, with no visible trends, periodicity, or stuck
segments. The density plots should show unimodal, roughly symmetric
posteriors for the regression coefficients, with the true values falling
within the central mass.

Common pathologies to watch for include: chains that remain separated in
different regions (indicating multimodality or insufficient warmup),
long stretches of identical values (indicating sticky proposals, common
in random-walk methods), and systematic drift (indicating the chain has
not yet reached stationarity). NUTS is generally resistant to the first
two pathologies but can exhibit the third if the warmup is too short for
the adaptation to converge.

### Posterior predictive fitted curve

Beyond scalar parameter estimates, the fitted model implies a
distribution over regression curves. Plotting the posterior predictive
mean function with credible bands provides a visual assessment of fit
quality and uncertainty.

``` r

# Posterior predictive draws at the observed covariate values
ppc_t <- predict(fit_t, Model = spec_t$Model, Data = Data_t)

# The 'Style = "Fitted"' ppc plot shows the posterior predictive mean for each observation against its observed value, with 95% credible intervals. Points falling on the diagonal indicate accurate prediction.

plot(ppc_t, Style = "Fitted", Data = Data_t)
```

The dark line shows the posterior predictive median and the nested
ribbons show the 50%, 80%, and 95% pointwise credible bands; the light
line traces the observed responses, sorted along the covariate. With \\N
= 300\\ observations, the bands should be narrow around the implied
regression curve and the observed values should fall predominantly
inside the 95% band. The bands widen near the edges of the covariate
range, reflecting the reduced information available at the boundaries of
the observed data.

### Comparing estimates against ground truth

A direct comparison between posterior summaries and the known true
values provides a concrete assessment of estimation quality.

``` r

truth <- c(beta0 = beta0_true, beta1 = beta1_true,
           sigma = sigma_true, nu = nu_true)

# Extract posterior means from Summary2 and derive nu from nu_raw
raw_means <- setNames(fit_t$Summary2[, "Mean"],
                      rownames(fit_t$Summary2))
post_means <- c(beta0 = raw_means["beta0"],
                beta1 = raw_means["beta1"],
                sigma = raw_means["sigma"],
                nu    = raw_means["nu_raw"] + 2)

comparison <- data.frame(
    parameter = names(truth),
    true_value = truth,
    posterior_mean = post_means,
    bias = post_means - truth
)

knitr::kable(comparison, digits = 3,
             caption = "Posterior means vs. ground truth for the Student-t model")
```

### Prior vs. posterior

A direct visual comparison of the prior and posterior distributions
shows how much the data updated the prior belief for each parameter. The
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md)
function overlays both densities in a faceted display, with optional
ground truth reference lines. When the posterior is much narrower than
the prior, the data are informative for that parameter; when the two
distributions are similar, the data carry little information and the
inference is prior-driven.

``` r

prior_vs_posterior(
    fit_t,
    prior = list(
        beta0  = function(x) dnorm(x, 0, 10),
        beta1  = function(x) dnorm(x, 0, 10),
        sigma  = function(x) 2 * dcauchy(x, 0, 2.5) * (x > 0),
        nu_raw = function(x) dexp(x, 0.1)
    ),
    ground_truth = c(beta0 = beta0_true, beta1 = beta1_true,
                     sigma = sigma_true, nu_raw = nu_true - 2)
)
```

For all four parameters the posterior concentrates tightly relative to
the prior, confirming that with \\N = 300\\ observations the data
dominate the inference. The regression coefficients show the most
dramatic contraction, as expected: the linear structure in the data pins
them down precisely. The scale parameter \\\sigma\\ and the tail
parameter \\\nu\_{\text{raw}}\\ show somewhat broader posteriors,
reflecting the inherent difficulty of estimating distributional shape
from a finite sample, but even these are far narrower than their priors.

## Posterior predictive checks

Before proceeding to formal diagnostics, a visual check of whether the
model can reproduce the observed data is invaluable. Posterior
predictive checks (PPCs) simulate new datasets from the posterior
predictive distribution \\p(\tilde{y} \mid y) = \int p(\tilde{y} \mid
\theta) p(\theta \mid y) d\theta\\ and compare them against the observed
data [\[11\]](#ref11). Systematic discrepancies between simulated and
observed data indicate model misspecification, even when the model
passes all convergence diagnostics.

The [`predict()`](https://rdrr.io/r/stats/predict.html) method generates
posterior predictive samples, returning a `demonoid.ppc` object whose
[`ppc_dens_overlay()`](https://robustecologies.github.io/lucifer/reference/ppc_dens_overlay.md)
method overlays the predicted distribution against the observations.

``` r

# Posterior predictive check for the Student-t model
# (ppc_t was generated above for the fitted-curve fan chart)
ppc_dens_overlay(ppc_t)
```

For a well-specified model, the observed data should appear as a typical
draw from the posterior predictive distribution. The PPC plot
superimposes the density of the observed \\y\\ against densities of
replicated datasets \\\tilde{y}\\. If the observed density falls outside
the range of the replicated densities, the model is failing to capture
some aspect of the data. For instance, if the observed data show heavier
tails than the replicated data, the error distribution is too thin (as
would happen with the normal model applied to Student-\\t\\ data).

``` r

# Compare PPC across models
ppc_normal <- predict(fit_normal, Model = spec_normal$Model, Data = Data_normal)
ppc_dens_overlay(ppc_normal)

ppc_laplace <- predict(fit_laplace, Model = spec_laplace$Model, Data = Data_laplace)
ppc_dens_overlay(ppc_laplace)
```

Comparing PPCs across models provides qualitative evidence that
complements the quantitative model comparison in later sections. The
normal model’s PPC will show replicated datasets that are too
concentrated near the mean, failing to reproduce the extreme
observations. The Student-\\t\\ model’s PPC should match the observed
distribution well. The Laplace model will show intermediate behavior.

## Diagnostics with Consort()

Convergence is not a property that can be proven; it can only be
diagnosed. The
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function implements modern convergence diagnostics following the
recommendations of Vehtari et al. [\[2\]](#ref2) and returns a
structured S3 object with pass/fail assessments, a continuous
convergence score, and actionable suggestions when problems are
detected.

``` r

cx <- Consort(fit_t)
print(cx)
```

The output reports split-\\\hat{R}\\ (threshold \\\< 1.01\\), bulk and
tail ESS (thresholds \\\geq 400\\ and \\\geq 200\\, respectively), MCSE
relative to the posterior standard deviation, divergent transitions, and
energy BFMI for gradient-based methods. A composite convergence score
between 0 and 100 summarizes all diagnostics. For a well-behaved
four-parameter regression with NUTS, the score should be near 100 with
all diagnostics passing.

When
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
detects problems, it provides specific suggestions. Low ESS triggers a
recommendation to increase `Iterations` or reduce `Thinning`; elevated
\\\hat{R}\\ suggests running more chains or extending the warmup;
divergent transitions prompt a step-size reduction or model
reparameterization. These suggestions are structured data, not just
console messages, which is what allows
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
to act on them automatically.

The distinction between bulk ESS and tail ESS deserves emphasis. Bulk
ESS measures how well the sampler explores the central region of the
posterior, which determines the reliability of posterior means and
medians. Tail ESS measures how well it explores the extremes, which
determines the reliability of credible intervals and tail probabilities.
A sampler can have high bulk ESS but low tail ESS if the chains mix well
near the mode but rarely visit the tails, a pattern that produces
apparently precise means but unreliable intervals. Consort reports both
and flags either if it falls below threshold.

``` r

plot(cx)
plot(cx, type = "scorecard")
```

## Model comparison with Arena()

Given multiple fitted models,
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
provides a unified comparison framework. It accepts a named list of fit
objects and computes DIC, WAIC, ESS per second, convergence scores, and
runtime for each. The comparison is meaningful even across different
inference engines: an MCMC fit can be compared against a variational fit
or a Laplace approximation on the same metrics.

``` r

arena <- Arena(x = list(
    normal  = fit_normal,
    student_t = fit_t,
    laplace = fit_laplace
))
print(arena)
summary(arena)
```

The printed output ranks models by their combined score, which weights
predictive accuracy (DIC, WAIC), computational efficiency (ESS/s), and
convergence quality. For data generated from a Student-\\t\\
distribution, the Student-\\t\\ model should rank first, with the
Laplace model second and the Gaussian model last. The DIC difference
between the Student-\\t\\ and Gaussian models quantifies the cost of
misspecifying the error distribution.

``` r

# Comparison dashboard
plot(arena)

# Forest plot of parameter estimates across models
plot(arena, type = "forest")
```

The forest plot displays point estimates and credible intervals from all
models for each shared parameter, revealing how error distribution
assumptions affect parameter inference. The regression coefficients
\\\beta_0\\ and \\\beta_1\\ should be similar across models (the mean
function is well identified by the data), but the scale parameter
\\\sigma\\ will differ because the Gaussian model must inflate
\\\sigma\\ to accommodate outliers that the Student-\\t\\ model handles
through its tail parameter.

This pattern, where location parameters are robust to error
specification but scale parameters are not, is a general phenomenon. It
arises because the mean function is determined primarily by the bulk of
the data, which all three error distributions model similarly, while the
scale parameter must account for the tails, where the distributions
differ fundamentally.

## Automated pipeline with Crucible()

For users who want the entire workflow automated,
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
wraps
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
into a single function call. It profiles the model, selects the most
promising methods, fits them with parallel chains, diagnoses
convergence, iteratively refines non-converged fits, and returns a
ranked comparison.

``` r

set.seed(666)
crucible <- Crucible(spec_t$Model, Data_t, IV_t,
                     n_methods = 5, Chains = 4, CPUs = 4)
print(crucible)
summary(crucible)
plot(crucible)
```

The `n_methods` argument controls how many inference algorithms are
tested.
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
selects the top-ranked methods from
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
across different families (at least one per family, where feasible),
fits each, and enters a refinement loop for any that fail convergence
diagnostics. The result is a comprehensive benchmark that would
otherwise require substantial manual effort.

The automation is valuable for two purposes. First, it provides a
defensible answer to the question “did you try other algorithms?” by
systematically testing the most promising alternatives. Second, it
surfaces unexpected performance differences: a variational method might
converge faster than NUTS for a particular posterior geometry, or a
slice sampler might outperform gradient-based methods when the
log-posterior has discontinuities.

The name evokes a severe test in which impurities are burned away:
models and algorithms enter the crucible, and only the most reliable
emerge. The metaphor extends the package’s naming convention, where
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
advises,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
diagnoses, and
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
compares, each role corresponding to a stage of the scientific workflow.

## Cross-validation

Model comparison via DIC or WAIC uses in-sample criteria that
approximate out-of-sample predictive performance. Cross-validation
provides a more direct estimate by evaluating predictions on held-out
data. The lucifer package implements three complementary approaches:
LOO-PSIS for efficient approximate leave-one-out cross-validation
[\[3\]](#ref3), exact K-fold cross-validation, and Bayesian model
stacking for combining models [\[4\]](#ref4).

### LOO-PSIS

Pareto smoothed importance sampling (PSIS) approximates the
leave-one-out predictive density from a single model fit by reweighting
posterior draws [\[3\]](#ref3).
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
requires an \\N \times S\\ matrix of pointwise log-likelihoods (N
observations, S posterior samples). For models built with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
the pointwise log-likelihoods are not monitored automatically; the user
must either write a manual Model function that returns per-observation
log-likelihoods in the `Monitor` vector, or compute them post-hoc. When
pointwise log-likelihoods are unavailable, K-fold cross-validation
(below) provides a robust alternative that works with any Model
function.

### K-fold cross-validation

K-fold cross-validation provides a direct, assumption-free estimate of
out-of-sample predictive performance. Each fold refits the model to a
training subset and evaluates predictions on the held-out observations.
It works with any Model function without requiring pointwise
log-likelihoods.

``` r

set.seed(22122)

# Extract a standalone Model function for serialization safety.
# model_spec closures carry their compilation environment, which
# may not survive PSOCK serialization. Resetting the environment
# to the package namespace avoids "object not found" errors on workers.
Model_t <- spec_t$Model
environment(Model_t) <- asNamespace("lucifer")

kf <- Kfold(Model_t, Data_t, fit = fit_t, K = 10,
            fit_args = list(Algorithm = "NUTS", Iterations = 10000),
            CPUs = 4)
print(kf)
```

The `fit` argument provides an initial fit to warm-start the
fold-specific fits, reducing the number of iterations needed for
convergence. The `CPUs` argument parallelizes across folds, which is
particularly important since K-fold requires \\K\\ separate model fits.

The computational cost of K-fold is approximately \\K\\ times that of a
single fit, which is why LOO-PSIS is preferred when it is reliable.
However, K-fold has two advantages: it makes no importance sampling
assumptions (so there are no \\\hat{k}\\ diagnostics to worry about),
and it naturally handles models where leave-one-out deletion
fundamentally changes the posterior geometry (e.g., when removing a
single observation from a small cluster eliminates that cluster). In
practice, the workflow is: try LOO-PSIS first; if any \\\hat{k} \>
0.7\\, fall back to K-fold for those observations or for the entire
comparison.

## Sensitivity analysis with RobustBayes()

A posterior that changes substantially under modest perturbations to the
prior, likelihood, or data is fragile, and conclusions drawn from it may
not be trustworthy.
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
provides a systematic framework for diagnosing three sources of
sensitivity [\[5\]](#ref5).

``` r

rb <- RobustBayes(fit_t, Model = spec_t$Model, Data = Data_t,
                  modules = c("power", "conflict", "influence"))
print(rb)
summary(rb)
```

The **power-scaling** module perturbs the prior and likelihood by
raising them to powers \\\alpha \in \[0.5, 2.0\]\\. If the posterior
mean of a parameter shifts substantially when the prior is strengthened
(\\\alpha \> 1\\) or weakened (\\\alpha \< 1\\), that parameter is
prior-sensitive, meaning the data are not sufficiently informative to
override the prior. Conversely, if the posterior is insensitive to prior
scaling, the inference is data-dominated. For a well-identified
regression with \\N = 300\\, all parameters should show low prior
sensitivity.

The **prior-data conflict** module detects whether the prior and the
likelihood pull the posterior in different directions. When the prior
mode and the data-preferred region of parameter space are far apart, the
posterior occupies an uncomfortable compromise position that neither the
prior nor the data support strongly. This module flags conflicts by
comparing prior predictive distributions against observed summary
statistics.

The **observation influence** module identifies observations that exert
disproportionate leverage on the posterior. It uses PSIS leave-one-out
weights to approximate the change in posterior summaries when each
observation is removed. Observations with high Pareto \\\hat{k}\\ values
are flagged as influential. In the presence of heavy-tailed data, the
Gaussian model will show more influential observations than the
Student-\\t\\ model, because the Gaussian cannot accommodate extreme
values through its distributional shape and must instead shift its
parameters.

``` r

plot(rb)
plot(rb, type = "influence")
```

The influence plot highlights individual observations and their effect
on parameter estimates. Comparing influence patterns across models
reveals whether a robust error specification (Student-\\t\\, Laplace)
successfully reduces the leverage of outlying observations.

A particularly instructive exercise is to run
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
on the Gaussian model and compare the results against the Student-\\t\\
model.

``` r

rb_normal <- RobustBayes(fit_normal, Model = spec_normal$Model, Data = Data_normal,
                         modules = c("power", "influence"))
plot(rb_normal, type = "influence")
```

The Gaussian model will show substantially more influential observations
than the Student-\\t\\ model, because extreme values that the \\t\\
distribution handles gracefully through its tail parameter must instead
be accommodated by shifting the Gaussian’s location and scale
parameters. This comparison provides both a diagnostic (the Gaussian is
fragile) and a motivation (the Student-\\t\\ is robust) for choosing
heavy-tailed error distributions.

The power-scaling analysis will also differ: the Gaussian model’s scale
parameter \\\sigma\\ may show elevated prior sensitivity because the
prior must compensate for the model’s inability to handle outliers
through the likelihood alone. In contrast, the Student-\\t\\ model
distributes the outlier accommodation between \\\sigma\\ (the scale) and
\\\nu\\ (the tail heaviness), reducing the burden on any single
parameter and resulting in lower sensitivity to prior perturbation.

### Direct prior sensitivity comparison

The
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
power-scaling approach perturbs prior strength continuously, but
sometimes the question is more concrete: how much would the conclusions
change under a genuinely different prior specification? The
[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md)
function addresses this by overlaying posterior densities obtained under
alternative priors, providing an immediate visual answer.

To illustrate, the Student-\\t\\ model is refitted with substantially
tighter priors: \\\text{Normal}(0, 1)\\ for the regression coefficients
(instead of \\\text{Normal}(0, 10)\\) and \\\text{HalfCauchy}(1)\\ for
the scale (instead of \\\text{HalfCauchy}(2.5)\\). If the posteriors are
nearly identical, the data overwhelm both priors and the analysis is
robust; if they diverge, at least one parameter is prior-sensitive and
the prior choice matters.

``` r

spec_t_tight <- model_spec("
    y ~ StudentT(mu, sigma, nu)
    mu = beta0 + beta1 * x
    beta0 ~ Normal(0, 1)
    beta1 ~ Normal(0, 1)
    sigma ~ HalfCauchy(1)
    nu_raw ~ Exponential(0.1)
    nu = nu_raw + 2
")

Data_t_tight <- spec_t_tight$data_template(y = y, x = x)
IV_t_tight <- spec_t_tight$initial_values(Data_t_tight)

set.seed(1232)
fit_t_tight <- lucifer(spec_t_tight$Model, Data_t_tight, IV_t_tight,
                       Iterations = 10000, Status = 5000, Thinning = 10,
                       Algorithm = "NUTS", Chains = 2, CPUs = 2)
```

``` r

prior_sensitivity(
    fits = list("Diffuse: N(0,10), HC(2.5)" = fit_t,
                "Tight: N(0,1), HC(1)"      = fit_t_tight),
    ground_truth = c(beta0 = beta0_true, beta1 = beta1_true,
                     sigma = sigma_true, nu_raw = nu_true - 2)
)
```

With \\N = 300\\ observations, the posteriors under both prior
specifications should be virtually indistinguishable for the regression
coefficients, confirming that the data are sufficiently informative to
render the prior irrelevant for location parameters. The scale and tail
parameters may show minor differences at the tails of the posterior,
reflecting the tighter prior’s gentle pull, but the central summaries
(means, medians) should agree closely. This is the hallmark of a
data-dominated analysis where the prior choice is defensible regardless
of which specification is adopted.

## Identifiability via DataCloning()

Identifiability, the question of whether the data can uniquely determine
the model parameters, is a prerequisite for meaningful inference that is
surprisingly often violated in complex models. Data cloning
[\[6\]](#ref6) [\[7\]](#ref7) provides a computationally elegant
diagnostic: the observed data are conceptually replicated \\K\\ times,
and as \\K\\ increases the posterior concentrates on the maximum
likelihood estimate at rate \\1/K\\. For identifiable parameters, the
posterior variance shrinks proportionally to \\1/K\\; for
non-identifiable parameters, the variance stabilizes or shrinks at a
different rate.

The theoretical basis rests on a Bernstein-von Mises argument. The
cloned posterior is \\\pi_K(\theta \mid y) \propto \[L(\theta \mid
y)\]^K \pi(\theta)\\, which converges to \\N(\hat{\theta}\_{\text{MLE}},
K^{-1} I(\hat{\theta})^{-1})\\ as \\K \to \infty\\ for identifiable
models. The eigenvalues of \\K \cdot \text{Var}\_K(\theta)\\ converge to
the eigenvalues of the inverse Fisher information
\\I(\hat{\theta})^{-1}\\. If any eigenvalue fails to stabilize, the
corresponding eigenvector direction is non-identifiable.

``` r

set.seed(3120)
dc <- DataCloning(spec_t$Model, Data_t, IV_t,
                  Iterations = 10000, Status = 2000, Thinning = 10,
                  Algorithm = "NUTS",
                  n_clones = c(1, 2, 5, 10),
                  Chains = 4, CPUs = 4)
print(dc)
summary(dc)
```

The [`summary()`](https://rdrr.io/r/base/summary.html) output reports
the posterior mean and scaled variance (\\K \times \text{Var}\\) at each
clone level. For identifiable parameters, the scaled variance should
stabilize across clone levels. The MLE estimates from the highest clone
level provide a useful frequentist benchmark against which to compare
the Bayesian posterior means from the \\K = 1\\ (uncloned) fit.

``` r

# Variance decay curves
plot(dc, type = 'diagnostics')
```

The eigenvalue plot is the most informative diagnostic. For a
four-parameter regression, all four eigenvalues should converge to
finite, positive limits as \\K\\ increases. An eigenvalue that remains
large or grows indicates a direction in parameter space along which the
likelihood is flat, meaning some parameter or combination of parameters
is not identifiable from the data.

Data cloning also serves as a bridge between Bayesian and frequentist
inference. At \\K = 1\\, the posterior is the standard Bayesian
posterior incorporating the prior. As \\K\\ increases, the prior’s
influence diminishes and the posterior concentrates on the MLE. The
scaled posterior covariance \\K \cdot \text{Var}\_K(\theta)\\ converges
to the inverse Fisher information, providing automatic standard error
estimates without computing second derivatives. This connection is
developed fully in the data cloning vignette
([`vignette("data_cloning")`](https://robustecologies.github.io/lucifer/articles/data_cloning.md))
and the frequentist bridge vignette
([`vignette("frequentist_bridge")`](https://robustecologies.github.io/lucifer/articles/frequentist_bridge.md)).

For models where identifiability is questionable, such as mixture models
with unknown numbers of components, overparameterized hierarchical
structures, or models with ridge-like likelihoods, data cloning provides
definitive evidence. An eigenvalue that fails to stabilize identifies a
non-estimable direction, and the corresponding eigenvector reveals which
parameter combinations are confounded. This information is difficult to
obtain by other means and can guide model reparameterization or
simplification.

## Hierarchical extension

The workflow extends naturally to hierarchical models, where parameters
are themselves drawn from a population-level distribution. This
structure enables partial pooling: group-level estimates borrow strength
from the population, shrinking extreme estimates toward the grand mean.
The effect is strongest for groups with small sample sizes, which is
precisely where raw estimates are least reliable.

### Simulating hierarchical data

The data-generating process creates \\J = 8\\ groups with \\N_j = 40\\
observations each. Each group has its own intercept \\\alpha_j\\ drawn
from a normal population distribution, while the slope and error
parameters are shared across groups.

``` r

set.seed(330)
J <- 8
N_per <- 40
N_total <- J * N_per
group <- rep(1:J, each = N_per)
x_h <- rnorm(N_total)

# Population-level parameters
mu_alpha_true <- 2.0
sigma_alpha_true <- 1.5
alpha_true <- rnorm(J, mu_alpha_true, sigma_alpha_true)
beta_true <- 1.5
sigma_h_true <- 1.0
nu_h_true <- 4

y_h <- alpha_true[group] + beta_true * x_h +
    sigma_h_true * rt(N_total, df = nu_h_true)
```

### Specifying and fitting the hierarchical model

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
syntax supports indexed parameters with the `j = 1,...,J` notation,
which declares a vector of group-level parameters with a shared prior.

``` r

spec_hier <- model_spec("
    y ~ StudentT(mu, sigma, nu)
    mu = alpha[group] + beta * x
    alpha[j] ~ Normal(mu_alpha, sigma_alpha), j = 1,...,J
    mu_alpha ~ Normal(0, 10)
    sigma_alpha ~ HalfCauchy(2.5)
    beta ~ Normal(0, 10)
    sigma ~ HalfCauchy(2.5)
    nu_raw ~ Exponential(0.1)
    nu = nu_raw + 2
")

Data_hier <- spec_hier$data_template(y = y_h, x = x_h, group = group, J = J)
IV_hier <- spec_hier$initial_values(Data_hier)
cspec_hier <- compile(spec_hier, Data_hier)

set.seed(1010011010)
fit_hier <- lucifer(spec_hier$Model, Data_hier, IV_hier,
                    Iterations = 30000, Status = 5000, Thinning = 10,
                    Algorithm = "NUTS", Chains = 4, CPUs = 4)
```

### Shrinkage

The defining feature of hierarchical models is shrinkage: group-level
estimates are pulled toward the population mean, with the degree of
shrinkage inversely proportional to the group sample size and directly
proportional to the population variance. Comparing raw group means (no
pooling) against posterior means (partial pooling) against the
population mean (complete pooling) reveals this phenomenon.

``` r

# Raw group means (no pooling)
raw_means <- tapply(Data_hier$y, Data_hier$group, mean)

# Posterior group intercepts (partial pooling)
hier_means <- setNames(fit_hier$Summary2[, "Mean"],
                       rownames(fit_hier$Summary2))
alpha_names <- paste0("alpha[", 1:J, "]")
post_alpha <- hier_means[alpha_names]

# Population mean (complete pooling)
mu_alpha_post <- hier_means["mu_alpha"]

shrinkage_df <- data.frame(
    group = 1:J,
    true = alpha_true,
    raw = raw_means,
    posterior = post_alpha,
    population = rep(mu_alpha_post, J)
)

knitr::kable(shrinkage_df, digits = 3,
             caption = "Shrinkage: raw group means vs. partial pooling estimates")
```

``` r

ggplot(shrinkage_df, aes(x = raw, y = posterior)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = mu_alpha_post, linetype = "dotted", color = "#D6604D") +
    geom_point(size = 3, color = "#2166AC") +
    geom_text(aes(label = group), nudge_y = 0.15, size = 3) +
    labs(x = "Raw group mean (no pooling)",
         y = "Posterior mean (partial pooling)",
         title = "Hierarchical shrinkage toward the population mean") +
    theme_minimal()
```

Groups whose raw means are extreme (far from the population mean) show
the most shrinkage, while groups near the center are barely affected.
The true group intercepts (known from simulation) will generally be
closer to the partial pooling estimates than to the raw means,
demonstrating that shrinkage improves estimation accuracy, not just
precision.

The amount of shrinkage is governed by the ratio of within-group
variance to between-group variance. When within-group variance is large
relative to between-group variance (as when \\\sigma / \sigma\_\alpha\\
is large), the raw group means are noisy and shrinkage is substantial.
When the ratio is small, the raw means are precise and shrinkage is
minimal. The hierarchical model learns this ratio from the data through
the hyperparameters \\\mu\_\alpha\\ and \\\sigma\_\alpha\\, adaptively
calibrating the degree of pooling. This is sometimes called “learning
the prior from the data,” though more precisely it is estimating the
population-level parameters that define the prior for the group-level
parameters.

The posterior estimate of \\\sigma\_\alpha\\ itself is informative: it
quantifies how much genuine variation exists between groups after
accounting for within-group noise. If \\\sigma\_\alpha\\ is estimated
near zero, the model is effectively pooling all groups toward a common
intercept. If \\\sigma\_\alpha\\ is large, the groups are genuinely
different and the model preserves their individuality. Comparing the
posterior of \\\sigma\_\alpha\\ against its true value (1.5 in this
simulation) provides another validation of the model’s accuracy.

``` r

summary(fit_hier)
plot(fit_hier)

cx_hier <- Consort(fit_hier)
print(cx_hier)
```

Hierarchical models with many group-level parameters can present
challenges for MCMC convergence, particularly when the group-level
variance \\\sigma\_\alpha\\ is small (near the boundary of its parameter
space) or when the number of groups is small. The “funnel” geometry that
arises when group-level parameters are tightly constrained by a small
\\\sigma\_\alpha\\ but loosely constrained by the data is a well-known
difficulty for Hamiltonian Monte Carlo [\[12\]](#ref12). Non-centered
parameterizations, where \\\alpha_j = \mu\_\alpha + \sigma\_\alpha \cdot
z_j\\ with \\z_j \sim N(0, 1)\\, can alleviate this problem by
separating the scale parameter from the individual effects. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL does not currently perform this reparameterization automatically, so
for challenging hierarchical models the manual approach may be
necessary.

``` r

# Non-centered parameterization (manual)
spec_hier_ncp <- model_spec("
    y ~ StudentT(mu, sigma, nu)
    mu = mu_alpha + sigma_alpha * z[group] + beta * x
    z[j] ~ Normal(0, 1), j = 1,...,J
    mu_alpha ~ Normal(0, 10)
    sigma_alpha ~ HalfCauchy(2.5)
    beta ~ Normal(0, 10)
    sigma ~ HalfCauchy(2.5)
    nu_raw ~ Exponential(0.1)
    nu = nu_raw + 2
")
```

The non-centered parameterization places the standard normal prior on
`z[j]` and constructs the group intercept as a deterministic function.
This eliminates the funnel geometry because the sampler explores the
\\z_j\\ parameters, which have unit-scale priors regardless of the value
of \\\sigma\_\alpha\\. The group intercepts \\\alpha_j = \mu\_\alpha +
\sigma\_\alpha z_j\\ are recovered as derived quantities. For problems
where the centered parameterization produces divergent transitions or
low ESS for \\\sigma\_\alpha\\, switching to the non-centered form often
resolves the issue.

## Frequentist bridge

The Bayesian posterior and the frequentist sampling distribution
converge under specific conditions: when the prior is flat relative to
the likelihood and the sample size is large enough for the Bernstein-von
Mises theorem to apply [\[8\]](#ref8). The lucifer package provides
explicit tools for computing frequentist quantities from Bayesian fits,
enabling direct comparison between the two paradigms.

### Frequentist summary

The
[`freq_summary()`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)
function extracts maximum likelihood estimates, standard errors (from
the inverse observed Fisher information), Wald confidence intervals, and
\\p\\-values from a fitted object, treating the posterior mode as the
MLE and the posterior curvature as the observed information.

``` r

fs <- freq_summary(fit_t)
print(fs)
plot(fs)
```

### Comparing confidence and credible intervals

The
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
function places Bayesian credible intervals alongside frequentist
confidence intervals and (when available) data cloning intervals, with
the true parameter values marked for reference.

``` r

# true_values must be a named vector whose names match the parameter
# names in the fit objects (beta0, beta1, sigma, nu_raw)
ci <- confint_compare(
    fits = list(Bayesian = fit_t, DataCloning = dc),
    true_values = c(beta0 = beta0_true, beta1 = beta1_true,
                    sigma = sigma_true, nu_raw = nu_true - 2)
)
print(ci)
```

For a well-identified model with weakly informative priors and \\N =
300\\ observations, the Bayesian and frequentist intervals will be
nearly identical, confirming that the prior has negligible influence.
The data cloning intervals, which are asymptotically exact frequentist
intervals, provide an additional reference point.

### Profile likelihood

The
[`profile_likelihood()`](https://robustecologies.github.io/lucifer/reference/profile_likelihood.md)
function computes the profile log-likelihood for each parameter by
maximizing the log-likelihood over all other parameters at a grid of
fixed values. The resulting curves reveal the shape of the likelihood
surface and produce likelihood-based confidence intervals that are valid
even when the normal approximation underlying Wald intervals fails.

``` r

pl <- profile_likelihood(fit_t, Model = spec_t$Model, Data = Data_t)
plot(pl)
```

The profile likelihood plot overlays the profile curve, the quadratic
Wald approximation, and the Bayesian marginal posterior for each
parameter. When all three agree, the inference is robust to the choice
of framework. When they disagree, particularly for parameters near
boundaries or in models with non-quadratic likelihoods, the profile
likelihood provides more reliable frequentist intervals than the Wald
approximation.

The frequentist bridge is not merely a curiosity. It serves three
practical purposes. First, for users transitioning from frequentist to
Bayesian methods, it provides a concrete verification that the Bayesian
machinery reproduces familiar results when priors are diffuse. This
builds confidence in the tools before introducing informative priors.
Second, for reviewers and collaborators who prefer frequentist
reporting, it allows the same model fit to produce both Bayesian
credible intervals and frequentist confidence intervals without
refitting. Third, it clarifies exactly where and why the two frameworks
diverge: when priors are informative, when the posterior is skewed or
multimodal, or when the sample size is too small for the Bernstein-von
Mises theorem to apply [\[8\]](#ref8). Understanding these boundaries is
essential for making a principled choice between paradigms.

## Summary and next steps

This vignette has traced a complete Bayesian workflow from data
simulation through model specification, algorithm selection, fitting,
convergence diagnostics, model comparison, cross-validation, sensitivity
analysis, identifiability testing, hierarchical extension, and
frequentist bridging. The core pattern is consistent across all of these
steps: write a Model function (or use
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
to generate one), prepare a Data list, and pass both to the appropriate
lucifer function.

The workflow is not strictly linear. In practice, an analyst iterates
between stages: a failed convergence diagnostic sends her back to
algorithm selection; a model comparison reveals that none of the
candidates fit well, prompting model revision; a sensitivity analysis
shows that the posterior is fragile, motivating a more careful prior
elicitation. The sections in this vignette are presented sequentially
for clarity, but the actual workflow is cyclical, and the lucifer
ecosystem is designed to support rapid iteration.

The key stages, in the order they typically appear in a real analysis,
are as follows.

``` r

# 1. Specify and prepare data
spec <- model_spec("...")
Data <- spec$data_template(...)
IV <- spec$initial_values(Data)

# 2. Profile and select algorithm
rx <- Prescribe(spec$Model, Data, IV)

# 3. Fit (multi-chain)
fit <- lucifer(spec$Model, Data, IV, Algorithm = "NUTS",
               Iterations = 20000, Chains = 4, CPUs = 4)

# 4. Diagnose convergence
cx <- Consort(fit)

# 5. Compare models
arena <- Arena(x = list(m1 = fit1, m2 = fit2, m3 = fit3))

# 6. Cross-validate
loo_result <- LOO(fit)

# 7. Assess sensitivity
rb <- RobustBayes(fit, Model = spec$Model, Data = Data)

# 8. Test identifiability
dc <- DataCloning(spec$Model, Data, IV, n_clones = c(1, 2, 5, 10),
                  Algorithm = "NUTS", Chains = 4, CPUs = 4)
```

Several other vignettes extend this foundation to specialized domains.

``` r

vignettes <- data.frame(
    Vignette = c(
        "model_spec", "mcmc-algorithms", "cross_validation",
        "robust_bayes", "data_cloning", "frequentist_bridge",
        "crucible", "prescribe_arena", "sde", "sbi", "node",
        "variational_bayes", "interoperability"
    ),
    Topic = c(
        "Declarative model specification DSL, indexed parameters, truncated distributions",
        "Survey of 82 MCMC algorithms with guidance on when to use each",
        "LOO-PSIS, K-fold, leave-future-out, stacking weights",
        "Prior sensitivity, prior-data conflict, observation influence",
        "Maximum likelihood via data cloning, identifiability testing",
        "Frequentist summaries, Wald/LR/score tests, profile likelihood",
        "Automated pipeline: Prescribe, fit, Consort, Arena in one call",
        "Algorithm selection profiling and multi-method benchmarking",
        "Stochastic differential equations, particle filters, exact likelihood",
        "Simulation-based inference: NPE, NLE, NRE, TSNPE",
        "Neural ODEs via Bayesian gradient matching",
        "ADVI, Pathfinder, normalizing flows",
        "Stan, JAGS, brms import/export"
    )
)
knitr::kable(vignettes, caption = "Companion vignettes for specialized topics")
```

The model specification vignette covers the full
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
grammar including indexed parameters, truncated distributions, and LaTeX
notation input. The MCMC algorithms vignette surveys all 82 MCMC methods
with guidance on when to use each. The SDE vignette demonstrates
stochastic differential equation models with particle filtering and
exact likelihood computation. The SBI vignette covers simulation-based
inference for models where the likelihood is intractable. The NODE
vignette introduces neural ordinary differential equations via Bayesian
gradient matching. The cross-validation vignette provides deeper
coverage of LOO-PSIS, K-fold, leave-future-out, and stacking. And the
robust Bayes vignette develops the full theory and practice of prior
sensitivity, likelihood sensitivity, and observation influence analysis.

Two general principles run through the entire workflow. First, always
simulate data with known ground truth before analyzing real data. If the
model cannot recover known parameters from simulated data, it will not
produce reliable inferences from real data, where the ground truth is
unknown. Second, never trust a single diagnostic or a single model. The
strength of the lucifer ecosystem lies in the multiplicity of
perspectives it provides: convergence via Consort, comparison via Arena,
prediction via LOO, sensitivity via RobustBayes, identifiability via
DataCloning. A result that survives all of these challenges is far more
credible than one that has been subjected to only a cursory check.

The complete workflow demonstrated in this vignette, from simulation
through model specification, fitting, diagnosis, comparison,
sensitivity, and identifiability, represents current best practice in
applied Bayesian statistics. Each component addresses a distinct source
of potential failure: convergence diagnostics guard against MCMC
artifacts, model comparison guards against structural misspecification,
sensitivity analysis guards against prior influence, and identifiability
testing guards against fundamental estimation impossibility. No single
tool catches all problems, but together they provide a comprehensive
safety net that substantially reduces the risk of publishing inferences
that are driven by computational artifacts or modeling assumptions
rather than by the data.

## References

**\[1\]** Gelman, A. (2006). *Prior distributions for variance
parameters in hierarchical models*. Bayesian Analysis, 1(3), 515-534.
[doi:10.1214/06-BA117A](https://doi.org/10.1214/06-BA117A)

**\[2\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., &
Burkner, P.-C. (2021). *Rank-normalization, folding, and localization:
an improved \\\hat{R}\\ for assessing convergence of MCMC (with
discussion)*. Bayesian Analysis, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

**\[3\]** Vehtari, A., Gelman, A., & Gabry, J. (2017). *Practical
Bayesian model evaluation using leave-one-out cross-validation and
WAIC*. Statistics and Computing, 27(5), 1413-1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[4\]** Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). *Using
stacking to average Bayesian predictive distributions (with
discussion)*. Bayesian Analysis, 13(3), 917-1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

**\[5\]** Kallioinen, N., Paananen, T., Burkner, P.-C., & Vehtari, A.
(2024). *Detecting and diagnosing prior and likelihood sensitivity with
power-scaling*. Statistics and Computing, 34, 57.
[doi:10.1007/s11222-023-10366-5](https://doi.org/10.1007/s11222-023-10366-5)

**\[6\]** Lele, S. R., Dennis, B., & Lutscher, F. (2007). *Data cloning:
easy maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods*. Ecology Letters, 10(7),
551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

**\[7\]** Lele, S. R., Nadeem, K., & Schmuland, B. (2010). *Estimability
and likelihood inference for generalized linear mixed models using data
cloning*. Journal of the American Statistical Association, 105(492),
1617-1625.
[doi:10.1198/jasa.2010.tm09757](https://doi.org/10.1198/jasa.2010.tm09757)

**\[8\]** van der Vaart, A. W. (1998). *Asymptotic Statistics*.
Cambridge University Press. ISBN: 978-0-521-78450-4.

**\[9\]** Hoffman, M. D. & Gelman, A. (2014). *The No-U-Turn Sampler:
adaptively setting path lengths in Hamiltonian Monte Carlo*. Journal of
Machine Learning Research, 15(47), 1593-1623.
<https://jmlr.org/papers/v15/hoffman14a.html>

**\[10\]** Girolami, M. & Calderhead, B. (2011). *Riemann manifold
Langevin and Hamiltonian Monte Carlo methods*. Journal of the Royal
Statistical Society: Series B (Statistical Methodology), 73(2), 123-214.
[doi:10.1111/j.1467-9868.2010.00765.x](https://doi.org/10.1111/j.1467-9868.2010.00765.x)

**\[11\]** Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.).
Chapman and Hall/CRC. ISBN: 978-1-4398-4095-5.

**\[12\]** Betancourt, M. (2017). *A conceptual introduction to
Hamiltonian Monte Carlo*. arXiv:1701.02434.
[doi:10.48550/arXiv.1701.02434](https://doi.org/10.48550/arXiv.1701.02434)
