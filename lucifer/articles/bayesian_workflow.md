# Bayesian workflow with lucifer

## Introduction

### From inference to workflow

Bayesian inference, the computation of \\p(\theta \mid y) \propto p(y
\mid \theta)\\ p(\theta)\\, is the formal core of Bayesian statistics,
but it is only one component of applied Bayesian data analysis. A
practitioner who runs a sampler and reports the posterior has completed
inference; a practitioner who also checked the prior implications,
validated the computation against simulated data, diagnosed convergence,
tested the model against the data it was supposed to explain, assessed
sensitivity to modeling choices, compared alternative specifications,
and propagated uncertainty into predictions has completed a workflow.
The distinction matters because each of these steps can reveal problems
that invalidate the inference, and the most informative failures are
typically discovered not during fitting but during the checks that
surround it.

Gelman et al. [\[1\]](#ref1) formalize this perspective in a
comprehensive treatment that draws on decades of applied experience.
Their central argument is that the tangled, nonlinear nature of real
data analysis should be acknowledged and systematized rather than hidden
behind a single call to a fitting function. The workflow they describe
is inherently iterative: the analyst specifies a model, checks its prior
implications, fits it, diagnoses the computation, evaluates the fit
against the data, and, finding the model inadequate, modifies it and
repeats the cycle. The final model is typically the survivor of several
such iterations, and the abandoned models along the way are not wasted
effort but essential scaffolding that builds understanding.

This iterative perspective has deep roots. Box [\[7\]](#ref7) introduced
the idea of a feedback loop between model criticism and model revision
in his 1980 paper on sampling and Bayes’ theorem, arguing that “all
models are wrong, but some are useful” and that the useful ones are
identified through systematic checking. Rubin [\[8\]](#ref8) developed
posterior predictive checking as a formal tool for this purpose, showing
how simulated replications from the fitted model can diagnose specific
forms of misfit. Berger [\[9\]](#ref9) emphasized the role of robustness
analysis in Bayesian statistics, arguing that conclusions should be
stable under reasonable perturbations of the prior and likelihood. More
recently, Betancourt [\[10\]](#ref10) has articulated a principled
workflow for probabilistic computation that integrates domain expertise,
prior calibration, and computational diagnostics into a coherent
methodology. Gabry et al. [\[6\]](#ref6) demonstrated how visualization
can serve as the connective tissue between these steps, making the
workflow tangible and actionable.

### Why not just fit the model?

The frequentist tradition, built on the Neyman-Pearson framework of
hypothesis testing and confidence intervals, separates estimation from
model checking in a fundamental way. A maximum likelihood estimate is
computed, its sampling distribution is characterized (analytically or
via the bootstrap), and goodness-of-fit tests are applied as a separate,
post-hoc diagnostic. The model is either accepted or rejected; there is
no formal mechanism for iterative revision within the framework itself.

The Bayesian approach differs in two respects that make a workflow both
possible and necessary. First, the generative model, which specifies a
joint distribution over parameters and data, \\p(\theta, y) =
p(\theta)\\ p(y \mid \theta)\\, enables simulation from the model at
every stage: before seeing data (prior predictive), after fitting
(posterior predictive), and under perturbations (sensitivity analysis).
These simulations provide a continuous diagnostic signal rather than a
binary accept/reject verdict. Second, the prior distribution is an
explicit modeling choice that must be calibrated against domain
knowledge, not a nuisance to be eliminated by asymptotics. This
calibration step, which has no analogue in frequentist practice, creates
a natural entry point for iterative refinement.

The cost of ignoring the workflow is well documented. Prior predictive
checks can reveal that apparently “non-informative” priors imply absurd
predictions (see Section 2.4 of [\[1\]](#ref1) and Simpson et
al. [\[11\]](#ref11) on penalized complexity priors). Convergence
diagnostics can detect computational failures that produce
plausible-looking but incorrect posteriors (see Section 5.1 of
[\[1\]](#ref1) on the “folk theorem of statistical computing”).
Posterior predictive checks can reveal model misspecification that
inflates or deflates uncertainty in ways that information criteria alone
cannot detect. Sensitivity analysis can reveal that conclusions depend
more on the prior than on the data, a situation that is invisible
without explicit perturbation. Each of these failures is invisible to an
analyst who only runs the sampler.

### The workflow as a graph

Textbook presentations of statistical workflow are often linear: define
a model, fit it, report the results. Real analyses are not linear. The
analyst moves forward and backward through multiple stages, sometimes
returning to the beginning when a fundamental assumption proves wrong,
sometimes making small adjustments within a single stage. The challenge
for software is to support this nonlinearity without imposing either too
much structure (which constrains the analyst) or too little (which
provides no guidance).

lucifer addresses this challenge through a layered architecture where
each workflow step is implemented by a specific function but no step
requires the others. The analyst can use
[`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
without ever fitting a model, or fit a model with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
without running
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
afterward. The functions are designed to compose but not to constrain.
The workflow diagram below makes this composability explicit by
organizing the functions into 14 nodes connected by directed arrows, but
the analyst is free to enter the diagram at any node and follow any
path.

The layout distinguishes two streams. The **analytical pipeline** (blue
and green nodes) represents the forward flow from scientific theory
through model specification, data collection, estimation, and
prediction. The **validation layer** (amber nodes) represents the
diagnostic and synthetic checks that provide feedback loops, potentially
sending the analyst back to modify the model or the computation. Every
node and every arrow corresponds to a specific lucifer function, and
this vignette demonstrates each one through a running example of
increasing complexity.

The diagram encodes three types of information flow through its arrow
styles. **Solid arrows** represent the deterministic forward flow from
theory through estimation to prediction: these steps execute once per
model specification. **Dashed arrows** represent the diagnostic pathways
where synthetic or real data are checked against model expectations.
**Dotted bidirectional arrows** represent the iterative feedback loops
where a diagnostic failure sends the analyst back to modify the
estimator. The bold dashed red loop from Estimates through Diagnostics
back to Estimator is the core refinement cycle that drives model
improvement.

### lucifer’s coverage of the workflow

lucifer implements every step of this workflow through a unified
interface. The package provides over 130 inference algorithms across
eight families, but inference is only one node in the workflow graph.
Prior predictive simulation
([`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)),
algorithm recommendation
([`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)),
convergence diagnostics
([`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)),
posterior predictive checks (40 styles via
[`plot.demonoid.ppc()`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md)),
cross-validation
([`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)),
sensitivity analysis
([`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)),
identifiability and estimability assessment
([`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)),
model comparison
([`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md),
[`BayesFactor()`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md)),
and automated orchestration
([`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md))
together cover the full workflow without leaving the package ecosystem.
The unifying contract across all of these, a single Model function
returning `LP`, `Dev`, `Monitor`, `yhat`, and `parm`, means that every
diagnostic, comparison, and visualization tool works identically
regardless of which inference engine produced the fit.

### Roadmap

The vignette follows a dose-response analysis through four models of
increasing complexity. At each stage, a diagnostic failure motivates
expansion to the next model, demonstrating the iterative character of
the workflow. The mapping between diagram nodes and lucifer functions is
summarized in the table at the end. The reader who works through this
vignette will have exercised every node and arrow in the diagram.

## Theory, estimand, and design

### Scientific context

The scientific problem is the analysis of a nonlinear dose-response
relationship across multiple experimental groups. Organisms in \\G = 4\\
populations are exposed to a continuous dose \\d\\ and a continuous
response \\y\\ is measured. The response saturates at high doses,
following Michaelis-Menten kinetics, a functional form that arises in
enzyme kinetics, receptor binding, and ecological resource uptake:

\\ \mu(d) = \frac{V\_{\max} \cdot d}{K_m + d} \\

where \\V\_{\max}\\ is the maximum response (the asymptotic plateau) and
\\K_m\\ is the dose at which the response reaches half its maximum. Both
parameters are strictly positive and have direct biological
interpretations, making them natural estimands. The residual variation
\\\sigma\\ captures measurement noise and unmodeled biological
variability.

The Michaelis-Menten function has two important properties that make it
a good vehicle for demonstrating the workflow. First, it is nonlinear in
the parameters, so that standard linear-model diagnostics are
insufficient and the full Bayesian machinery (prior predictive checking,
posterior predictive checking, sensitivity analysis) becomes essential
rather than optional. Second, the saturation behavior encodes genuine
scientific knowledge: the response cannot grow without bound, and any
model that predicts unbounded growth is incompatible with the theory,
regardless of what the data say. This tension between data and theory is
exactly the kind of constraint that the prior predictive check in
Section 4 is designed to detect.

### Defining the estimand

The estimand, the quantity the analysis aims to estimate, must be
specified before any model is written. This discipline, emphasized by
Rubin [\[8\]](#ref8) and operationalized in the potential outcomes
framework, forces the analyst to distinguish the scientific question
from the statistical model. Two analyses can use the same model but
target different estimands, and two estimands can require different
models even when applied to the same data. Separating the estimand from
the estimator is important because it prevents the common mistake of
letting computational convenience drive scientific conclusions: if the
easiest model to fit does not identify the estimand, the solution is to
change the model, not to redefine the question.

In dose-response analysis, the distinction is especially consequential.
A pharmacologist interested in the \\\text{ED}\_{50}\\ (the dose
producing half-maximal response) has a different estimand from an
ecologist interested in the maximum sustainable yield (which depends on
\\V\_{\max}\\ relative to the organism’s metabolic rate). Both might use
the same Michaelis-Menten model, but they would evaluate model adequacy
by different criteria and might reach different conclusions about
whether M2 or M3 is preferable.

The estimand is defined at two levels. At the population level, the
quantities of interest are the population-average parameters
\\V\_{\max}^{(\text{pop})}\\ and \\K_m^{(\text{pop})}\\ and the
between-group variability \\\tau\_{V}\\ and \\\tau\_{K}\\. At the group
level, the group-specific parameters \\V\_{\max,g}\\ and \\K\_{m,g}\\
characterize how each population deviates from the population average.
In lucifer’s model contract, the estimand is encoded in two places: the
`parm.names` vector in the Data list names the parameters to be
estimated, and the `mon.names` vector names additional quantities to
track during sampling (derived parameters, predictions, likelihood
components).

### Design and target contexts

The design specifies \\N = 30\\ observations per group at doses evenly
spaced on \\(0.1, 10)\\, for a total of \\N\_{\text{tot}} = 120\\
observations. The choice of a balanced design with evenly spaced doses
simplifies the example but is not required by the methodology; lucifer’s
model contract makes no assumptions about the design structure.

The target contexts for prediction include new doses within the observed
range (interpolation), doses beyond the observed range (extrapolation,
where the saturation behavior becomes critical), and predictions for
unobserved groups drawn from the same population of groups (hierarchical
prediction). The distinction between interpolation and extrapolation is
important: a model that fits the observed data well may extrapolate
poorly if its functional form is wrong, and the prior predictive check
is the natural tool for detecting this kind of misspecification before
the data are seen.

## Synthetic data

Working with synthetic data where the ground truth is known is the
foundation of computational validation. This idea, formalized by Talts
et al. [\[2\]](#ref2) as simulation-based calibration and by Cook,
Gelman, and Rubin (2006) as the “self-recovering” property, provides the
strongest available test of the correctness of a Bayesian computation.
Any model that cannot recover known parameters from its own generative
process cannot be trusted to provide reliable inference from real data.
This principle applies regardless of model complexity: even simple
models should be validated against simulated data before being applied.

``` r

library(lucifer)
library(ggplot2)

set.seed(666)

# Design
G <- 4                                        # groups
N_per <- 30                                   # observations per group
doses <- seq(0.1, 10, length.out = N_per)

# Population-level truth
Vmax_pop <- 8.0;  Km_pop <- 2.0
tau_Vmax <- 1.0;  tau_Km <- 0.5

# Group-level truth (drawn from the hierarchical prior)
Vmax_true <- rnorm(G, Vmax_pop, tau_Vmax)
Km_true   <- pmax(rnorm(G, Km_pop, tau_Km), 0.1)
sigma_true <- 0.6
nu_true    <- 4                               # Student-t errors

# Assemble observations
group <- rep(1:G, each = N_per)
dose  <- rep(doses, G)
mu    <- Vmax_true[group] * dose / (Km_true[group] + dose)
y     <- mu + sigma_true * rt(G * N_per, df = nu_true)

# Collect ground truth
truth <- list(
    Vmax_pop = Vmax_pop, Km_pop = Km_pop,
    tau_Vmax = tau_Vmax, tau_Km = tau_Km,
    Vmax = Vmax_true, Km = Km_true,
    sigma = sigma_true, nu = nu_true
)
```

The data-generating process uses Student-\\t\\ errors with \\\nu = 4\\
degrees of freedom. This choice is deliberate: it produces occasional
extreme residuals that a Gaussian error model cannot accommodate,
creating a genuine modeling decision that the workflow must detect. The
three diagnostic iterations in this vignette (Sections 4, 7, and 8) each
discover a different inadequacy that the heavier-than-Gaussian tails
expose.

The scatter reveals two features that will drive the iterative workflow.
The four groups show visible differences in their saturation levels and
half-maximal doses, indicating that a single-population model (M1) will
systematically mispredict some groups. Several observations deviate
substantially from the mean curves, reflecting the heavy tails of the
\\t_4\\ distribution. These features are by construction, but in
practice they would emerge from exploratory data analysis before any
model is specified.

An important methodological point deserves emphasis. The decision to
generate synthetic data before fitting any model is not a preliminary
step to be rushed through. It is a fundamental component of the workflow
that serves multiple purposes simultaneously. First, it verifies that
the estimator can recover known parameters, establishing a necessary
condition for trusting real-data results. Second, it reveals the
sensitivity of estimation to sample size, design structure, and error
distribution, information that is invisible when only real data are
available. Third, it calibrates the analyst’s expectations: if the
posterior intervals for \\V\_{\max}\\ under the synthetic DGP are \\\pm
0.5\\, then intervals of \\\pm 5\\ on the real data suggest either a
radically different signal-to-noise ratio or a model problem. Gelman et
al. [\[1\]](#ref1) describe this as “fake-data simulation” and recommend
it as the first step in any serious analysis.

``` r

# Summary statistics by group
library(ggplot2)

df_summary <- do.call(rbind, lapply(1:G, function(g) {
    idx <- group == g
    data.frame(
        Group = factor(g),
        N = sum(idx),
        Mean_y = mean(y[idx]),
        SD_y = sd(y[idx]),
        Range_dose = paste0("[", round(min(dose[idx]), 1), ", ",
                            round(max(dose[idx]), 1), "]"),
        Vmax_true = Vmax_true[g],
        Km_true = Km_true[g]
    )
}))

knitr::kable(df_summary, digits = 2,
    caption = "Summary of synthetic data by group.",
    col.names = c("Group", "N", "Mean(y)", "SD(y)",
                  "Dose range", "True Vmax", "True Km"))
```

## Model specification

A central idea in Bayesian workflow is modular construction: the model
is assembled from replaceable components (likelihood, link function,
prior, hierarchical structure), and each component can be swapped
independently. Gelman et al. [\[1\]](#ref1) frame this as thinking of
“components as placeholders” that allow the analyst to go back and
generalize as needed. In lucifer, modularity is supported at two levels:
the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL compiles a probabilistic notation into a Model function that the
analyst can inspect and extend, and the hand-written Model interface
allows arbitrary R code with the only constraint being the five-element
return list.

The vignette develops four models, each extending the previous one by
replacing a single component. The first two use the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL; the last two use hand-written Model functions to demonstrate the
full flexibility of the interface when the DSL’s scope is exceeded.

### M0: linear regression (deliberately misspecified)

The simplest possible model ignores both the nonlinear dose-response
relationship and the group structure. Starting with a deliberately
misspecified model is a standard workflow technique (Section 5.2 of
[\[1\]](#ref1)) that serves two purposes. First, it provides a baseline
against which every subsequent improvement can be measured, making the
value of each model expansion quantifiable rather than qualitative.
Second, its failure at the prior predictive stage demonstrates that
prior predictive checking can detect fundamental misspecification before
any data are fit, saving the computational cost of fitting a model that
is known a priori to be wrong. In frequentist practice, the linear model
would typically be fit first and rejected based on residual diagnostics
after the fact; in the Bayesian workflow, the rejection occurs before
fitting, which is both faster and more principled.

``` r

spec_m0 <- model_spec("
    y ~ Normal(mu, sigma)
    mu = a + b * dose
    a ~ Normal(0, 10)
    b ~ Normal(0, 5)
    sigma ~ HalfCauchy(2.5)
")

Data_m0 <- spec_m0$data_template(y = y, dose = dose)
IV_m0   <- spec_m0$initial_values(Data_m0)
```

### M1: Michaelis-Menten (correct functional form, no groups)

The Michaelis-Menten equation captures the saturation behavior that the
linear model misses. It represents a genuine application of domain
knowledge: the analyst knows from biochemistry that the response must
saturate, and encodes this knowledge in the functional form of the mean
function rather than letting the data “discover” the nonlinearity
through a flexible basis. This is an example of what Box [\[7\]](#ref7)
calls using theory to constrain the model space, reducing the effective
dimensionality of the problem and improving identifiability. The
deterministic assignment `mu = Vmax * dose / (Km + dose)` defines the
nonlinear mean function, and the priors are chosen to be weakly
informative on the biologically relevant scale.

``` r

spec_m1 <- model_spec("
    y ~ Normal(mu, sigma)
    mu = Vmax * dose / (Km + dose)
    Vmax ~ Normal(10, 5)
    Km ~ Normal(2, 2)
    sigma ~ HalfCauchy(2.5)
")

Data_m1 <- spec_m1$data_template(y = y, dose = dose)
IV_m1   <- spec_m1$initial_values(Data_m1)
```

The prior on \\V\_{\max}\\ is centered at 10 with standard deviation 5,
covering a wide range of plausible maximum responses. The prior on
\\K_m\\ is centered at 2 with standard deviation 2. Both priors are
deliberately vague relative to the ground truth, ensuring that the
posterior is driven primarily by the data. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
system automatically enforces the positivity constraints implied by the
distribution support: \\V\_{\max}\\ and \\K_m\\ are constrained to
\\\[0, \infty)\\ via the
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)
function in the generated code.

### M2: hierarchical Michaelis-Menten

When the single-population model M1 fails posterior predictive checks
due to unmodeled group heterogeneity (Section 7), the natural expansion
is a hierarchical model with group-level parameters drawn from a
population distribution. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL does not currently support observation-level indexing into
group-level parameter vectors (e.g., `Vmax[group[i]]`), so M2 is written
as a hand-written Model function. This transition from declarative to
manual specification is common in practice: the DSL handles standard
models efficiently, and the manual interface provides unlimited
flexibility for everything else.

``` r

Model_m2 <- function(parm, Data) {
    ## -- Extract and constrain parameters -- ##
    G <- Data$G
    Vmax_pop <- parm[1]
    Km_pop   <- parm[2]
    tau_Vmax <- interval(parm[3], 1e-100, Inf);  parm[3] <- tau_Vmax
    tau_Km   <- interval(parm[4], 1e-100, Inf);  parm[4] <- tau_Km
    sigma    <- interval(parm[5], 1e-100, Inf);  parm[5] <- sigma
    Vmax_g   <- interval(parm[6:(5 + G)], 1e-100, Inf)
    parm[6:(5 + G)] <- Vmax_g
    Km_g     <- interval(parm[(6 + G):(5 + 2 * G)], 1e-100, Inf)
    parm[(6 + G):(5 + 2 * G)] <- Km_g

    ## -- Deterministic transformation -- ##
    mu <- Vmax_g[Data$group] * Data$dose / (Km_g[Data$group] + Data$dose)

    ## -- Log-likelihood -- ##
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))

    ## -- Log-prior -- ##
    # Population-level
    LP_prior <- dnormv(Vmax_pop, 10, 25, log = TRUE) +
                dnormv(Km_pop, 2, 4, log = TRUE) +
                dhalfcauchy(tau_Vmax, 2.5, log = TRUE) +
                dhalfcauchy(tau_Km, 2.5, log = TRUE) +
                dhalfcauchy(sigma, 2.5, log = TRUE)
    # Group-level
    LP_prior <- LP_prior +
        sum(dnorm(Vmax_g, Vmax_pop, tau_Vmax, log = TRUE)) +
        sum(dnorm(Km_g, Km_pop, tau_Km, log = TRUE))

    ## -- Assemble -- ##
    LP  <- LL + LP_prior
    yhat <- mu

    Modelout <- list(
        LP = LP, Dev = -2 * LL, Monitor = LP,
        yhat = yhat, parm = parm
    )
    return(Modelout)
}
```

The Data list for M2 requires the group indicator and the number of
groups, in addition to the dose and response vectors.

``` r

parm_names_m2 <- c("Vmax_pop", "Km_pop", "tau_Vmax", "tau_Km", "sigma",
                   paste0("Vmax_", 1:G), paste0("Km_", 1:G))

Data_m2 <- list(
    y = y, dose = dose, group = group, G = G,
    mon.names  = "LP",
    parm.names = parm_names_m2
)

IV_m2 <- c(8, 2, 1, 0.5, 1,              # population + sigma
           rep(8, G), rep(2, G))           # group-level
```

### M3: robust hierarchical (Student-t errors)

When M2 fails posterior predictive checks due to heavy-tailed residuals
(Section 8), the remedy is to replace the Gaussian likelihood with a
Student-\\t\\ distribution. The modified Model function adds a
degrees-of-freedom parameter \\\nu\\ and uses
[`dst()`](https://robustecologies.github.io/lucifer/reference/dist.Student.t.md)
for the log-likelihood.

``` r

Model_m3 <- function(parm, Data) {
    G <- Data$G
    Vmax_pop <- parm[1]
    Km_pop   <- parm[2]
    tau_Vmax <- interval(parm[3], 1e-100, Inf);  parm[3] <- tau_Vmax
    tau_Km   <- interval(parm[4], 1e-100, Inf);  parm[4] <- tau_Km
    sigma    <- interval(parm[5], 1e-100, Inf);  parm[5] <- sigma
    nu       <- interval(parm[6], 2, Inf);       parm[6] <- nu
    Vmax_g   <- interval(parm[7:(6 + G)], 1e-100, Inf)
    parm[7:(6 + G)] <- Vmax_g
    Km_g     <- interval(parm[(7 + G):(6 + 2 * G)], 1e-100, Inf)
    parm[(7 + G):(6 + 2 * G)] <- Km_g

    mu <- Vmax_g[Data$group] * Data$dose / (Km_g[Data$group] + Data$dose)

    LL <- sum(dst(Data$y, mu = mu, sigma = sigma, nu = nu, log = TRUE))

    LP_prior <- dnormv(Vmax_pop, 10, 25, log = TRUE) +
                dnormv(Km_pop, 2, 4, log = TRUE) +
                dhalfcauchy(tau_Vmax, 2.5, log = TRUE) +
                dhalfcauchy(tau_Km, 2.5, log = TRUE) +
                dhalfcauchy(sigma, 2.5, log = TRUE) +
                dexp(nu - 2, rate = 0.1, log = TRUE) +
                sum(dnorm(Vmax_g, Vmax_pop, tau_Vmax, log = TRUE)) +
                sum(dnorm(Km_g, Km_pop, tau_Km, log = TRUE))

    LP   <- LL + LP_prior
    yhat <- mu

    Modelout <- list(
        LP = LP, Dev = -2 * LL, Monitor = LP,
        yhat = yhat, parm = parm
    )
    return(Modelout)
}
```

``` r

parm_names_m3 <- c("Vmax_pop", "Km_pop", "tau_Vmax", "tau_Km", "sigma",
                   "nu", paste0("Vmax_", 1:G), paste0("Km_", 1:G))

Data_m3 <- list(
    y = y, dose = dose, group = group, G = G,
    mon.names  = "LP",
    parm.names = parm_names_m3
)

IV_m3 <- c(8, 2, 1, 0.5, 1, 5,           # population + sigma + nu
           rep(8, G), rep(2, G))
```

The prior on \\\nu\\ deserves comment. The reparameterization
`nu = interval(parm[6], 2, Inf)` constrains the degrees of freedom to be
at least 2, guaranteeing a finite variance. The prior \\\nu - 2 \sim
\text{Exponential}(0.1)\\ has a mean of 10 (so the prior mean of \\\nu\\
is 12), gently favoring moderate tail behavior while allowing the data
to push \\\nu\\ toward heavier tails if warranted. The true value \\\nu
= 4\\ lies well within the support of this prior.

The choice to constrain \\\nu \> 2\\ rather than \\\nu \> 0\\ is a
deliberate modeling decision. Student-\\t\\ distributions with \\\nu
\leq 1\\ (Cauchy and heavier) have no finite mean, and distributions
with \\\nu \leq 2\\ have no finite variance. For a dose-response
analysis where the mean and variance of the response are quantities of
interest, allowing \\\nu \leq 2\\ would create a model where the
estimand is undefined under the model’s own assumptions. The constraint
\\\nu \> 2\\ avoids this pathology while still allowing very heavy tails
(\\\nu = 2.1\\ has kurtosis \\\> 100\\). This kind of principled
constraint, motivated by the estimand rather than computational
convenience, exemplifies the interaction between the Theory and
Estimator nodes in the workflow diagram.

### Inspecting the compiled model

The `model_spec` object provides methods for inspection. The
[`code()`](https://robustecologies.github.io/lucifer/reference/code.md)
method reveals the generated R source for the Model function, which is
useful for learning and for extracting a starting point when manual
modifications are needed. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method renders
a directed acyclic graph of the dependency structure.

``` r

print(spec_m1)
code(spec_m1)
plot(spec_m1)
```

For computationally demanding models,
[`compile()`](https://robustecologies.github.io/lucifer/reference/compile.md)
generates C++ source from the intermediate representation and compiles
it via Rcpp, providing 10-100x speedups. The compiled model obeys the
same contract as the R-level model and is interchangeable in every
function that accepts a Model.

``` r

cspec_m1 <- compile(spec_m1, Data_m1)
```

## Prior predictive checking

Prior predictive checking simulates data from the model using only the
prior, without conditioning on observations. The resulting prior
predictive distribution, \\p(\tilde{y}) = \int p(\tilde{y} \mid
\theta)\\ p(\theta)\\ d\theta\\, reflects the joint implications of the
prior and the likelihood for observable quantities. If this distribution
places substantial mass on impossible outcomes (negative concentrations,
probabilities outside \\\[0,1\]\\, response values larger than any
physical system could produce), the prior is miscalibrated and should be
revised before any data are seen.

This step has no analogue in frequentist practice. A maximum likelihood
estimate is computed without any prior, so there is no prior predictive
distribution to check. The Bayesian requirement to specify a prior
creates both a burden (an additional modeling choice that must be
justified) and an opportunity (a tool for incorporating domain knowledge
and detecting misspecification before fitting). Simpson et
al. [\[11\]](#ref11) develop this idea further with penalized complexity
(PC) priors, which are constructed by specifying a base model (the
simplest scientifically defensible model) and penalizing deviations from
it. The PC prior framework provides a principled way to calibrate priors
that is more systematic than “weakly informative” recommendations alone.

This step corresponds to the bidirectional red arrow between the
Estimator and Prior_predictive nodes in the workflow diagram: the check
may send the analyst back to modify the model. The bidirectionality is
important. It is not enough to check the prior once; after modifying the
prior in response to a failed check, the analyst should re-check to
verify that the modification resolved the problem without introducing
new ones.

### M0: linear model fails prior predictive

The linear model `y = a + b \cdot d` predicts unbounded growth as dose
increases. Under the prior \\b \sim \text{Normal}(0, 5)\\, roughly half
of prior draws produce negative slopes (declining response with dose)
and the other half produce positive slopes that grow without bound.
Neither behavior is consistent with the known saturation of biological
dose-response relationships.

``` r

rprior_m0 <- function(n) {
    cbind(
        a     = rnorm(n, 0, 10),
        b     = rnorm(n, 0, 5),
        sigma = abs(rcauchy(n, 0, 2.5))
    )
}

prior_predictive_check(spec_m0$Model, Data_m0, rprior = rprior_m0, n = 300,
                       type = "density")
```

The density overlay reveals that the prior predictive distribution is
extremely diffuse, spanning values from \\-100\\ to \\+100\\ while the
observed data occupy a narrow range near \\(0, 10)\\. More critically,
the prior predictive mean function is linear, which cannot reproduce the
characteristic saturation plateau visible in the data. This is a prior
predictive failure driven not by the prior scale but by the functional
form of the model.

**Decision:** the linear model is abandoned in favor of M1. This is the
first pass through the Estimator \\\leftrightarrow\\ Prior_predictive
loop.

### M1: Michaelis-Menten passes prior predictive

The Michaelis-Menten model, with its asymptotic plateau at \\V\_{\max}\\
and inflection near \\K_m\\, produces prior predictive curves that are
qualitatively consistent with the data.

``` r

rprior_m1 <- function(n) {
    cbind(
        Vmax  = abs(rnorm(n, 10, 5)),
        Km    = abs(rnorm(n, 2, 2)),
        sigma = abs(rcauchy(n, 0, 2.5))
    )
}

# Density overlay
prior_predictive_check(spec_m1$Model, Data_m1, rprior = rprior_m1, n = 300,
                       type = "density")

# Ribbon: prior predictive dose-response bands
prior_predictive_check(spec_m1$Model, Data_m1, rprior = rprior_m1, n = 300,
                       type = "ribbon")

# Test statistic: does the prior predict a plausible range for sd(y)?
prior_predictive_check(spec_m1$Model, Data_m1, rprior = rprior_m1, n = 300,
                       type = "stat", stat_fun = sd)
```

The ribbon plot shows that the 95% prior predictive band covers the
observed data without being excessively wide. The test statistic check
confirms that the observed standard deviation of \\y\\ falls within the
central mass of the prior predictive distribution of `sd(y)`, indicating
that the prior on \\\sigma\\ is compatible with the variability in the
data. The prior is weakly informative: it excludes absurd values without
pulling inference toward any particular plausible value.

### M2: hierarchical prior predictive

For the hierarchical model, the prior predictive check must account for
the additional variability introduced by the group-level distribution.
The `rprior` function samples population-level parameters, then draws
group-level parameters from the hierarchical prior, and finally
evaluates the model.

``` r

rprior_m2 <- function(n) {
    Vmax_pop <- rnorm(n, 10, 5)
    Km_pop   <- rnorm(n, 2, 2)
    tau_V    <- abs(rcauchy(n, 0, 2.5))
    tau_K    <- abs(rcauchy(n, 0, 2.5))
    sig      <- abs(rcauchy(n, 0, 2.5))
    Vmax_g   <- matrix(NA, n, G)
    Km_g     <- matrix(NA, n, G)
    for (g in seq_len(G)) {
        Vmax_g[, g] <- pmax(rnorm(n, Vmax_pop, tau_V), 1e-2)
        Km_g[, g]   <- pmax(rnorm(n, Km_pop, tau_K), 1e-2)
    }
    cbind(Vmax_pop, Km_pop, tau_V, tau_K, sig, Vmax_g, Km_g)
}

prior_predictive_check(Model_m2, Data_m2, rprior = rprior_m2, n = 200,
                       type = "intervals")
```

The intervals plot shows, for each observation, the 50% and 95% prior
predictive intervals sorted by median. The observed values (dark points)
should fall within these bands if the prior is compatible with the data.
For M2, the prior predictive intervals are wider than those for M1 due
to the additional hierarchical variability, which is expected and
appropriate.

## Fake-data validation and synthetic estimates

Before fitting to real data, the model should be validated against data
where the truth is known. This step corresponds to the
Synthetic_estimates node in the workflow diagram: the synthetic sample
(Section 2) is fed through the estimator to produce estimates that can
be compared against the known parameters.

### The role of synthetic validation in the workflow

The relationship between synthetic validation and real-data inference is
subtle. A model that passes synthetic validation (recovering the true
parameters from data generated by its own DGP) has demonstrated
computational correctness but not model adequacy. The model could be
perfectly implemented yet scientifically wrong. Conversely, a model that
fails synthetic validation is almost certainly inappropriate for real
data: if it cannot recover parameters from its own assumptions, any
apparent fit to real data is coincidental. Synthetic validation is
therefore a necessary but not sufficient condition for trustworthy
inference, and it should be completed before any real data are analyzed.

This asymmetry explains why the workflow diagram places synthetic
validation (the Synthetic_sample \\\to\\ Synthetic_estimates pathway)
before the real-data fitting. The analyst proceeds to real data only
after establishing that the computational machinery works under
controlled conditions, following the principle of Talts et
al. [\[2\]](#ref2) that the correctness of the inference algorithm
should be verified independently of the correctness of the model.

### Algorithm recommendation

[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
profiles the model to characterize its computational properties and
recommends an inference strategy. It evaluates dimensionality, gradient
availability, constraint structure, evaluation speed, posterior
conditioning, and multimodality risk, then scores each algorithm family.

``` r

rx_m1 <- Prescribe(spec_m1$Model, Data_m1, IV_m1)
print(rx_m1)
```

For a 3-parameter nonlinear model with smooth gradients,
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
typically recommends NUTS as the primary method, with variational Bayes
and Laplace approximation as fast alternatives for exploration. The
recommendation includes ready-to-paste code for the top-ranked
algorithm.

### Fit fast, fail fast

An important intermediate goal in Bayesian workflow is the ability to
fail fast when fitting bad models. The time spent on (near-)perfect
inference for a model that will ultimately be abandoned is wasted.
Gelman et al. [\[1\]](#ref1) illustrate this with a compelling analogy:
if fitting a model to fake data fails, there is no point in fitting it
to real data. The same principle applies at a coarser level: if a cheap
approximation (Laplace, variational) reveals that the model is
fundamentally problematic (multimodal, non-identified, numerically
unstable), a full MCMC run will only confirm the problem more slowly and
expensively.

Before committing to a full MCMC run, a quick Laplace approximation can
verify that the model is computationally tractable and that the
posterior is at least unimodal near the starting values. If the Laplace
approximation fails to converge or produces a Hessian that is not
positive definite, the model almost certainly has structural problems
that MCMC will only obscure with chain-by-chain variability.

``` r

la_m1 <- LaplaceApproximation(spec_m1$Model, IV_m1, Data_m1,
                              Method = "SPG", Iterations = 500)
print(la_m1)
```

A variational Bayes fit provides a richer approximation and can reveal
multimodality or challenging geometry that the Laplace approximation
misses.

``` r

vb_m1 <- VariationalBayes(spec_m1$Model, IV_m1, Data_m1,
                           Covar = NULL, Iterations = 2000,
                           Method = "ADVI.mf")
print(vb_m1)
```

### Full MCMC fit

With the quick checks confirming tractability, the full MCMC fit can
proceed. NUTS (the No-U-Turn Sampler of Hoffman and Gelman, 2014) is the
default recommendation for differentiable models with moderate
dimensionality. It adapts the trajectory length automatically by
detecting U-turns in the Hamiltonian dynamics, eliminating the need to
tune the number of leapfrog steps that standard HMC requires. lucifer’s
implementation supports multiple chains via parallel PSOCK workers, with
the `Chains` and `CPUs` arguments controlling the parallelism. Running 4
chains is the standard recommendation from [\[3\]](#ref3), providing
enough between-chain variability for reliable \\\hat{R}\\ computation
while keeping the total wall-clock time manageable.

``` r

fit_m1 <- lucifer(spec_m1$Model, Data_m1, IV_m1,
                  Iterations = 20000, Algorithm = "NUTS",
                  Thinning = 10, Chains = 4, CPUs = 4)
```

### Parameter recovery

The caterpillar plot compares posterior estimates against the known
ground truth. For a correctly specified model fitted to data from its
own generative process, the true parameter values should fall within the
95% credible intervals.

``` r

gt_m1 <- c(Vmax = Vmax_pop, Km = Km_pop, sigma = sigma_true)
caterpillar.plot(fit_m1, ground_truth = gt_m1)
```

For M1, which ignores group structure, the population-level estimates
\\V\_{\max}\\ and \\K_m\\ will be close to the population-average truth
but with inflated uncertainty, because the unmodeled group heterogeneity
appears as excess residual variance.

### Calibration via coverage simulation

[`coverage_sim()`](https://robustecologies.github.io/lucifer/reference/coverage_sim.md)
provides a more rigorous check than a single parameter recovery plot. It
repeatedly simulates data from the DGP, fits the model, and computes the
empirical coverage of the credible intervals. If the model and inference
are well calibrated, the 95% intervals should cover the true values
approximately 95% of the time.

``` r

dgp_m1 <- function() {
    d <- seq(0.1, 10, length.out = 30)
    mu <- Vmax_pop * d / (Km_pop + d)
    y_sim <- mu + sigma_true * rnorm(30)
    list(y = y_sim, dose = d,
         mon.names = "LP",
         parm.names = c("Vmax", "Km", "sigma"))
}

fit_method_m1 <- function(Model, Data, IV) {
    lucifer(spec_m1$Model, Data_m1, IV_m1, Iterations = 5000,
            Algorithm = "NUTS", Thinning = 5, Chains = 3)
}

cov_m1 <- coverage_sim(
    dgp = dgp_m1,
    fit_methods = list(NUTS = fit_method_m1),
    true_parms = c(Vmax = Vmax_pop, Km = Km_pop, sigma = sigma_true),
    n_sims = 50, seed = 666, verbose = FALSE
)
print(cov_m1)
plot(cov_m1)
```

Simulation-based calibration (SBC; Talts et al. [\[2\]](#ref2)) provides
an even stronger test by checking that the rank of the true parameter
value within the posterior samples is uniformly distributed. lucifer
provides
[`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md)
for simulation-based inference objects; for general MCMC fits,
[`coverage_sim()`](https://robustecologies.github.io/lucifer/reference/coverage_sim.md)
tests the same calibration property through interval coverage, which is
the frequentist dual of the SBC rank uniformity condition.

## Convergence diagnostics and the refinement loop

Once a model is fitted, the first question is whether the computation
succeeded. This question is distinct from whether the model is correct:
convergence diagnostics assess the quality of the approximation to the
target posterior, not the quality of the posterior itself. A perfectly
converged MCMC chain for a misspecified model produces an accurate
representation of the wrong distribution, while a poorly converged chain
for a well-specified model produces an inaccurate representation of the
right one. Both situations are problematic, but they require different
remedies. Convergence diagnostics address the second; posterior
predictive checks and sensitivity analysis (Sections 7-8) address the
first.

Betancourt [\[10\]](#ref10) emphasizes that computational diagnostics
should not be an afterthought but an integral part of the workflow. A
model whose computation does not converge is providing no information,
and running it longer or with different random seeds is not a solution
when the problem is structural (the “folk theorem of statistical
computing” from Section 5.1 of [\[1\]](#ref1)). The right response to
persistent convergence failure is to investigate the model itself: is
the posterior multimodal? Is the geometry pathological (funnels,
ridges)? Are the data insufficient to identify some parameters? These
questions are addressed by the diagnostic tools in this section and by
the model modification tools in Section 8.

### Automated assessment with Consort

[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
evaluates a fitted object against a checklist of convergence criteria
tailored to the inference family. For MCMC, it checks the acceptance
rate, Monte Carlo standard error relative to posterior standard
deviation, bulk and tail effective sample sizes, split rank-normalized
\\\hat{R}\\, and (for NUTS) divergent transitions. If any criterion
fails,
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
returns a structured suggestion with a rationale, recommended action,
and ready-to-paste code.

``` r

dx_m1 <- Consort(fit_m1)
print(dx_m1)
is.appeased(fit_m1)
```

The
[`is.appeased()`](https://robustecologies.github.io/lucifer/reference/is.appeased.md)
function provides a boolean summary: `TRUE` if all criteria pass,
`FALSE` otherwise. It is the simplest entry point for automated
workflows where convergence must be verified programmatically.

### Individual diagnostics

When
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
raises concerns, the individual diagnostics help locate the source.

``` r

# Split rank-normalized R-hat (Vehtari et al. 2021)
Rhat(fit_m1)

# Bulk and tail effective sample size
ESS.bulk(fit_m1)
ESS.tail(fit_m1)

# Monte Carlo standard error as fraction of posterior SD
MCSE(fit_m1)

# Integrated autocorrelation time
IAT(fit_m1)

# Automatic burn-in detection via BMK stationarity test
burnin(fit_m1)
```

The split rank-normalized \\\hat{R}\\ of Vehtari et al. [\[3\]](#ref3)
is the modern replacement for the classical Gelman-Rubin diagnostic. It
splits each chain in half, rank-normalizes the combined draws, and
computes the between/within variance ratio. Values below 1.01 indicate
adequate mixing. The bulk ESS measures efficiency in the center of the
distribution, while the tail ESS measures efficiency in the tails, which
is critical for credible interval estimation. A minimum bulk ESS of 400
and tail ESS of 200 are the thresholds used by
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md).

### Visualization

MCMC diagnostic plots complement the numerical summaries by revealing
patterns that scalar summaries can miss: multimodal marginals,
correlations between parameters, and slow mixing in specific parameters.

``` r

# Trace plots: stationarity and mixing
plot(fit_m1, type = "trace")

# Posterior densities by chain
plot(fit_m1, type = "density")

# Pairwise scatter: posterior correlations
plot(fit_m1, type = "pairs")

# Rank histograms: mixing quality (should be uniform)
plot(fit_m1, type = "rank")

# Energy diagnostics: exploration efficiency
plot(fit_m1, type = "energy")
```

The rank histogram deserves emphasis. For well-mixed chains, the rank
distribution within each chain should be approximately uniform across
the combined sample. Deviations from uniformity indicate that some
chains are exploring different regions, which is a sign of poor mixing
that \\\hat{R}\\ may not detect when the deviation is subtle.

### The refinement loop: when diagnostics fail

To demonstrate the feedback arrow from Diagnostics to Estimator,
consider a deliberately bad fit: M2 with too few iterations and an
inappropriate algorithm.

``` r

fit_m2_bad <- lucifer(Model_m2, Data_m2, IV_m2,
                      Iterations = 1000, Algorithm = "RWM",
                      Thinning = 1, Chains = 1)

dx_bad <- Consort(fit_m2_bad)
```

[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
detects the convergence failure and returns a suggestion. The suggestion
includes a rationale explaining why the current fit is inadequate, a
recommended action (e.g., increase iterations, switch algorithm, add
chains), and ready-to-paste code implementing the recommendation. The
analyst applies the suggestion and refits.

``` r

# The suggestion from Consort typically recommends:
# 1. Switch to a gradient-based algorithm (NUTS or HMCDA)
# 2. Increase iterations
# 3. Run multiple chains

fit_m2 <- lucifer(Model_m2, Data_m2, IV_m2,
                  Iterations = 3000, Algorithm = "NUTS",
                  Thinning = 15, Chains = 4, CPUs = 4)

dx_m2 <- Consort(fit_m2)
is.appeased(fit_m2)
```

This cycle of fit, diagnose, modify, refit is the core of the
Diagnostics \\\to\\ Estimator arrow. In practice, it may take several
iterations to find an algorithm and tuning configuration that produces
adequate convergence.

### NUTS-specific diagnostics

For NUTS fits,
[`check_nuts()`](https://robustecologies.github.io/lucifer/reference/check_nuts.md)
provides additional diagnostics: divergent transitions (indicating
regions where the leapfrog integrator’s error becomes large), tree depth
saturation (indicating the integrator is hitting its maximum depth), and
energy statistics.

``` r

check_nuts <- check_nuts(fit_m2)
print(check_nuts)
plot(check_nuts)
```

Divergent transitions are the most informative diagnostic. A single
divergent transition may be benign, but a pattern of divergences
concentrated in a specific region of parameter space indicates
problematic geometry (funnels, ridges, or cusps) that should be resolved
through reparameterization or model modification rather than by tuning
the sampler.

### Automated pipeline with Crucible

The
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
function automates the entire fit-diagnose-refine loop. It calls
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
to profile the model, selects the top-ranked algorithms, fits each with
multiple chains, evaluates convergence with
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md),
applies suggestions and refits when needed, and finally compares all
surviving fits with
[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md).
The result is a ranked list of methods with convergence metadata.

``` r

cruc_m2 <- Crucible(Model_m2, Data_m2, IV_m2,
                    n_methods = 5, max_rounds = 3,
                    Chains = 4, CPUs = 4,
                    families = c("MCMC", "VB", "Laplace"))
print(cruc_m2)
plot(cruc_m2)
```

[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
is the “push-button” approach to inference when the analyst wants to
explore multiple methods without manually managing the fit-diagnose
loop. It embodies the idea from [\[1\]](#ref1) that “modeling as
software development” requires both a manual craft mode (writing models,
inspecting diagnostics, making decisions) and an automated testing mode
(running suites of algorithms, checking convergence, comparing results).
The manual approach demonstrated above provides more control and deeper
understanding;
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
provides speed and breadth. In practice, the two approaches are
complementary:
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
quickly identifies which algorithm families work for a given model, and
the analyst then tunes the best performer manually for the final
production run.

The parallel to software testing is apt. A software engineer does not
skip unit tests because integration tests exist, nor does an integration
test replace understanding the code. Similarly,
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
does not replace the analyst’s judgment about model adequacy (which
requires posterior predictive checks and domain knowledge), but it does
replace the tedious process of manually trying half a dozen algorithms
and comparing their convergence properties.

## Posterior predictive checking

Posterior predictive checking, introduced by Rubin [\[8\]](#ref8) and
developed systematically by Gelman, Meng, and Stern (1996), uses the
fitted model to generate replicated datasets \\y^{\text{rep}} \sim
p(y^{\text{rep}} \mid y)\\ and compares them to the observed data. The
idea is elegantly simple: if the model is adequate, data simulated from
the posterior predictive distribution should be “similar” to the
observed data in ways that the analyst cares about. Systematic
discrepancies indicate model misfit, and the specific nature of the
discrepancy (which test statistic fails, which observations are poorly
predicted, which subgroups are systematically biased) provides
actionable guidance for model modification.

The power of posterior predictive checking lies in the choice of test
quantities. A single global discrepancy (e.g., the deviance) can detect
gross misfit but cannot localize the problem. Subgroup-specific checks
(e.g., predictions by group, by covariate level, or by spatial location)
can identify the specific aspect of the data that the model fails to
capture. Gabry et al. [\[6\]](#ref6) provide a comprehensive taxonomy of
visual PPC methods and demonstrate how each one targets a different type
of misfit. lucifer implements 40 PPC styles that cover this taxonomy,
from marginal density overlays to calibration checks via LOO-PIT.

This step corresponds to the bidirectional red arrow between the
Estimator and Post_predictive nodes: if the check fails, the model must
be modified.

### Generating posterior predictive samples

The [`predict()`](https://rdrr.io/r/stats/predict.html) method evaluates
the Model function at each posterior sample to obtain a matrix of
predicted values \\\hat{y}\\ with dimensions \\N\_{\text{obs}} \times
S\\, where \\S\\ is the number of posterior samples. Each column of this
matrix is a complete replicated dataset generated under one posterior
draw of the parameters: the \\s\\-th column contains \\\hat{y}^{(s)} =
f(\theta^{(s)}, X)\\ where \\\theta^{(s)}\\ is the \\s\\-th posterior
sample and \\f\\ is the Model function’s `yhat` output. The collection
of these replicated datasets constitutes the posterior predictive
distribution, from which any desired summary statistic, interval, or
visualization can be derived.

``` r

ppc_m1 <- predict(fit_m1, spec_m1$Model, Data_m1)
```

### Visual checks

lucifer provides 40 PPC plot styles through
[`plot.demonoid.ppc()`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md).
The most informative for the present analysis are the density overlay,
fitted-vs-observed scatter, predictive intervals, and test statistics.

``` r

# Density overlay: marginal distribution of y vs. y_rep
plot(ppc_m1, Style = "Density Overlay")
```

The density overlay superposes the density of each replicated dataset
(light lines) on the observed density (dark line). For M1, the overall
shape should be reasonable since the Michaelis-Menten function captures
the saturation behavior, but the width and shape of the replicated
densities may differ from the observed if the group structure or heavy
tails are not captured.

``` r

# Fitted vs. observed: posterior predictive median vs. y
plot(ppc_m1, Style = "Fitted")
```

The fitted-vs-observed scatter reveals systematic prediction errors.
Points far from the 1:1 line indicate observations that the model
consistently mispredicts.

``` r

# Predictive intervals: per-observation 50% and 95% bands
plot(ppc_m1, Style = "Intervals")
```

``` r

# Test statistic: is the model capturing the variance structure?
plot(ppc_m1, Style = "Stat", stat_fun = sd)
```

### Detecting group heterogeneity

The critical PPC for M1 is a grouped check. By passing a `Group`
argument, the predictive distribution is split by group, revealing
whether the single-population model systematically over- or
underpredicts specific groups.

``` r

plot(ppc_m1, Style = "Violin Grouped", Data = Data_m1,
     Group = factor(group))
```

The grouped violin plot shows that the posterior predictive distribution
for each group is centered at the same location (the population-average
prediction), while the observed data for each group are shifted relative
to this center. Groups with higher true \\V\_{\max}\\ are systematically
underpredicted; groups with lower true \\V\_{\max}\\ are systematically
overpredicted. This pattern is the signature of unmodeled group
heterogeneity and motivates the expansion to M2.

### LOO-PIT calibration

The leave-one-out probability integral transform (LOO-PIT) provides what
is arguably the most rigorous calibration check available for Bayesian
predictive distributions. The PIT of observation \\y_i\\ is defined as
\\u_i = P(y^{\text{rep}}\_i \leq y_i \mid y\_{-i})\\, the cumulative
predictive probability evaluated at the observed value, computed using
the leave-one-out predictive distribution. Under a well-calibrated
predictive distribution, the PIT values should be uniformly distributed
on \\(0, 1)\\. This is a consequence of the probability integral
transform theorem: if \\Y\\ has CDF \\F\\, then \\F(Y) \sim
\text{Uniform}(0, 1)\\.

Deviations from uniformity have specific interpretations. U-shaped PIT
histograms (excess mass near 0 and 1) indicate underdispersion: the
model’s predictive intervals are too narrow and observations fall in the
tails more often than expected. Inverse-U shaped distributions (excess
mass near 0.5) indicate overdispersion: the model’s predictive intervals
are too wide. Skewed distributions indicate systematic bias in the
predictions. Each pattern suggests a different remedy: underdispersion
may require heavier-tailed errors or additional variance components;
overdispersion may indicate overparameterization; bias may indicate
missing covariates or a wrong functional form. Gabry et
al. [\[6\]](#ref6) provide a comprehensive visual guide to interpreting
PIT diagnostics.

``` r

loo_m1 <- LOO(fit_m1)
plot(ppc_m1, Style = "LOO-PIT", loo = loo_m1)
```

### Posterior predictive summary

The [`summary()`](https://rdrr.io/r/base/summary.html) method computes
18 discrepancy measures including the Bayesian predictive information
criterion, L-criterion, concordance correlation, discrepancy statistics,
and replicate-level p-values.

``` r

summary(ppc_m1)
```

### PPC for M2: detecting heavy tails

After fitting M2 with the hierarchical structure, the grouped
heterogeneity is resolved, but a new problem appears: the Gaussian error
distribution cannot accommodate the occasional extreme residuals
produced by the Student-\\t\\ DGP.

``` r

ppc_m2 <- predict(fit_m2, Model_m2, Data_m2)

# The grouped violin now shows group-specific predictions
plot(ppc_m2, Style = "Violin Grouped", Data = Data_m2,
     Group = factor(group))

# But the stat check with max reveals heavy tails
plot(ppc_m2, Style = "Stat", stat_fun = function(x) max(abs(x)))
```

The test statistic `max(|y|)` is sensitive to extreme values. If the
observed maximum exceeds the bulk of the replicated maxima, the model’s
tails are too thin, specifically, the Gaussian error distribution
assigns too little probability to the extreme observations that the
Student-\\t\\ DGP produces. The observed `max(|y|)` falls in the upper
tail of the prior predictive distribution of the same statistic,
yielding a Bayesian p-value below 0.05. This motivates the expansion
from M2 to M3.

## Model modification and sensitivity analysis

The iterative workflow has now produced two failures that motivate model
expansion: M1 cannot capture group heterogeneity, and M2 cannot capture
heavy-tailed residuals. Model M3, the robust hierarchical model with
Student-\\t\\ errors, addresses both issues. Before fitting M3, a
sensitivity analysis on M2 provides additional insight into which
aspects of the inference are driven by the prior versus the data.

Sensitivity analysis occupies a distinctive position in the Bayesian
workflow. Unlike prior predictive and posterior predictive checks, which
test whether the model is consistent with domain knowledge and data
respectively, sensitivity analysis tests whether the conclusions are
robust to the analyst’s choices. Berger [\[9\]](#ref9) distinguishes
between “local” sensitivity (small perturbations of the prior) and
“global” sensitivity (comparison across a class of priors), arguing that
both are necessary for a thorough assessment. lucifer’s
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
framework implements both types: power-scaling provides local
sensitivity by raising the prior to a power \\\alpha\\ and tracking how
the posterior shifts, while cross-model comparison provides global
sensitivity by fitting the model under different prior specifications
and measuring divergence.

The practical value of sensitivity analysis is that it separates
conclusions that are supported by the data from conclusions that are
driven by the prior. In a well-identified model with informative data,
the posterior should be insensitive to reasonable perturbations of the
prior. When sensitivity is detected, the analyst faces a genuine
decision: either tighten the prior (adding information) or acknowledge
the uncertainty (reporting the range of conclusions). Both responses are
scientifically honest; what is dishonest is ignoring sensitivity and
presenting a single posterior as if it were definitive.

### Prior-to-posterior learning

The
[`prior_vs_posterior()`](https://robustecologies.github.io/lucifer/reference/prior_vs_posterior.md)
function overlays prior and posterior densities for selected parameters,
revealing how much the data have informed each parameter. Parameters
whose posterior closely resembles their prior are weakly identified by
the data; parameters whose posterior has contracted substantially are
strongly identified.

``` r

rprior_m2_fn <- function(n) {
    cbind(
        Vmax_pop = rnorm(n, 10, 5),
        Km_pop   = rnorm(n, 2, 2),
        tau_Vmax = abs(rcauchy(n, 0, 2.5)),
        tau_Km   = abs(rcauchy(n, 0, 2.5)),
        sigma    = abs(rcauchy(n, 0, 2.5))
    )
}

gt_m2 <- c(Vmax_pop = Vmax_pop, Km_pop = Km_pop,
           tau_Vmax = tau_Vmax, tau_Km = tau_Km,
           sigma = sigma_true)

prior_vs_posterior(fit_m2, prior = rprior_m2_fn,
                   Parms = c("Vmax_pop", "Km_pop", "tau_Vmax",
                             "tau_Km", "sigma"),
                   ground_truth = gt_m2)
```

### Prior sensitivity

[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md)
compares posteriors obtained under different prior specifications for
the same model and data. This is the visual counterpart of formal
sensitivity analysis: if the posteriors are similar across priors, the
inference is data-driven; if they differ, the prior is influential.

``` r

# Refit M2 with tighter hyperpriors
IV_m2_tight <- IV_m2
fit_m2_tight <- lucifer(Model_m2, Data_m2, IV_m2_tight,
                        Iterations = 30000, Algorithm = "NUTS",
                        Thinning = 15, Chains = 4, CPUs = 4)

prior_sensitivity(list("Wide priors" = fit_m2,
                       "Tight priors" = fit_m2_tight),
                  Parms = c("Vmax_pop", "Km_pop", "sigma"),
                  ground_truth = gt_m2)
```

### RobustBayes: comprehensive sensitivity

[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
provides a unified framework for five complementary sensitivity
analyses. Three modules are particularly relevant at this stage of the
workflow.

``` r

# Power-scaling sensitivity: which parameters are prior-driven?
rb_m2 <- RobustBayes(fit_m2, Model_m2, Data_m2, modules = "power")
print(rb_m2)
plot(rb_m2, type = "power_density")
```

The power-scaling analysis (Kallioinen et al. [\[4\]](#ref4)) raises the
prior to a power \\\alpha\\ and reweights the posterior via PSIS
importance sampling. Parameters whose posterior shifts substantially as
\\\alpha\\ varies are prior-sensitive. For the hierarchical model, the
group-level variance parameters \\\tau_V\\ and \\\tau_K\\ are often the
most prior-sensitive, because the data provide limited information about
between-group variability when \\G\\ is small.

``` r

# Observation influence: which data points drive the posterior?
rb_inf <- RobustBayes(fit_m2, Model_m2, Data_m2, modules = "influence")
plot(rb_inf, type = "influence")
```

The influence module uses PSIS-LOO importance weights to estimate the
posterior shift when each observation is removed. High-influence
observations (Pareto-\\k \> 0.7\\) warrant inspection: they may be
genuine extreme values (as in this case, where they arise from the
Student-\\t\\ DGP), data errors, or indicators of model
misspecification. In this analysis, the high-influence observations are
the ones with the largest residuals, which are exactly the observations
that motivate the switch from Gaussian to Student-\\t\\ errors.

``` r

# Prior-data conflict: are the priors compatible with the data?
Data_m2_pgf <- Data_m2
Data_m2_pgf$PGF <- function(Data) {
    parm <- rprior_m2_fn(1)
    mu <- parm[1] * Data$dose / (parm[2] + Data$dose)
    rnorm(length(mu), mu, parm[5])
}

rb_con <- RobustBayes(fit_m2, Model_m2, Data_m2_pgf, modules = "conflict")
print(rb_con)
plot(rb_con)
```

### Fitting M3

With the evidence from PPC (heavy tails) and `RobustBayes`
(high-influence extreme observations), the expansion to Student-\\t\\
errors is motivated. M3 is fitted and verified.

``` r

fit_m3 <- lucifer(Model_m3, Data_m3, IV_m3,
                  Iterations = 40000, Algorithm = "NUTS",
                  Thinning = 20, Chains = 4, CPUs = 4)

# Convergence check
dx_m3 <- Consort(fit_m3)
print(dx_m3)
is.appeased(fit_m3)
```

``` r

# PPC: the heavy-tail problem should be resolved
ppc_m3 <- predict(fit_m3, Model_m3, Data_m3)
plot(ppc_m3, Style = "Stat", stat_fun = function(x) max(abs(x)))
plot(ppc_m3, Style = "Density Overlay")
plot(ppc_m3, Style = "Violin Grouped", Data = Data_m3,
     Group = factor(group))
```

The `max(|y|)` statistic now falls comfortably within the replicated
distribution, confirming that the Student-\\t\\ errors accommodate the
extreme observations. The grouped violin plot shows appropriate
group-specific predictions. The model passes all posterior predictive
checks, completing the third and final pass through the diagnostic loop.

## Identifiability and data cloning

Before declaring a model adequate, the workflow should address a
question that sensitivity analysis alone cannot fully resolve: are the
parameters structurally identifiable from the data? A parameter is
identifiable if the likelihood function has a unique maximum; it is
estimable if the posterior concentrates on that maximum as the data
become increasingly informative. Sensitivity to the prior, detected by
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
in the previous section, can arise from two distinct causes: (a) the
data genuinely provide little information about the parameter (weak
identification), or (b) the parameter is structurally non-identifiable
under the model, so no amount of data would resolve the ambiguity.
Distinguishing these cases is critical because the remedies differ: weak
identification calls for more data or a more informative prior, while
structural non-identifiability calls for model reparameterization or
simplification.

Data cloning (Lele, Dennis and Lutscher 2007 [\[12\]](#ref12)) provides
a formal diagnostic for this distinction. The technique raises the
likelihood to the \\K\\-th power, \\\pi_K(\theta \mid y) \propto
\[L(\theta \mid y)\]^K \\ \pi(\theta)\\, which is mathematically
equivalent to observing \\K\\ independent copies of the dataset. As \\K
\to \infty\\, the posterior mean converges to the MLE, the scaled
posterior covariance \\K \cdot \text{Var}\_K(\theta)\\ converges to the
inverse Fisher information, and the posterior becomes asymptotically
normal. Three diagnostics emerge from this construction. First, for
identifiable parameters, the posterior variance must shrink at rate
\\1/K\\; failure to shrink indicates non-identifiability. Second, the
scaled variance \\K \cdot \text{Var}\_K(\theta)\\ must stabilize;
oscillation indicates that the asymptotic regime has not been reached.
Third, when the analysis is repeated under multiple distinct priors
(Campbell and Lele 2014 [\[13\]](#ref13)), identifiable parameters yield
posterior means that converge to the same MLE regardless of the prior,
while non-identifiable parameters remain prior-dependent even at large
\\K\\.

This step corresponds to the DataClone node in the workflow diagram. It
connects bidirectionally to the Estimator (the same model and sampler
are used, with a likelihood-power wrapper) and feeds into Diagnostics
(the convergence and estimability results inform model modification). In
the broader workflow, data cloning occupies a unique position: it is the
only diagnostic that simultaneously addresses prior sensitivity (by
checking whether conclusions depend on the prior at large \\K\\) and
structural identifiability (by checking whether the posterior
concentrates at all). Neither
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
nor
[`prior_sensitivity()`](https://robustecologies.github.io/lucifer/reference/prior_sensitivity.md)
can detect structural non-identifiability, because they perturb the
prior while holding \\K = 1\\ fixed; data cloning varies \\K\\ while
optionally varying the prior, separating the two effects.

### Data cloning on M2: hierarchical identifiability

The hierarchical model M2 has 11 parameters, including the between-group
variance components \\\tau_V\\ and \\\tau_K\\ that
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
flagged as prior-sensitive. The question is whether this sensitivity
reflects weak identification (too few groups, \\G = 4\\, to pin down the
variance components) or structural non-identifiability (the model cannot
distinguish \\\tau_V\\ from \\\sigma\\ given the data structure). Data
cloning answers this directly.

``` r

dc_m2 <- DataCloning(Model_m2, Data_m2, IV_m2,
                     n_clones = c(1, 2, 4, 8, 16),
                     Algorithm = "NUTS",
                     Iterations = 10000, Thinning = 10,
                     Chains = 3, CPUs = 3,
                     n_priors = 3,
                     update_prior = TRUE,
                     verbose = TRUE)

print(dc_m2)
summary(dc_m2)
```

The `n_priors = 3` argument triggers the Campbell-Lele estimability
test. Since no `prior_generator` is supplied,
[`DataCloning()`](https://robustecologies.github.io/lucifer/reference/DataCloning.md)
automatically creates location-shifted alternative priors from the \\K =
1\\ posterior: the alternatives are centered at \\\bar{\theta}\_1 \pm 3
\\ \text{SD}\_1\\ with standard deviation equal to the \\K = 1\\
posterior SD. This ensures that the alternative priors genuinely pull
the posterior away from the MLE at low \\K\\, creating a detectable
prior effect that should disappear at high \\K\\ for estimable
parameters. The `update_prior = TRUE` flag enables warm starting,
passing the posterior mean and covariance from clone \\K\\ as the
initial values and proposal covariance for clone \\K+1\\, which
substantially improves sampling efficiency at high \\K\\ where the
posterior becomes very narrow.

### Interpreting the diagnostics

``` r

# Individual diagnostic plots
plot(dc_m2, type = "convergence")   # posterior means across K
plot(dc_m2, type = "eigenvalue")    # largest eigenvalue ratio vs 1/K
plot(dc_m2, type = "scaled_var")    # K * Var(theta) stabilization
plot(dc_m2, type = "density")       # posterior density evolution across K
```

The convergence plot shows how each parameter’s posterior mean evolves
as \\K\\ increases. For well-identified parameters
(\\V\_{\max}^{(\text{pop})}\\, \\K_m^{(\text{pop})}\\, \\\sigma\\), the
means stabilize quickly and converge toward the MLE. For the variance
components (\\\tau_V\\, \\\tau_K\\), convergence may be slower,
reflecting the weaker identification that comes from having only \\G =
4\\ groups, but the means should still converge.

The scaled variance plot is the most informative diagnostic. Each panel
shows \\K \cdot \text{Var}\_K(\theta_j)\\ against \\K\\. Stabilization
(a horizontal line) confirms that the asymptotic regime has been reached
and the parameter is identifiable. A parameter whose scaled variance
continues to grow with \\K\\ is non-identifiable; one whose scaled
variance shrinks toward zero faster than \\1/K\\ indicates a boundary
issue (the MLE lies at the boundary of the parameter space, as can
happen with variance components near zero).

The eigenvalue diagnostic computes the standardized largest eigenvalue
ratio \\\lambda_S(K) = \lambda\_{\max}(\Sigma_K) /
\lambda\_{\max}(\Sigma_1)\\, which should decrease at rate \\1/K\\ for
identifiable models. Departures from this rate indicate that some linear
combination of parameters is poorly identified, even when individual
parameters appear stable.

### Estimability across priors

When `n_priors > 1`, the
[`summary()`](https://rdrr.io/r/base/summary.html) method reports the
Campbell-Lele ANOVA results for each parameter. The test classifies
parameters as “estimable” (posterior mean converges to the same value
across priors), “weakly estimable” (convergence is detectable but slow),
or “non-estimable” (posterior mean depends on the prior even at high
\\K\\).

``` r

# Estimability classification
summary(dc_m2)$estimability

# Profile likelihood for key parameters (requires Model and Data)
plot(dc_m2, type = "profile", Model = Model_m2, Data = Data_m2,
     n_grid = 80)
```

For the dose-response model, we expect all parameters to be estimable
because the Michaelis-Menten function is identifiable given sufficient
dose coverage and the hierarchical structure is identified when \\G \>
1\\. The variance components \\\tau_V\\ and \\\tau_K\\ should be
classified as estimable but with slower convergence, confirming that the
prior sensitivity detected by
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
reflects weak rather than structural non-identifiability. This
distinction has practical consequences: it means that collecting data
from more groups would reduce the prior dependence, whereas if the
parameters were structurally non-identifiable, additional groups would
not help.

### Data cloning on M3: verifying the final model

A brief data cloning check on M3 confirms that the addition of the
degrees-of-freedom parameter \\\nu\\ does not introduce identifiability
problems.

``` r

dc_m3 <- DataCloning(Model_m3, Data_m3, IV_m3,
                     n_clones = c(1, 2, 4, 8),
                     Algorithm = "NUTS",
                     Iterations = 10000, Thinning = 10,
                     Chains = 3, CPUs = 3,
                     update_prior = TRUE,
                     verbose = TRUE)

print(dc_m3)
plot(dc_m3, type = "scaled_var")
```

The scaled variance should stabilize for all 12 parameters of M3,
including \\\nu\\. The degrees-of-freedom parameter is typically
well-identified when the data contain genuine outliers (as they do here,
with Student-\\t_4\\ errors), because the tails of the data distribution
provide strong information about the tail weight. In models where the
errors are close to Gaussian (\\\nu \to \infty\\), the
degrees-of-freedom parameter becomes weakly identified and the scaled
variance convergence is slower, which would be visible in the data
cloning diagnostic.

### From data cloning to the MLE: a Bayesian bridge to frequentist inference

An elegant by-product of data cloning is that it provides MLE estimates
and asymptotic standard errors through Bayesian MCMC, without requiring
direct optimization of the likelihood. The posterior mean at high \\K\\
converges to the MLE, and \\K \cdot \text{Var}\_K(\theta)\\ converges to
the inverse Fisher information. This connection is useful in several
ways. First, it allows comparison of the Bayesian posterior (at \\K =
1\\) with the frequentist MLE (at high \\K\\), revealing where the prior
has a substantive effect on the inference. Second, it provides a
robustness check: if the Bayesian estimates at \\K = 1\\ differ
substantially from the data cloning MLE, the prior is influential and
the analyst should investigate whether this influence is scientifically
justified. Third, for models where direct MLE computation is difficult
(hierarchical models, models with latent variables), data cloning
provides a computationally straightforward alternative that reuses the
same MCMC infrastructure already employed for Bayesian inference.

``` r

# Compare Bayesian posterior (K=1) with data cloning MLE
cat("Data cloning MLE:\n")
print(dc_m2$MLE)

cat("\nAsymptotic SE:\n")
print(dc_m2$SE)

cat("\nWald 95% CI:\n")
print(dc_m2$Wald.CI)

cat("\nBayesian posterior mean (K=1):\n")
print(colMeans(fit_m2$Posterior2))
```

## Model comparison

With three fitted models (M1, M2, M3), the next step is systematic
comparison. The goal is not simply to pick the “best” model but to
understand how inferences change across models and to assess the
robustness of conclusions to modeling choices. Gelman et
al. [\[1\]](#ref1) frame this as exploring a “topology of models” rather
than performing a tournament with a single winner. A conclusion that
holds across M1, M2, and M3 is more credible than one that depends on
M3’s specific distributional assumptions. A conclusion that reverses
between M2 and M3 reveals that the error distribution matters for the
estimand, which is itself a valuable finding.

Vehtari et al. [\[5\]](#ref5) provide the theoretical foundation for
cross-validation-based model comparison in the Bayesian context, showing
that LOO-PSIS provides a computationally efficient and statistically
principled estimate of out-of-sample predictive accuracy. Unlike
marginal likelihood-based Bayes factors, which are sensitive to the
prior and compare models on their ability to predict the data before it
is observed, cross-validation compares models on their ability to
predict each observation after seeing all the others. This difference is
consequential: cross-validation rewards models that predict well for the
current data, while Bayes factors reward models that are probable under
the prior.

This section corresponds to the Estimates node in the workflow diagram,
where inferences from multiple models are compared.

### LOO-PSIS and WAIC

Leave-one-out cross-validation via Pareto-smoothed importance sampling
(LOO-PSIS; Vehtari et al. [\[5\]](#ref5)) estimates the out-of-sample
predictive accuracy without refitting. The WAIC provides a similar
estimate using a different decomposition.

``` r

loo_m1 <- LOO(fit_m1)
loo_m2 <- LOO(fit_m2)
loo_m3 <- LOO(fit_m3)

waic_m1 <- WAIC(fit_m1)
waic_m2 <- WAIC(fit_m2)
waic_m3 <- WAIC(fit_m3)

print(loo_m3)
plot(loo_m3)
```

The `plot.loo()` method visualizes the Pareto-\\k\\ diagnostic for each
observation. Observations with \\k \> 0.7\\ are poorly approximated by
importance sampling and may require exact cross-validation via
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md).
The proportion and distribution of high-\\k\\ observations also serves
as a model criticism tool: a model that poorly predicts many
observations has worse predictive performance than one that occasionally
struggles with a few.

### Ranked comparison

[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
ranks models by expected log pointwise predictive density (ELPD) and
reports the difference and its standard error for each pair. The best
model appears first; a difference larger than twice its standard error
provides moderate evidence of meaningful predictive superiority. The
standard error accounts for the correlation between pointwise LOO
contributions, which arises because the same data are used to fit and
evaluate the model (with one observation held out each time). This
correlation-adjusted standard error is essential for honest model
comparison: naive standard errors can be orders of magnitude too small,
leading to overconfident model selection. Sivula, Magnusson, and Vehtari
(2020) provide further discussion of the statistical properties of
LOO-based model comparison and the conditions under which the ELPD
difference is approximately normally distributed.

``` r

comp <- loo_compare(M1 = loo_m1, M2 = loo_m2, M3 = loo_m3)
print(comp)
plot(comp)
```

For this analysis, M3 should rank first (correctly specified), followed
by M2 (correct structure but wrong error distribution), and M1 (missing
group structure and wrong errors). The magnitudes of the differences
reflect the severity of each misspecification.

### K-fold cross-validation

When Pareto-\\k\\ diagnostics flag too many observations, exact
cross-validation via
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)
provides a more reliable estimate at the cost of refitting the model
\\K\\ times.

``` r

kf_m2 <- Kfold(Model_m2, Data_m2, fit_m2, K = 10)
print(kf_m2)
```

### Bayes factors

[`BayesFactor()`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md)
computes Bayes factor matrices from log marginal likelihood estimates.
Unlike LOO and WAIC, Bayes factors compare models on the basis of their
marginal likelihoods, which penalize complexity through the prior
predictive probability of the observed data.

``` r

bf <- BayesFactor(fit_m1, fit_m2, fit_m3)
print(bf)
```

### Model averaging and stacking

When no single model is clearly superior, model averaging combines
predictions across models using weights that reflect each model’s
predictive accuracy.
[`stacking_weights()`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md)
implements three weighting schemes: Bayesian stacking (which minimizes
cross-validation error), pseudo-BMA (which uses LOO-based marginal
likelihood proxies), and pseudo-BMA+ (which adds Bayesian bootstrap
correction).

``` r

w <- stacking_weights(loo_m1, loo_m2, loo_m3)
print(w)
```

Stacking is recommended over BMA in M-open settings where the true model
is not in the candidate set (Yao et al. 2018). In an M-closed setting,
where the true model is one of the candidates, BMA is theoretically
optimal. In practice, the M-open case is overwhelmingly more common:
real data are generated by processes more complex than any model in the
candidate set, and the question is not “which model is true?” but “which
model predicts best?” Stacking answers the latter question by finding
the convex combination of predictive distributions that minimizes
leave-one-out cross-validation error. In this vignette the true model
(M3, or rather its DGP) is approximately in the candidate set, so
stacking and BMA should agree in assigning most weight to M3.

### Cross-method benchmarking with Arena

[`Arena()`](https://robustecologies.github.io/lucifer/reference/Arena.md)
compares the performance of different inference methods on the same
model, complementing the cross-model comparison above. It evaluates
efficiency (ESS per second), accuracy (KL divergence and Wasserstein
distance against a reference), and reliability (convergence
diagnostics).

``` r

# Compare NUTS, VB, and Laplace on M2
vb_m2 <- VariationalBayes(Model_m2, IV_m2, Data_m2,
                           Covar = NULL, Iterations = 5000,
                           Method = "ADVI.mf")
la_m2 <- LaplaceApproximation(Model_m2, IV_m2, Data_m2,
                              Method = "SPG", Iterations = 500)

arena_m2 <- Arena(
    x = list(NUTS = fit_m2, VB = vb_m2, Laplace = la_m2),
    Model = Model_m2, Data = Data_m2
)
print(arena_m2)
plot(arena_m2)
```

The Arena plot shows each method’s position on the efficiency-accuracy
tradeoff. NUTS typically achieves the highest accuracy but lowest speed;
VB is fast but may underestimate posterior variance; Laplace is fastest
but limited to unimodal posteriors. The relative positions depend on the
model and data, which is why the comparison should be problem-specific
rather than based on general recommendations.

This cross-method comparison addresses the concern raised by Gelman et
al. [\[1\]](#ref1) about approximate algorithms: “there is no
one-size-fits-all approximate inference algorithm.” A variational
approximation that is accurate for a simple regression may be
dangerously inaccurate for a hierarchical model with funnel geometry. By
running both exact (NUTS) and approximate (VB, Laplace) methods on the
same problem and comparing the results, the analyst can assess whether
the computational shortcut is acceptable for the specific model at hand.
If the VB and NUTS posteriors agree closely, VB can be used for rapid
model exploration; if they disagree, NUTS is necessary for the final
inference.

## Effects, predictions, and target contexts

The final stage of the workflow extracts substantive conclusions from
the best-fitting model (M3). This is where the analytical pipeline
reaches the rightmost nodes in the workflow diagram: Estimates, Effects,
Predictions, and Target Contexts. The distinction between “estimates”
and “effects” is important. Estimates are posterior summaries of model
parameters; effects are scientifically interpretable quantities derived
from those parameters, often involving contrasts, transformations, or
propagation through the model’s deterministic structure. The Bayesian
framework handles this distinction naturally: any function of the
posterior samples is itself a posterior, with full uncertainty
quantification inherited from the joint posterior.

### Posterior summaries

The posterior distribution of M3 provides estimates for all parameters,
including the population-level kinetic parameters, the between-group
variability, the residual scale, and the degrees of freedom.

``` r

summary(fit_m3)
```

Credible intervals can be computed as highest posterior density (HPD)
intervals or as equal-tailed quantile intervals. HPD intervals are the
shortest intervals containing a given probability mass; quantile
intervals are symmetric in probability but not necessarily in width.

``` r

p.interval(fit_m3, HPD = TRUE, prob = 0.95)
p.interval(fit_m3, HPD = FALSE, prob = 0.95)
```

The caterpillar plot provides a visual summary of all parameter
estimates with their 50% and 95% credible intervals. The `ground_truth`
argument overlays the true parameter values for validation.

``` r

gt_m3 <- c(Vmax_pop = Vmax_pop, Km_pop = Km_pop,
           tau_Vmax = tau_Vmax, tau_Km = tau_Km,
           sigma = sigma_true, nu = nu_true,
           setNames(Vmax_true, paste0("Vmax_", 1:G)),
           setNames(Km_true, paste0("Km_", 1:G)))

caterpillar.plot(fit_m3, ground_truth = gt_m3)
```

### Derived quantities and effects

The posterior samples can be used to compute derived quantities with
full uncertainty propagation. For instance, the effective dose at 50% of
maximum response (\\\text{ED}\_{50}\\) is simply \\K_m\\ for the
Michaelis-Menten model, but for more complex dose-response models, it
would need to be computed from the posterior samples.

``` r

posterior <- fit_m3$Posterior2
Vmax_pop_post <- posterior[, "Vmax_pop"]
Km_pop_post   <- posterior[, "Km_pop"]

# ED50 = Km (by definition for Michaelis-Menten)
ED50_post <- Km_pop_post

# Population-level response at dose = 5
response_at_5 <- Vmax_pop_post * 5 / (Km_pop_post + 5)

# Between-group coefficient of variation for Vmax
cv_Vmax <- posterior[, "tau_Vmax"] / posterior[, "Vmax_pop"]

cat("ED50: ", sprintf("%.2f (%.2f, %.2f)",
    median(ED50_post),
    quantile(ED50_post, 0.025),
    quantile(ED50_post, 0.975)), "\n")

cat("Response at dose 5: ", sprintf("%.2f (%.2f, %.2f)",
    median(response_at_5),
    quantile(response_at_5, 0.025),
    quantile(response_at_5, 0.975)), "\n")

cat("CV(Vmax): ", sprintf("%.1f%% (%.1f%%, %.1f%%)",
    100 * median(cv_Vmax),
    100 * quantile(cv_Vmax, 0.025),
    100 * quantile(cv_Vmax, 0.975)), "\n")
```

### Predictions for observed groups

Predictions at new doses for observed groups use the group-specific
posterior estimates.

``` r

dose_new <- seq(0, 15, length.out = 100)
n_post <- nrow(posterior)

pred_list <- list()
for (g in seq_len(G)) {
    Vmax_g <- posterior[, paste0("Vmax_", g)]
    Km_g   <- posterior[, paste0("Km_", g)]
    sigma_g <- posterior[, "sigma"]
    nu_g   <- posterior[, "nu"]

    # Posterior predictive mean and intervals
    mu_mat <- outer(Vmax_g, dose_new, function(v, d)
        v * d / (Km_g + d))
    pred_list[[g]] <- data.frame(
        dose = dose_new,
        group = factor(g),
        med  = apply(mu_mat, 2, median),
        lo   = apply(mu_mat, 2, quantile, 0.025),
        hi   = apply(mu_mat, 2, quantile, 0.975)
    )
}
df_pred <- do.call(rbind, pred_list)

ggplot() +
    geom_point(data = df_syn, aes(x = dose, y = y, color = group),
               alpha = 0.5, size = 1.5) +
    geom_ribbon(data = df_pred, aes(x = dose, ymin = lo, ymax = hi,
                fill = group), alpha = 0.15) +
    geom_line(data = df_pred, aes(x = dose, y = med, color = group),
              linewidth = 0.8) +
    scale_color_brewer(palette = "Set2", name = "Group") +
    scale_fill_brewer(palette = "Set2", name = "Group") +
    labs(x = "Dose", y = "Response",
         title = "Posterior predictive dose-response curves",
         subtitle = paste0("Hierarchical Student-t model (M3), ",
                           "median and 95% credible band"),
         caption = paste0("Predictions at 100 dose values for each ",
                          "observed group. Bands reflect both parameter ",
                          "uncertainty and residual variability.")) +
    theme_minimal() +
    theme(plot.caption = element_text(color = "grey40", size = 8, hjust = 0))
```

### Predictions for a new group

Hierarchical prediction for an unobserved group uses the
population-level distribution. New group parameters are drawn from the
population prior, incorporating both the uncertainty in the
population-level parameters and the between-group variability.

``` r

Vmax_new <- rnorm(n_post, posterior[, "Vmax_pop"], posterior[, "tau_Vmax"])
Km_new   <- pmax(rnorm(n_post, posterior[, "Km_pop"],
                       posterior[, "tau_Km"]), 0.01)

mu_new <- outer(Vmax_new, dose_new, function(v, d) v * d / (Km_new + d))
df_new <- data.frame(
    dose = dose_new,
    med  = apply(mu_new, 2, median),
    lo   = apply(mu_new, 2, quantile, 0.025),
    hi   = apply(mu_new, 2, quantile, 0.975)
)

ggplot(df_new, aes(x = dose)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.3) +
    geom_line(aes(y = med), color = "#2166AC", linewidth = 0.8) +
    labs(x = "Dose", y = "Response",
         title = "Prediction for a new (unobserved) group",
         subtitle = "Population-level hierarchical prediction with 95% band",
         caption = paste0("Wider bands than group-specific predictions ",
                          "because new-group uncertainty includes both ",
                          "population parameter uncertainty and between-group ",
                          "variability.")) +
    theme_minimal() +
    theme(plot.caption = element_text(color = "grey40", size = 8, hjust = 0))
```

The prediction band for a new group is wider than the group-specific
bands because it incorporates both the uncertainty in the
population-level parameters and the between-group variability. This is
the hallmark of hierarchical prediction: the model honestly reports that
predicting a new, unobserved group is more uncertain than predicting an
observed one.

This honest reporting of uncertainty is one of the strongest practical
arguments for the Bayesian workflow. A frequentist prediction interval
can account for sampling variability and parameter uncertainty, but it
cannot naturally account for between-group variability in a hierarchical
setting without bootstrapping or other simulation-based approaches. The
Bayesian hierarchical model provides this naturally: the
population-level posterior encodes the distribution of plausible
group-level parameters, and sampling from this distribution generates
predictions with appropriate uncertainty. The result is that a
policy-maker who needs to know “what would happen in a new, untested
population” receives an interval that reflects all relevant sources of
uncertainty, not just the sampling variability of the observed groups.

## Summary: workflow nodes and lucifer functions

The table below maps every node and arrow in the workflow diagram to the
lucifer functions exercised in this vignette.

``` r

df_map <- data.frame(
    Node = c("Theory", "Estimand", "Design", "Target Contexts",
             "Sample", "Synthetic sample", "Estimator",
             "Prior predictive", "Synthetic estimates",
             "Post predictive", "Estimates", "Diagnostics",
             "Data cloning", "Effects", "Predictions"),
    Stream = c(rep("Black", 5), "Red", "Black", rep("Red", 3),
               "Black", "Red", "Red", "Black", "Black"),
    Functions = c(
        "Narrative; scientific model",
        "mon.names, parm.names in Data list",
        "Data list structure; CenterScale()",
        "predict() with new Data list",
        "Data$y, Data$dose, Data$group",
        "set.seed(); simulation from DGP",
        "model_spec(), Prescribe(), lucifer(), VB(), LA(), Crucible()",
        "prior_predictive_check() (4 types)",
        "caterpillar.plot(ground_truth=), coverage_sim()",
        "predict(), plot.demonoid.ppc() (40 styles)",
        "summary(), p.interval(), caterpillar.plot()",
        "Consort(), Rhat(), ESS(), RobustBayes(), LOO(), check_nuts()",
        "DataCloning(); identifiability, estimability, MLE bridge",
        "Posterior contrasts; derived quantities from Posterior2",
        "predict() with new doses; hierarchical prediction"
    ),
    Section = c("1", "1", "1", "1,11", "2", "2", "3,5,6",
                "4", "5", "7", "5,11", "6,8,10", "9", "11", "11"),
    stringsAsFactors = FALSE
)

knitr::kable(df_map, col.names = c("Diagram node", "Stream",
    "lucifer functions", "Sections"),
    caption = "Mapping between workflow diagram nodes and lucifer functions.")
```

The arrows in the diagram correspond to the transitions between
sections. The three red-loop iterations are:

1.  **Prior predictive loop (Section 4):** M0 \\\to\\
    [`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
    \\\to\\ fails (no saturation) \\\to\\ M1 \\\to\\
    [`prior_predictive_check()`](https://robustecologies.github.io/lucifer/reference/prior_predictive_check.md)
    \\\to\\ passes.

2.  **Posterior predictive loop (Section 7):** M1 \\\to\\
    [`predict()`](https://rdrr.io/r/stats/predict.html) \\\to\\
    `plot(Style = "Violin Grouped")` \\\to\\ group heterogeneity
    detected \\\to\\ expand to M2.

3.  **Sensitivity and PPC loop (Section 8):** M2 \\\to\\
    `plot(Style = "Stat")` \\\to\\ heavy tails detected \\\to\\
    `RobustBayes(modules = "influence")` \\\to\\ high-leverage
    observations confirmed \\\to\\ expand to M3 \\\to\\ all checks pass.

The workflow is iterative, not linear. Each pass through the diagnostic
loop refines the model, and the comparison tools in Section 9 quantify
the improvement. The final model M3 is not the endpoint of the workflow
but a provisional best model that future data or domain knowledge may
motivate revising. The function
[`is.appeased()`](https://robustecologies.github.io/lucifer/reference/is.appeased.md)
returning `TRUE` means the computation has converged and the model
passes its own internal consistency checks. It does not mean the model
is correct.

### Reflections on the workflow

The progression from M0 through M3 illustrates a general pattern that
Gelman et al. [\[1\]](#ref1) describe as “fitting many models on the way
to the final one.” The abandoned models M0 and M1 are not failures but
necessary scaffolding that builds understanding. M0 was discarded before
fitting because its functional form contradicted scientific knowledge
(no saturation), a decision enabled by prior predictive checking. M1 was
fit successfully but revealed group heterogeneity through posterior
predictive checking, motivating the hierarchical expansion to M2. M2 was
fit successfully and captured the group structure but failed to
accommodate heavy-tailed residuals, motivating the robust expansion to
M3. At each transition, the diagnostic that triggered the expansion was
a specific, identifiable discrepancy between the model’s predictions and
the data or domain knowledge, not a vague sense that the model could be
“better.”

This pattern, where each iteration addresses a specific identified
failure rather than a general aspiration to improve, is the hallmark of
a disciplined workflow. Box’s dictum that “all models are wrong, but
some are useful” [\[7\]](#ref7) is sometimes invoked to justify fitting
a single model without checking. The workflow perspective inverts this:
precisely because all models are wrong, we need systematic tools to
discover how they are wrong and whether the specific ways in which they
are wrong matter for the estimand at hand. A model that is “wrong” in
ways that do not affect the target quantity is useful; a model that is
wrong in ways that bias the estimand is not, regardless of its formal
complexity or fit statistics.

The Bayesian approach to this iterative refinement has two structural
advantages over its frequentist counterpart. First, the generative model
provides a natural simulation engine at every stage: prior predictive
checks simulate from \\p(\tilde{y})\\ before seeing data, posterior
predictive checks simulate from \\p(\tilde{y} \mid y)\\ after fitting,
and sensitivity analysis perturbs \\p(\theta)\\ or \\p(y \mid \theta)\\
to assess robustness. These simulations are not ad-hoc diagnostics
bolted onto the inference but natural consequences of the probabilistic
model itself. Second, the explicit treatment of uncertainty, propagated
through every derived quantity and prediction, means that the analyst is
always confronted with what the model does not know, not just what it
estimates. A frequentist confidence interval covers the true parameter
with specified probability under repeated sampling from the same design;
a Bayesian credible interval quantifies the remaining uncertainty
conditional on the data actually observed. The workflow perspective
makes this distinction operational: every prediction in Section 10
carries uncertainty bands that reflect both parameter uncertainty and
between-group variability, and the hierarchical prediction for a new
group is honestly wider than the prediction for an observed group.

The workflow also provides a natural framework for what Gelman et
al. [\[1\]](#ref1) call “modeling as software development.” Just as
software is never finished but only released in versions that pass their
test suites, a Bayesian model is never final but only provisionally
accepted after passing its diagnostic suite. The difference between a
model that has been checked and one that has not is not that the checked
model is necessarily correct, but that the ways in which it might be
wrong have been enumerated, investigated, and either resolved or
documented. This disciplined honesty, more than any single algorithm or
diagnostic, is the contribution of the Bayesian workflow to applied
statistics.

## References

**\[1\]** Gelman, A., Vehtari, A., Simpson, D., Margossian, C.C.,
Carpenter, B., Yao, Y., Kennedy, L., Gabry, J., Burkner, P.-C., and
Modrak, M. (2020). *Bayesian workflow*. arXiv:2011.01808.

**\[2\]** Talts, S., Betancourt, M., Simpson, D., Vehtari, A., and
Gelman, A. (2018). *Validating Bayesian inference algorithms with
simulation-based calibration*. arXiv:1804.06788.

**\[3\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and
Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an
improved \\\hat{R}\\ for assessing convergence of MCMC (with
discussion). *Bayesian Analysis*, 16(2), 667-718.

**\[4\]** Kallioinen, N., Paananen, T., Burkner, P.-C., and Vehtari, A.
(2024). Detecting and diagnosing prior and likelihood sensitivity with
power-scaling. *Statistics and Computing*, 34(57).

**\[5\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical
Bayesian model evaluation using leave-one-out cross-validation and WAIC.
*Statistics and Computing*, 27(5), 1413-1432.

**\[6\]** Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and
Gelman, A. (2019). Visualization in Bayesian workflow. *Journal of the
Royal Statistical Society: Series A*, 182(2), 389-402.

**\[7\]** Box, G.E.P. (1980). Sampling and Bayes’ inference in
scientific modelling and robustness (with discussion). *Journal of the
Royal Statistical Society: Series A*, 143(4), 383-430.

**\[8\]** Rubin, D.B. (1984). Bayesianly justifiable and relevant
frequency calculations for the applied statistician. *Annals of
Statistics*, 12(4), 1151-1172.

**\[9\]** Berger, J.O. (1994). An overview of robust Bayesian analysis
(with discussion). *Test*, 3(1), 5-124.

**\[10\]** Betancourt, M. (2020). *Towards a principled Bayesian
workflow*.
betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html.

**\[11\]** Simpson, D., Rue, H., Riebler, A., Martins, T.G., and Sorbye,
S.H. (2017). Penalising model component complexity: a principled,
practical approach to constructing priors. *Statistical Science*, 32(1),
1-28.

**\[12\]** Lele, S.R., Dennis, B. and Lutscher, F. (2007). Data cloning:
easy maximum likelihood estimation for complex ecological models using
Bayesian Markov chain Monte Carlo methods. *Ecology Letters*, 10,
551-563.

**\[13\]** Campbell, D. and Lele, S. (2014). An ANOVA test for parameter
estimability using data cloning with application to statistical
inference for dynamic systems. *Computational Statistics and Data
Analysis*, 70, 257-267.
