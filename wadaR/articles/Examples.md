# A bestiary of models to appaise lucifer

  

## Overview

The **lucifer** package provides a comprehensive environment for
Bayesian inference, supporting over 82 MCMC algorithms, variational
Bayes, Laplace approximation, sequential Monte Carlo, and
simulation-based inference. This vignette collects more than 100 worked
examples spanning the full range of models that lucifer can handle, from
simple linear regression to state-space models, spatial processes, and
mixture models.

Each example is self-contained and includes: a brief description of the
model and its purpose; the mathematical form expressed in standard
probabilistic notation; simulated ground truth data with known parameter
values; the `Model` function; initial values; and a fitting-and-recovery
section that uses
[`Prescribe()`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
to select an appropriate algorithm, fits the model with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and evaluates parameter recovery through posterior predictive checks.
For models expressible within the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
declarative DSL, the equivalent specification is also shown alongside
the manual `Model` function so that users can compare both approaches.

All code blocks in this vignette are shown but not evaluated during the
build; they are designed to be copied into an R session and run
interactively. Initial values are usually hard-coded, though the
parameter-generating function (PGF) is also specified. It is recommended
to generate initial values with the `GIV` function according to the
user-specified PGF.

Notation follows these conventions: Greek letters represent parameters;
lower-case letters represent indices; lower-case boldface letters
represent scalars or vectors; probability distributions are represented
with calligraphic font (\\\mathcal{N}\\, \\\mathcal{HC}\\, etc.);
upper-case letters represent index limits; and upper-case boldface
letters represent matrices.

  

## Linear regression and extensions

### ANCOVA

This example is essentially the same as the two-way ANOVA (see section
[link](#anova.two.way)), except that a continuous covariate
\\\textbf{X}\_{,3}\\ has been added, and its parameter is \\\delta\\.
The model combines two crossed factors with sum-to-zero constraints and
a covariate effect, with hierarchical priors on the group-level variance
components. The analysis of covariance combines Fisher’s ANOVA framework
(Fisher [\[41\]](#ref41)) with regression adjustment, following the
Bayesian formulation in Gelman et al. [\[48\]](#ref48). ANCOVA is
routinely applied in clinical trials to adjust treatment comparisons for
baseline covariates such as pre-treatment severity scores.

#### Form

\\\textbf{y}\_i \sim \mathcal{N}(\mu_i, \sigma^2_1)\\ \\\mu_i = \alpha +
\beta\[\textbf{X}\_{i,1}\] + \gamma\[\textbf{X}\_{i,2}\] + \delta
\textbf{X}\_{i,3}, \quad i=1,\dots,N\\ \\\epsilon_i = \textbf{y}\_i -
\mu_i\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\beta_j \sim
\mathcal{N}(0, \sigma^2_2), \quad j=1,\dots,J\\ \\\beta_J = -
\sum^{J-1}\_{j=1} \beta_j\\ \\\gamma_k \sim \mathcal{N}(0, \sigma^2_3),
\quad k=1,\dots,K\\ \\\gamma_K = - \sum^{K-1}\_{k=1} \gamma_k\\ \\\delta
\sim \mathcal{N}(0, 1000)\\ \\\sigma_m \sim \mathcal{HC}(25), \quad
m=1,\dots,3\\

#### model_spec() notation

The ANCOVA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires sum-to-zero constraints on the factor effects
(\\\beta_J = -\sum\_{j=1}^{J-1} \beta_j\\ and \\\gamma_K =
-\sum\_{k=1}^{K-1} \gamma_k\\). These constraints involve deriving the
last level from the sampled parameters through an imperative operation
on a parameter subset, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42102)
N <- 200
J <- 5  # Number of levels in factor 1
K <- 3  # Number of levels in factor 2

#### True parameter values
true.alpha <- 1.5
true.beta <- c(-1.2, 0.8, 0.3, -0.5)      # J-1 free effects
true.beta <- c(true.beta, -sum(true.beta))  # beta[J] from sum-to-zero
true.gamma <- c(0.9, -0.4)                  # K-1 free effects
true.gamma <- c(true.gamma, -sum(true.gamma)) # gamma[K] from sum-to-zero
true.delta <- 0.75
true.sigma <- c(0.5, 1.0, 0.8) # sigma[1]=obs, sigma[2]=beta sd, sigma[3]=gamma sd

#### Generate design matrix
X <- cbind(
    rcat(N, rep(1/J, J)),  # Factor 1 (J levels)
    rcat(N, rep(1/K, K)),  # Factor 2 (K levels)
    runif(N, -2, 2)        # Continuous covariate
)

#### Generate response from the true model
mu.true <- true.alpha + true.beta[X[,1]] + true.gamma[X[,2]] +
    true.delta * X[,3]
y <- mu.true + rnorm(N, 0, true.sigma[1])

#### Assemble Data list
mon.names <- c("LP", paste0("beta[", J, "]"),
    paste0("gamma[", K, "]"), "s.beta", "s.gamma", "s.epsilon")
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,J-1),
    gamma=rep(0,K-1), delta=0, sigma=rep(0,3)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(Data$J-1)
    gamma <- rnorm(Data$K-1)
    delta <- rnorm(1)
    sigma <- runif(3)
    return(c(alpha, beta, gamma, delta, sigma))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha,
    pos.beta=pos.beta, pos.gamma=pos.gamma, pos.delta=pos.delta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    beta <- c(beta, -sum(beta)) #Sum-to-zero constraint
    gamma <- parm[Data$pos.gamma]
    gamma <- c(gamma, -sum(gamma)) #Sum-to-zero constraint
    delta <- parm[Data$pos.delta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnorm(beta, 0, sigma[2], log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma[3], log=TRUE))
    delta.prior <- dnormv(delta, 0, 1000, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- alpha + beta[Data$X[,1]] + gamma[Data$X[,2]] +
        delta*Data$X[,3]
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Variance Components
    s.beta <- sd(beta)
    s.gamma <- sd(gamma)
    s.epsilon <- sd(Data$y - mu)
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + delta.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, beta[Data$J],
        gamma[Data$K], s.beta, s.gamma, s.epsilon),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,(J-1)), rep(0,(K-1)), 0, rep(1,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: intercept
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")

#### Parameter recovery: factor 1 effects (J-1 free + derived J-th)
post.beta.free <- fit$Posterior2[, pos.beta]
post.betaJ <- -rowSums(post.beta.free)
for (j in 1:(J-1))
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(mean(post.beta.free[, j]), 3), "\n")
cat("beta[", J, "] -- true:", true.beta[J],
    " post. mean:", round(mean(post.betaJ), 3), " (derived)\n")

#### Parameter recovery: factor 2 effects (K-1 free + derived K-th)
post.gamma.free <- fit$Posterior2[, pos.gamma]
post.gammaK <- -rowSums(post.gamma.free)
for (k in 1:(K-1))
    cat("gamma[", k, "] -- true:", true.gamma[k],
        " post. mean:", round(mean(post.gamma.free[, k]), 3), "\n")
cat("gamma[", K, "] -- true:", true.gamma[K],
    " post. mean:", round(mean(post.gammaK), 3), " (derived)\n")

#### Parameter recovery: covariate effect
cat("delta -- true:", true.delta,
    " post. mean:", round(fit$Summary2[pos.delta, "Mean"], 3), "\n")

#### Parameter recovery: scale parameters
for (m in 1:3)
    cat("sigma[", m, "] -- true:", true.sigma[m],
        " post. mean:", round(fit$Summary2[pos.sigma[m], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### ANOVA, One-Way

When \\J=2\\, this is a Bayesian form of a t-test. The group effects
\\\beta_j\\ are hierarchically centred around zero with a shared
standard deviation \\\sigma_2\\ estimated from the data, and a
sum-to-zero constraint identifies the overall intercept \\\alpha\\
separately from the group contrasts. The Bayesian one-way ANOVA with
sum-to-zero constraints follows the parameterization discussed in Gelman
et al. [\[48\]](#ref48). One-way ANOVA is a workhorse in agricultural
science for comparing crop yields across different fertilizer
treatments.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2_1)\\ \\\mu_i = \alpha +
\beta\[\textbf{x}\_i\], \quad i=1,\dots,N\\ \\\alpha \sim \mathcal{N}(0,
1000)\\ \\\beta_j \sim \mathcal{N}(0, \sigma^2_2), \quad j=1,\dots,J\\
\\\beta_J = - \displaystyle\sum^{J-1}\_{j=1} \beta_j\\ \\\sigma\_{1:2}
\sim \mathcal{HC}(25)\\

#### model_spec() notation

The one-way ANOVA cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a sum-to-zero constraint on the group effects
(\\\beta_J = -\sum\_{j=1}^{J-1} \beta_j\\). This constraint derives the
last level from the sampled parameters through an imperative operation
that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42103)
N <- 200
J <- 4

#### True parameter values
true.alpha <- 2.0
true.beta <- c(-1.5, 0.5, 0.3)            # J-1 free effects
true.beta <- c(true.beta, -sum(true.beta)) # beta[J] from sum-to-zero = 0.7
true.sigma <- c(0.8, 1.0)                  # sigma[1]=obs, sigma[2]=beta sd

#### Generate design
x <- rcat(N, rep(1/J, J))

#### Generate response from the true model
mu.true <- true.alpha + true.beta[x]
y <- mu.true + rnorm(N, 0, true.sigma[1])

#### Assemble Data list
mon.names <- c("LP", paste0("beta[", J, "]"))
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,J-1), sigma=rep(0,2)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(Data$J-1)
    sigma <- runif(2)
    return(c(alpha, beta, sigma))
    }
Data <- list(J=J, N=N, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.sigma=pos.sigma, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    beta <- c(beta, -sum(beta)) #Sum-to-zero constraint
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnorm(beta, 0, sigma[2], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- alpha + beta[Data$x]
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,beta[Data$J]),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,(J-1)), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery: intercept
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")

#### Parameter recovery: group effects (J-1 free + derived J-th)
post.beta.free <- fit$Posterior2[, pos.beta]
post.betaJ <- -rowSums(post.beta.free)
for (j in 1:(J-1))
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(mean(post.beta.free[, j]), 3), "\n")
cat("beta[", J, "] -- true:", true.beta[J],
    " post. mean:", round(mean(post.betaJ), 3), " (derived)\n")

#### Parameter recovery: scale parameters
for (m in 1:2)
    cat("sigma[", m, "] -- true:", true.sigma[m],
        " post. mean:", round(fit$Summary2[pos.sigma[m], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### ANOVA, Two-Way

In this representation, \\\sigma^m\\ are the superpopulation variance
components, `s.beta` and `s.gamma` are the finite-population
within-variance components of the factors or treatments, and `s.epsilon`
is the finite-population between-variance component. The two-way ANOVA
extends Fisher’s original design framework (Fisher, 1925) to crossed
factors, with the Bayesian formulation following Gelman et
al. [\[48\]](#ref48). Two-way factorial designs are common in
manufacturing for simultaneously evaluating the effects of temperature
and pressure on product quality.

#### Form

\\\textbf{y}\_i \sim \mathcal{N}(\mu_i, \sigma^2_1)\\ \\\mu_i = \alpha +
\beta\[\textbf{X}\_{i,1}\] + \gamma\[\textbf{X}\_{i,2}\], \quad
i=1,\dots,N\\ \\\epsilon_i = \textbf{y}\_i - \mu_i\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\beta_j \sim \mathcal{N}(0, \sigma^2_2), \quad
j=1,\dots,J\\ \\\beta_J = - \sum^{J-1}\_{j=1} \beta_j\\ \\\gamma_k \sim
\mathcal{N}(0, \sigma^2_3), \quad k=1,\dots,K\\ \\\gamma_K = -
\sum^{K-1}\_{k=1} \gamma_k\\ \\\sigma_m \sim \mathcal{HC}(25), \quad
m=1,\dots,3\\

#### model_spec() notation

The two-way ANOVA cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires sum-to-zero constraints on both factor effects
(\\\beta_J = -\sum\_{j=1}^{J-1} \beta_j\\ and \\\gamma_K =
-\sum\_{k=1}^{K-1} \gamma_k\\). These constraints derive the last level
of each factor from the sampled parameters through imperative operations
that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42104)
N <- 300
J <- 5  # Number of levels in factor (treatment) 1
K <- 3  # Number of levels in factor (treatment) 2

#### True parameter values
true.alpha <- 1.8
true.beta <- c(-1.0, 0.6, 0.2, -0.3)       # J-1 free effects
true.beta <- c(true.beta, -sum(true.beta))   # beta[J] from sum-to-zero = 0.5
true.gamma <- c(0.7, -0.4)                   # K-1 free effects
true.gamma <- c(true.gamma, -sum(true.gamma)) # gamma[K] from sum-to-zero = -0.3
true.sigma <- c(0.6, 0.9, 0.7)  # sigma[1]=obs, sigma[2]=beta sd, sigma[3]=gamma sd

#### Generate design matrix
X <- cbind(
    rcat(N, rep(1/J, J)),  # Factor 1 (J levels)
    rcat(N, rep(1/K, K))   # Factor 2 (K levels)
)

#### Generate response from the true model
mu.true <- true.alpha + true.beta[X[,1]] + true.gamma[X[,2]]
y <- mu.true + rnorm(N, 0, true.sigma[1])

#### Assemble Data list
mon.names <- c("LP", paste0("beta[", J, "]"),
    paste0("gamma[", K, "]"), "s.beta", "s.gamma", "s.epsilon")
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,J-1), gamma=rep(0,K-1),
    sigma=rep(0,3)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(Data$J-1)
    gamma <- rnorm(Data$K-1)
    sigma <- runif(3)
    return(c(alpha, beta, gamma, sigma))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.gamma=pos.gamma, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    beta <- c(beta, -sum(beta)) #Sum-to-zero constraint
    gamma <- parm[Data$pos.gamma]
    gamma <- c(gamma, -sum(gamma)) #Sum-to-zero constraint
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnorm(beta, 0, sigma[2], log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma[3], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- alpha + beta[Data$X[,1]] + gamma[Data$X[,2]]
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Variance Components
    s.beta <- sd(beta)
    s.gamma <- sd(gamma)
    s.epsilon <- sd(Data$y - mu)
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, beta[Data$J],
        gamma[Data$K], s.beta, s.gamma, s.epsilon),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,(J-1)), rep(0,(K-1)), rep(1,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery: intercept
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")

#### Parameter recovery: factor 1 effects (J-1 free + derived J-th)
post.beta.free <- fit$Posterior2[, pos.beta]
post.betaJ <- -rowSums(post.beta.free)
for (j in 1:(J-1))
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(mean(post.beta.free[, j]), 3), "\n")
cat("beta[", J, "] -- true:", true.beta[J],
    " post. mean:", round(mean(post.betaJ), 3), " (derived)\n")

#### Parameter recovery: factor 2 effects (K-1 free + derived K-th)
post.gamma.free <- fit$Posterior2[, pos.gamma]
post.gammaK <- -rowSums(post.gamma.free)
for (k in 1:(K-1))
    cat("gamma[", k, "] -- true:", true.gamma[k],
        " post. mean:", round(mean(post.gamma.free[, k]), 3), "\n")
cat("gamma[", K, "] -- true:", true.gamma[K],
    " post. mean:", round(mean(post.gammaK), 3), " (derived)\n")

#### Parameter recovery: scale parameters
for (m in 1:3)
    cat("sigma[", m, "] -- true:", true.sigma[m],
        " post. mean:", round(fit$Summary2[pos.sigma[m], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Errors-in-Variables Regression **\[NEW\]**

Errors-in-variables regression corrects the attenuation bias that arises
when covariates are measured with noise by introducing latent true
values \\x^\star_i\\ as parameters, jointly estimating the regression
coefficients and the unobserved predictors from both the noisy covariate
measurements and the response. The measurement error variance
\\\sigma^2_x\\ is assumed known, while the latent \\x^\star_i\\ receive
a hierarchical normal prior centred at \\\mu_x\\ with spread \\\tau_x\\,
which regularises the estimation and provides partial pooling across
observations. The errors-in-variables framework was formalised by Fuller
[\[126\]](#ref126), and the Bayesian treatment with measurement error
models was developed comprehensively by Carroll et
al. [\[127\]](#ref127). Errors-in-variables regression is critical in
nutritional epidemiology, where dietary intake measured by food
frequency questionnaires contains substantial measurement error that
biases diet-disease associations toward the null.

#### Form

\\y_i \sim \mathcal{N}(\alpha + \beta x^\star_i, \sigma^2_y), \quad
i=1,\dots,N\\ \\x^{\text{obs}}\_i \sim \mathcal{N}(x^\star_i,
\sigma^2_x), \quad \sigma_x \text{ known}\\ \\x^\star_i \sim
\mathcal{N}(\mu_x, \tau^2_x), \quad i=1,\dots,N\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\beta \sim \mathcal{N}(0, 1000)\\ \\\sigma_y
\sim \mathcal{HC}(25)\\ \\\mu_x \sim \mathcal{N}(0, 1000)\\ \\\tau_x
\sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the latent true covariates \\x^\star_i\\ are treated as
parameters with their own hierarchical prior and an auxiliary
measurement model for \\x^{\text{obs}}\_i\\, requiring a joint
likelihood over both the response and the observed covariates that falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42211)
N <- 100

#### True parameter values
true.alpha <- 1.0
true.beta <- 2.0
true.sigma.y <- 0.5
true.mu.x <- 3.0
true.tau.x <- 1.0
true.sigma.x <- 0.3   # Known measurement error SD

#### Simulate latent and observed data
x.star.true <- rnorm(N, true.mu.x, true.tau.x)
x.obs <- rnorm(N, x.star.true, true.sigma.x)
y <- rnorm(N, true.alpha + true.beta * x.star.true, true.sigma.y)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(
    alpha=0, beta=0, sigma.y=0, mu.x=0, tau.x=0, x.star=rep(0,N)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma.y <- grep("sigma.y", parm.names)
pos.mu.x <- grep("mu.x", parm.names)
pos.tau.x <- grep("tau.x", parm.names)
pos.x.star <- grep("x.star", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(1)
    sigma.y <- runif(1)
    mu.x <- rnorm(1, mean(Data$x.obs), 1)
    tau.x <- runif(1)
    x.star <- Data$x.obs + rnorm(Data$N, 0, 0.1)
    return(c(alpha, beta, sigma.y, mu.x, tau.x, x.star))
    }
Data <- list(N=N, x.obs=x.obs, y=y, sigma.x=true.sigma.x,
    PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.beta=pos.beta, pos.sigma.y=pos.sigma.y,
    pos.mu.x=pos.mu.x, pos.tau.x=pos.tau.x, pos.x.star=pos.x.star)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    sigma.y <- interval(parm[Data$pos.sigma.y], 1e-100, Inf)
    parm[Data$pos.sigma.y] <- sigma.y
    mu.x <- parm[Data$pos.mu.x]
    tau.x <- interval(parm[Data$pos.tau.x], 1e-100, Inf)
    parm[Data$pos.tau.x] <- tau.x
    x.star <- parm[Data$pos.x.star]
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- dnormv(beta, 0, 1000, log=TRUE)
    sigma.y.prior <- dhalfcauchy(sigma.y, 25, log=TRUE)
    mu.x.prior <- dnormv(mu.x, 0, 1000, log=TRUE)
    tau.x.prior <- dhalfcauchy(tau.x, 25, log=TRUE)
    x.star.prior <- sum(dnorm(x.star, mu.x, tau.x, log=TRUE))
    ### Log-Likelihood
    LL.y <- sum(dnorm(Data$y, alpha + beta * x.star, sigma.y, log=TRUE))
    LL.x <- sum(dnorm(Data$x.obs, x.star, Data$sigma.x, log=TRUE))
    LL <- LL.y + LL.x
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.y.prior +
        mu.x.prior + tau.x.prior + x.star.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, alpha + beta * x.star, sigma.y), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 1, 0.5, mean(x.obs), 1, x.obs)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha   -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("beta    -- true:", true.beta,
    " post. mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("sigma.y -- true:", true.sigma.y,
    " post. mean:", round(fit$Summary2[pos.sigma.y, "Mean"], 3), "\n")
cat("mu.x    -- true:", true.mu.x,
    " post. mean:", round(fit$Summary2[pos.mu.x, "Mean"], 3), "\n")
cat("tau.x   -- true:", true.tau.x,
    " post. mean:", round(fit$Summary2[pos.tau.x, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Generalized Least Squares (GLS) **\[NEW\]**

Generalized least squares extends ordinary linear regression by
modelling correlated residuals through a structured covariance matrix,
here an AR(1) correlation \\R\_{ij} = \rho^{\|i-j\|}\\ scaled by the
error variance \\\sigma^2\\, so that the response vector follows a
multivariate normal distribution with non-diagonal covariance.
Estimating \\\rho\\ alongside the regression coefficients accounts for
the serial dependence that would otherwise inflate standard errors and
distort inference in time-ordered or spatially-ordered data. The GLS
estimator was introduced by Aitken [\[151\]](#ref151) as the
minimum-variance linear unbiased estimator under a known covariance
structure. GLS with AR(1) errors is used in econometrics to correct for
autocorrelation when modelling macroeconomic time series such as GDP
growth on lagged predictors.

#### Form

\\\textbf{y} \sim \mathcal{N}\_N(\textbf{X}\beta, \sigma^2 \textbf{R})\\
\\R\_{ij} = \rho^{\|i-j\|}, \quad i,j=1,\dots,N\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\ \\\rho \sim \mathcal{U}(0, 1)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the likelihood is a multivariate normal with a structured
AR(1) correlation matrix whose off-diagonal entries depend on the
estimated parameter \\\rho\\, requiring construction of the full \\N
\times N\\ covariance matrix inside the model function. The manual
`Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42214)
N <- 50
J <- 2   # Intercept + 1 covariate

#### True parameter values
true.beta <- c(2.0, 3.0)
true.sigma <- 1.0
true.rho <- 0.7

#### Generate design matrix
x1 <- sort(runif(N))
X <- cbind(1, x1)

#### Generate correlated response
R.true <- true.rho^abs(outer(1:N, 1:N, "-"))
V.true <- true.sigma^2 * R.true
L <- chol(V.true)
y <- as.vector(X %*% true.beta + crossprod(L, rnorm(N)))

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, rho=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.rho <- grep("rho", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    rho <- runif(1, 0.1, 0.9)
    return(c(beta, sigma, rho))
    }
Data <- list(N=N, J=J, X=X, y=y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma,
    pos.rho=pos.rho)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    rho <- interval(parm[Data$pos.rho], 1e-100, 1 - 1e-100)
    parm[Data$pos.rho] <- rho
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    rho.prior <- dunif(rho, 0, 1, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(tcrossprod(Data$X, t(beta)))
    R <- rho^abs(outer(1:Data$N, 1:Data$N, "-"))
    Sigma <- sigma^2 * R
    LL <- dmvn(Data$y, mu, Sigma, log=TRUE)
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior + rho.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=as.vector(rmvn(1, mu, Sigma)), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0, 1, 0.5)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta: ", true.beta, "\n")
cat("Post. mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True sigma:", true.sigma, "\n")
cat("Post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
cat("True rho:  ", true.rho, "\n")
cat("Post. mean:", round(fit$Summary2[pos.rho, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.rho]   <- true.rho[seq_along(pos.rho)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Laplace Regression

Laplace (double-exponential) regression replaces the normal likelihood
with a Laplace distribution, yielding the Bayesian analogue of median
regression; the heavier tails and sharp peak at the mode make it more
robust to outliers than standard normal regression while still admitting
a closed-form log-likelihood. It has been claimed that it should be
surprising that the normal distribution became the standard, when the
Laplace distribution usually fits better and has wider tails
[\[65\]](#ref65). Another popular alternative is to use the
t-distribution (see Robust Regression in section [link](#robust.reg)),
though it is more computationally expensive to estimate, because it has
three parameters. The Laplace distribution has only two parameters,
location and scale like the normal distribution, and is computationally
easier to fit. Laplace regression uses the asymmetric Laplace
distribution for quantile estimation, as developed by Yu and Moyeed
[\[116\]](#ref116) in a Bayesian framework. Bayesian quantile regression
is applied in economics to study how the determinants of wages differ
across the income distribution.

#### Form

\\\textbf{y} \sim \mathcal{L}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The Laplace distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry as `Laplace(location, scale)`. The linear predictor is
expressed through a deterministic node, producing a compact
specification equivalent to the manual Model function.

``` r

spec <- model_spec("
  y ~ Laplace(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42135)
N <- 300
J <- 4

#### True parameter values
true.beta <- c(1.5, -0.8, 0.6, 0.3)  # intercept + 3 predictors
true.sigma <- 0.5                      # scale parameter

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rlaplace(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dlaplace(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rlaplace(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression

Linear regression is the most fundamental Bayesian regression model and
the natural starting point for any probabilistic modelling workflow. The
likelihood is normal with mean given by a linear predictor \\\mu =
\textbf{X}\beta\\, vague normal priors are placed on the regression
coefficients \\\beta\\, and a half-Cauchy prior on the residual standard
deviation \\\sigma\\ provides a weakly informative constraint that keeps
the scale parameter positive while remaining diffuse enough to let the
data dominate inference. Bayesian linear regression dates to Lindley and
Smith [\[72\]](#ref72), extending the classical least squares framework
of Gauss [\[44\]](#ref44) with prior distributions on the coefficients.
Linear regression is the foundation of empirical science, applied
universally from predicting blood pressure from patient characteristics
to forecasting energy consumption.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

This model is fully compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
which compresses the Data list, Model function, parameter names, and
initial values into a single declaration.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

Rather than relying on a pre-packaged dataset, we simulate data from
known parameter values so that posterior recovery can be verified
against the ground truth. The intercept is 2.5 with three predictor
coefficients of \\-1.3\\, \\0.8\\, and \\0.4\\, and a residual standard
deviation of 1.2.

``` r

set.seed(42148)
N <- 200
J <- 4
true.beta <- c(2.5, -1.3, 0.8, 0.4)
true.sigma <- 1.2
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- rnorm(N, X %*% true.beta, true.sigma)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare estimates to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression with Power Priors

Power priors [\[59\]](#ref59) are a class of informative priors for
incorporating relevant historical data into the analysis of a current
dataset. Both the current data (\\\textbf{y}\\, \\\textbf{X}\\) and
historical data (\\\textbf{y}\_h\\, \\\textbf{X}\_h\\) contribute to the
posterior, but the historical likelihood is raised to a power \\\alpha
\in \[0,1\]\\ that controls how much influence the historical evidence
exerts. When \\\alpha = 0\\ the historical data is ignored entirely;
when \\\alpha = 1\\ it contributes as fully as the current data. In this
example \\\alpha\\ is fixed at 0.5. Power priors were introduced by
Ibrahim and Chen [\[59\]](#ref59) to formally incorporate historical
data into the current analysis through a discounting parameter. Power
priors are applied in clinical trials to leverage results from previous
studies while controlling the degree of historical borrowing.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\textbf{y}\_h \sim
\mathcal{N}(\mu_h, \sigma^2)^\alpha\\ \\\mu = \textbf{X}\beta\\ \\\mu_h
= \textbf{X}\_h\beta\\ \\\alpha = 0.5\\ \\\beta_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the power-prior construction requires two separate likelihood
terms (current and historical) with the historical term raised to an
exponent \\\alpha\\. The DSL does not support multiple likelihood
contributions or power-weighted log-densities.

#### Ground truth and data

We simulate both a current dataset and a historical dataset from the
same true parameters. The intercept and four predictor coefficients are
\\\beta = (1.0, -2.0, 0.5, 1.5, -0.3)\\, the residual standard deviation
is 0.4, and the power discount is \\\alpha = 0.5\\.

``` r

set.seed(42155)
N <- 100
J <- 5
true.beta <- c(1.0, -2.0, 0.5, 1.5, -0.3)
true.sigma <- 0.4
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
Xh <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- rnorm(N, X %*% true.beta, true.sigma)
yh <- rnorm(N, Xh %*% true.beta, true.sigma)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(alpha=0.5, J=J, PGF=PGF, X=X, Xh=Xh, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y,
    yh=yh)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    muh <- tcrossprod(Data$Xh, t(beta))
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(Data$alpha*dnorm(Data$yh, muh, sigma, log=TRUE) +
        dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare estimates to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression with Zellner’s g-Prior

Zellner’s g-prior places a multivariate normal prior on \\\beta\\ with
covariance proportional to \\g\sigma^2(\textbf{X}^T\textbf{X})^{-1}\\,
automatically scaling the prior to the geometry of the design matrix.
The scalar \\g\\ controls the prior’s informativeness: small \\g\\
concentrates the prior near zero while large \\g\\ yields a diffuse
prior that lets the data dominate. A hyper-g prior on \\g\\ allows the
data to determine the appropriate level of shrinkage. See the `dzellner`
documentation for details. Zellner’s g-prior (Zellner
[\[122\]](#ref122)) provides a data-dependent prior covariance for
regression coefficients, with the choice of \\g\\ controlling shrinkage
toward zero. The g-prior is popular in Bayesian model averaging for
variable selection in genomic studies with many candidate predictors.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta \sim \mathcal{N}\_J(0, g \sigma^2 (\textbf{X}^T
\textbf{X})^{-1})\\ \\g \sim \mathcal{HG}(\alpha), \quad \alpha = 3\\
\\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Zellner g-prior (`dzellner`) requires passing the design
matrix \\\textbf{X}\\ as a prior parameter and the hyper-g prior on
\\g\\ involves a non-standard distribution. The DSL does not currently
support data-dependent prior covariance structures.

#### Ground truth and data

We simulate data with an intercept and three predictors (\\J = 4\\) from
known coefficients \\\beta = (1.5, -0.7, 1.1, 0.3)\\ and \\\sigma =
1.0\\. The g-prior parameter \\g\\ and its hyper-g hyperprior are
estimated jointly with \\\beta\\ and \\\sigma\\.

``` r

set.seed(42156)
N <- 200
J <- 4
true.beta <- c(1.5, -0.7, 1.1, 0.3)
true.sigma <- 1.0
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- rnorm(N, X %*% true.beta, true.sigma)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), g0=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.g <- grep("g0", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    g0 <- runif(1)
    sigma <- runif(1)
    return(c(beta, g0, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.g=pos.g,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.g] <- g <- interval(parm[Data$pos.g], 1e-100, Inf)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperprior
    g.prior <- dhyperg(g, alpha=3, log=TRUE)
    ### Log-Prior
    beta.prior <- dzellner(beta, g, sigma, Data$X, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + g.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(1, J), rep(1, 2))
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare estimates to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")
cat("Posterior mean g:", round(mean(fit$Posterior1[, pos.g]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression, Frequentist

By eliminating prior probabilities, this example demonstrates that
lucifer can perform maximum-likelihood estimation in a purely
frequentist framework. The log-posterior equals the log-likelihood, so
the sampler targets the MLE rather than a proper posterior. Although
frequentism is not endorsed here, the example is useful for benchmarking
and for users transitioning from likelihood-based workflows. This
frequentist-in-Bayesian-clothing formulation uses flat priors to
approximate maximum likelihood estimation, following the equivalence
discussed in Gelman et al. [\[48\]](#ref48). Flat-prior linear
regression is useful as a pedagogical bridge between classical and
Bayesian approaches, demonstrating posterior convergence to the MLE.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\

This model has no priors; `LP` equals `LL` directly.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the DSL requires at least one prior distribution on every
parameter. A purely frequentist model with `LP = LL` (no priors) cannot
be expressed in the declarative grammar.

#### Ground truth and data

We simulate data with a known intercept and four predictor coefficients
so that the MLE can be verified against the true values. The residual
standard deviation is set to 0.5.

``` r

set.seed(42149)
N <- 200
J <- 5
true.beta <- c(1.8, -0.9, 1.5, -0.3, 0.7)
true.sigma <- 0.5
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- rnorm(N, X %*% true.beta, true.sigma)
mon.names <- "LL"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    Modelout <- list(LP=LL, Dev=-2*LL, Monitor=LL,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare MLEs to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Polynomial Regression

Polynomial regression extends linear regression by including powers of a
single predictor as additional covariates, allowing the model to capture
nonlinear relationships without requiring a separate nonlinear function
specification. The design matrix contains columns \\x, x^2, \dots, x^P\\
where \\P\\ is the polynomial order. For a more robust extension to
estimating nonlinear relationships between \\\textbf{y}\\ and
\\\textbf{x}\\, see penalized spline regression in section
[link](#penalized.spline). Polynomial regression extends linear
regression with powered terms, with Bayesian regularisation preventing
overfitting as discussed in Gelman et al. [\[48\]](#ref48). Polynomial
regression is used in chemistry to model calibration curves relating
instrument readings to analyte concentrations.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}
\beta\\ \\\textbf{X}\_{i,d} = \textbf{x}^{d-1}\_i, \quad
d=1,\dots,(D+1)\\ \\\textbf{X}\_{i,1} = 1\\ \\\beta_d \sim
\mathcal{N}(0, 1000), \quad d=1,\dots,(D+1)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec

Polynomial regression is structurally identical to linear regression
once the polynomial design matrix \\\textbf{X}\\ is constructed. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL expresses it directly with the pre-built \\\textbf{X}\\.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = X %*% beta
  beta[d] ~ NormalV(0, 1000), d = 1,...,D+1
  sigma ~ HalfCauchy(25)
")
```

#### Ground truth and data

``` r

set.seed(42162)
N <- 200
D <- 2 #Degree of polynomial (quadratic)
true.beta  <- c(1.0, 2.5, -0.8)  # intercept, linear, quadratic
true.sigma <- 0.5

#### Generate predictor and polynomial design matrix
x <- rnorm(N, 0, 1)
X <- cbind(1, matrix(x, N, D))
for (d in 2:D) {X[, d + 1] <- X[, d + 1]^d}

#### Generate response
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- mu.true + rnorm(N, 0, true.sigma)

mon.names <- "LP"
parm.names <- as.parm.names(list(beta = rep(0, D + 1), sigma = 0))
pos.beta  <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta  <- rnorm(1 + Data$D)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(D = D, N = N, PGF = PGF, mon.names = mon.names,
    parm.names = parm.names, pos.beta = pos.beta, pos.sigma = pos.sigma,
    x = x, y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    X <- matrix(Data$x, Data$N, Data$D)
    for (d in 2:Data$D) {X[,d] <- X[,d]^ d}
    X <- cbind(1,X)
    mu <- tcrossprod(X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,D+1), 1)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "sigma"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Quantile Regression

Quantile regression estimates conditional quantiles of the response
distribution rather than the conditional mean, using an asymmetric
Laplace likelihood parameterized by the quantile level \\q \in (0,1)\\.
This approach is particularly useful when the tails of the distribution
or the effects of covariates at different quantiles are of substantive
interest. Bayesian quantile regression was introduced by Yu and Moyeed
[\[116\]](#ref116) using the asymmetric Laplace likelihood. Quantile
regression is applied in education to study how school quality
differentially affects students at different points of the achievement
distribution.

#### Form

\\\textbf{y} \sim \mathcal{N}(\phi, \sigma^2)\\ \\\phi = \frac{(1 -
2P)}{P(1 - P)} \zeta + \mu\\ \\\mu = \textbf{X} \beta\\ \\\sigma =
\frac{P (1 - P) \tau}{2 \zeta}\\ \\\beta \sim \mathcal{N}(0, 1000)\\
\\\tau \sim \mathcal{HC}(25)\\ \\\zeta \sim \mathcal{EXP}(\tau)\\ where
\\P\\ is the user-specified quantile in \\(0,1)\\.

#### model_spec

The asymmetric Laplace mixture representation used here, which
introduces per-observation latent variables \\\zeta_i\\, is not
expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The manual Model function is therefore the recommended approach.

#### Ground truth and data

``` r

set.seed(42152)
N <- 300
J <- 4
true.beta <- c(2.0, 1.5, -0.8, 0.4)
P <- 0.5 #Quantile in (0,1)

#### Design matrix with intercept
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 1), N, J - 1))
for (j in 2:J) X[, j] <- CenterScale(X[, j])

#### Generate response with heteroscedastic errors
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- mu.true + rnorm(N, 0, 1.5)

mon.names <- "LP"
parm.names <- as.parm.names(list(beta = rep(0, J), tau = 0, zeta = rep(0, N)))
pos.beta <- grep("beta", parm.names)
pos.tau  <- grep("tau", parm.names)
pos.zeta <- grep("zeta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    tau  <- runif(1)
    zeta <- rexp(Data$N, tau)
    return(c(beta, tau, zeta))
    }
Data <- list(J = J, N = N, P = P, PGF = PGF, X = X, mon.names = mon.names,
    parm.names = parm.names, pos.beta = pos.beta, pos.tau = pos.tau,
    pos.zeta = pos.zeta, y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    zeta <- interval(parm[Data$pos.zeta], 1e-100, Inf)
    parm[Data$pos.zeta] <- zeta
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    zeta.prior <- sum(dexp(zeta, tau, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    phi <- (1 - 2*Data$P) / (Data$P*(1 - Data$P))*zeta + mu
    sigma <- (Data$P*(1 - Data$P)*tau) / (2*zeta)
    LL <- sum(dnorm(Data$y, phi, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + tau.prior + zeta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(phi), phi, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1, rep(1,N))
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "tau"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Robust Regression

Robust regression replaces the normal likelihood with a Student-t
distribution, introducing a degrees-of-freedom parameter \\\nu\\ that
controls tail weight. When \\\nu\\ is small the model accommodates
outliers far better than ordinary normal regression; as \\\nu \to
\infty\\ it converges to standard least squares. As an alternative
approach to robust regression, consider Laplace regression (see section
[link](#laplace.reg)). Robust regression using the Student-t likelihood
was developed by Lange et al. [\[70\]](#ref70), providing resistance to
outliers by estimating the degrees of freedom. Student-t regression is
applied in astronomy where outlying observations from instrument errors
or transient phenomena would distort normal-theory estimates.

#### Form

\\\textbf{y} \sim \mathrm{t}(\mu, \sigma^2, \nu)\\ \\\mu =
\textbf{X}\beta\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\sigma \sim \mathcal{HC}(25)\\ \\\nu \sim
\mathcal{HC}(25)\\

#### model_spec

The `StudentT` distribution is registered in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL, so this model can be expressed directly.

``` r

spec <- model_spec("
  y ~ StudentT(nu, mu, sigma)
  mu = X %*% beta
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  sigma ~ HalfCauchy(25)
  nu ~ HalfCauchy(25)
")
```

#### Ground truth and data

``` r

set.seed(42153)
N <- 200
J <- 4
true.beta  <- c(2.0, -1.5, 0.8, 0.3)
true.sigma <- 1.0
true.nu    <- 4

#### Design matrix with intercept
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 1), N, J - 1))

#### Generate response with heavy-tailed errors
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- mu.true + true.sigma * rt(N, df = true.nu)

mon.names <- "LP"
parm.names <- as.parm.names(list(beta = rep(0, J), sigma = 0, nu = 0))
pos.beta  <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.nu    <- grep("nu", parm.names)
PGF <- function(Data) {
    beta  <- rnorm(Data$J)
    sigma <- runif(1)
    nu    <- runif(1)
    return(c(beta, sigma, nu))
    }
Data <- list(J = J, PGF = PGF, X = X, mon.names = mon.names,
    parm.names = parm.names, pos.beta = pos.beta, pos.sigma = pos.sigma,
    pos.nu = pos.nu, y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[1:Data$J]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    nu.prior <- dhalfcauchy(nu, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dst(Data$y, mu, sigma, nu, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior + nu.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rst(length(mu), mu, sigma, nu), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1, 5)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "sigma", "nu"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.nu]    <- true.nu[seq_along(pos.nu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Weighted Regression

Weighted regression accommodates heteroscedastic data by scaling the
normal likelihood of each observation by a known weight \\w_i\\; the
effective variance for observation \\i\\ is \\\sigma^2 / w_i\\, so
observations with larger weights exert more influence on the posterior.
Here, weights are applied to a standard linear regression to demonstrate
the mechanism. Weighted regression addresses heteroscedastic errors by
assigning known weights to observations, placed in a Bayesian framework
by Gelman et al. [\[48\]](#ref48). Weighted regression is used in
meta-analysis to combine study-level estimates where each study
contributes in proportion to its precision.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\text{LL} = \sum\_{i=1}^{N} w_i \log \mathcal{N}(\textbf{y}\_i \mid
\mu_i, \sigma^2)\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

Weighted regression cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires observation-level scaling of the log-likelihood
contributions (\\w_i \log p(y_i \mid \mu_i, \sigma)\\), which is a
custom modification to the likelihood evaluation that falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42158)
N <- 200
J <- 4

#### True parameter values
true.beta <- c(3.0, -1.2, 0.7, 0.4)  # intercept + 3 predictors
true.sigma <- 1.0

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Create heterogeneous weights (some observations down-weighted)
w <- rep(1, N)
w[sample(N, 20)] <- 0.1     # 20 low-influence observations
w <- w * (N / sum(w))         # normalise so sum(w) = N

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, w=w,
    y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(Data$w * dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Generalized linear models: binary and binomial

### Binary Log-Log Link Mixture

A weighted mixture of the log-log and complementary log-log link
functions is used, where \\\alpha\\ is the mixing weight. Since both
link functions are asymmetric (unlike the symmetric logit and probit
links), it may be unknown a priori which direction of asymmetry best
describes the data; the mixture lets the data decide by estimating
\\\alpha\\ alongside the regression coefficients. Log-log link mixtures
extend standard binary regression by accommodating asymmetric response
curves, following Pregibon [\[94\]](#ref94) and the Bayesian mixture
framework of Frühwirth-Schnatter [\[42\]](#ref42). In reliability
engineering, log-log link mixtures model failure probabilities when the
population contains a subgroup of items resistant to the stress factor.

#### Form

\\\textbf{y} \sim \mathcal{BERN}(\eta)\\ \\\eta = \alpha
\exp(-\exp(\mu)) + (1 - \alpha) (1 - \exp(-\exp(\mu)))\\ \\\mu =
\textbf{X} \beta\\ \\\alpha \sim \mathcal{U}(0, 1)\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the likelihood involves a weighted mixture of two
nonstandard link functions (`invloglog` and `invcloglog`) applied to the
linear predictor. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
framework does not support user-defined link mixtures or the
`invloglog`/`invcloglog` transformations. The manual `Model` function
below handles this directly.

#### Ground truth and data

``` r

set.seed(42117)
N <- 100
J <- 3

#### True parameter values
true.alpha <- 0.6
true.beta <- c(-0.5, 1.0, -0.7)

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu.true <- tcrossprod(X, t(true.beta))
eta.true <- true.alpha * invloglog(mu.true) +
    (1 - true.alpha) * invcloglog(mu.true)
y <- rbern(N, eta.true)

#### Assemble Data list
mon.names <- c("LP", "alpha")
parm.names <- as.parm.names(list(beta=rep(0,J), logit.alpha=0))
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    logit.alpha <- rnorm(1)
    return(c(beta, logit.alpha))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[Data$J+1] <- alpha <- interval(parm[Data$J+1], -700, 700)
    beta <- parm[1:Data$J]
    ### Log-Prior
    alpha.prior <- dunif(alpha, 0, 1, log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    eta <- alpha*invloglog(mu) + (1-alpha)*invcloglog(mu)
    LL <- sum(dbern(Data$y, eta, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,alpha),
        yhat=rbern(length(eta), eta), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 0)
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
post.beta <- fit$Summary2[1:J, "Mean"]
post.alpha <- mean(fit$Monitor[, 2])
cat("True beta: ", true.beta, "\n")
cat("Posterior mean:", round(post.beta, 3), "\n")
cat("True alpha:", true.alpha, " Posterior mean:", round(post.alpha, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("alpha", parm.names)] <- true.alpha[seq_along(grep("alpha", parm.names))]
ground_truth[grep("beta", parm.names)]  <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binary Logit

Binary logistic regression models a binary outcome as a Bernoulli trial
whose success probability is linked to a linear predictor through the
inverse-logit (expit) function. This is the standard approach for
classification problems where the response takes values 0 or 1, and the
log-odds of success are assumed linear in the covariates. The logistic
regression model was introduced by Berkson [\[8\]](#ref8) and its
Bayesian treatment systematised by Albert and Chib [\[2\]](#ref2).
Binary logit models are ubiquitous in medicine for predicting disease
diagnosis from clinical biomarkers.

#### Form

\\\textbf{y} \sim \mathcal{BERN}(\eta)\\ \\\eta = \frac{1}{1 +
\exp(-\mu)}\\ \\\mu = \textbf{X} \beta\\ \\\beta_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\

#### model_spec

The binary logit model is expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The inverse-logit link appears as a deterministic statement using
any valid R expression.

``` r

spec <- model_spec("
  y ~ Bernoulli(eta)
  eta = 1 / (1 + exp(-(X %*% beta)))
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
")
```

#### Data

``` r

set.seed(42116)
N <- 500
J <- 3
true.beta <- c(0.5, -1.2, 0.8)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
eta <- 1 / (1 + exp(-mu))
y <- rbinom(N, 1, eta)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    return(beta)
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[1:Data$J]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    eta <- invlogit(mu)
    LL <- sum(dbern(Data$y, eta, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbern(length(eta), eta), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0,J)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binary Probit

Binary probit regression models a binary outcome using the standard
normal CDF as the link function, mapping the linear predictor
\\\textbf{X}\beta\\ to a probability \\p = \Phi(\textbf{X}\beta)\\. The
probit link is the natural choice when the latent variable generating
the binary outcome is assumed to be normally distributed. The probit
model dates to Bliss [\[14\]](#ref14), with the modern Bayesian data
augmentation approach due to Albert and Chib [\[2\]](#ref2). Probit
models are standard in bioassay for estimating the dose at which 50% of
organisms exhibit a response (LD50).

#### Form

\\\textbf{y} \sim \mathcal{BERN}(\textbf{p})\\ \\\textbf{p} =
\phi(\mu)\\ \\\mu = \textbf{X} \beta \in \[-10,10\]\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ where \\\phi\\ is the CDF of
the standard normal distribution, and \\J\\=3.

#### model_spec() notation

The binary probit model is expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The normal CDF link appears as a deterministic statement.

``` r

spec <- model_spec("
  y ~ Bernoulli(p)
  p = pnorm(X %*% beta)
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
")
```

#### Ground truth and data

``` r

set.seed(42118)
N <- 500
J <- 3
true.beta <- c(0.3, -0.9, 0.6)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
mu <- pmin(pmax(mu, -10), 10)
p <- pnorm(mu)
y <- rbinom(N, 1, p)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    return(beta)
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[1:Data$J]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    mu <- interval(mu, -10, 10, reflect=FALSE)
    p <- pnorm(mu)
    LL <- sum(dbern(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbern(length(p), p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0,J)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binary Robit

The robit model extends binary probit by replacing the normal CDF with
the Student-t CDF as the link function, introducing a degrees-of-freedom
parameter \\\nu\\ that controls tail thickness. When \\\nu\\ is small
the link has heavier tails than probit, making the model more robust to
outlying covariate values; as \\\nu \to \infty\\ the robit converges to
probit. The robit model replaces the normal latent variable in probit
regression with a Student-t distribution, as proposed by Liu
[\[75\]](#ref75), providing robustness to outlying observations. Robit
models are applied in credit scoring where misclassified defaulters
create influential outliers that distort standard logistic predictions.

#### Form

\\\textbf{y} \sim \mathcal{BERN}(\textbf{p})\\ \\\textbf{p} =
\textbf{T}\_\nu(\mu)\\ \\\mu = \textbf{X} \beta \in \[-10,10\]\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\nu \sim
\mathcal{U}(1 \times 10^{-100}, 1000)\\ where \\\textbf{T}\_\nu\\ is the
CDF of the standard t-distribution with \\\nu\\ degrees of freedom.

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Student-t CDF link (`pst`) is a custom distribution function
not available in the DSL.

#### Ground truth and data

``` r

set.seed(42119)
N <- 500
J <- 3
true.beta <- c(0.5, -1.0, 0.7)
true.nu <- 7
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
mu <- pmin(pmax(mu, -10), 10)
p <- pt(mu, df = true.nu)
y <- rbinom(N, 1, p)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), nu=0))
pos.beta <- grep("beta", parm.names)
pos.nu <- grep("nu", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    nu <- runif(1, 1, 100)
    return(c(beta, nu))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.nu=pos.nu, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, 1000)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    nu.prior <- dunif(nu, 1e-100, 1000, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    mu <- interval(mu, -10, 10, reflect=FALSE)
    p <- pst(mu, nu=nu)
    LL <- sum(dbern(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + nu.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbern(length(p), p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 5)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.nu]   <- true.nu[seq_along(pos.nu)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binomial Logit

Binomial logistic regression extends the binary logit model to grouped
data where each observation records the number of successes \\y_i\\ out
of \\n_i\\ trials, with the success probability linked to a linear
predictor through the logistic function. Binomial logistic regression
extends binary logit to grouped data and was formalised within the
generalised linear model framework by Nelder and Wedderburn
[\[89\]](#ref89). Binomial logit models are used in entomology to
analyse dose-mortality curves from insecticide resistance bioassays.

#### Form

\\\textbf{y} \sim \mathcal{BIN}(\textbf{p}, \textbf{n})\\ \\\textbf{p} =
\frac{1}{1 + \exp(-\mu)}\\ \\\mu = \textbf{X} \beta\\ \\\beta_j \sim
\mathcal{N}(0,1000), \quad j=1,\dots,J\\

#### model_spec() notation

The binomial logit model is expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL using the inverse-logit deterministic transformation.

``` r

spec <- model_spec("
  y ~ Binomial(n, p)
  p = 1 / (1 + exp(-(X %*% beta)))
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
")
```

#### Ground truth and data

``` r

set.seed(42120)
N <- 200
J <- 3
true.beta <- c(-0.5, 1.0, -0.7)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
p <- 1 / (1 + exp(-mu))
n <- sample(20:100, N, replace = TRUE)
y <- rbinom(N, n, p)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
PGF <- function(Data) return(rnorm(Data$J))
Data <- list(J=J, N=N, PGF=PGF, X=X, n=n, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[1:Data$J]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    p <- invlogit(mu)
    LL <- sum(dbinom(Data$y, Data$n, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(length(p), Data$n, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0,J)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binomial Probit

Binomial probit regression models grouped binary data using the normal
CDF link; each observation is the count of successes from a known number
of trials with success probability \\p = \Phi(\textbf{X}\beta)\\. The
binomial probit model uses a normal CDF link for grouped binary data,
with Bayesian estimation following Albert and Chib [\[2\]](#ref2). In
sensory science, binomial probit models estimate detection thresholds
from triangle test data where panels of judges classify stimuli.

#### Form

\\\textbf{y} \sim \mathcal{BIN}(\textbf{p}, \textbf{n})\\ \\\textbf{p} =
\phi(\mu)\\ \\\mu = \textbf{X} \beta \in \[-10,10\]\\ \\\beta_j \sim
\mathcal{N}(0,1000), \quad j=1,\dots,J\\ where \\\phi\\ is the CDF of
the standard normal distribution.

#### model_spec() notation

The binomial probit model is expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL using the normal CDF as a deterministic transformation.

``` r

spec <- model_spec("
  y ~ Binomial(n, p)
  p = pnorm(X %*% beta)
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
")
```

#### Ground truth and data

``` r

set.seed(42121)
N <- 200
J <- 3
true.beta <- c(-0.3, 0.8, -0.5)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
mu <- pmin(pmax(mu, -10), 10)
p <- pnorm(mu)
n <- sample(20:100, N, replace = TRUE)
y <- rbinom(N, n, p)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
PGF <- function(Data) return(rnorm(Data$J))
Data <- list(J=J, N=N, PGF=PGF, X=X, n=n, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[1:Data$J]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    mu <- interval(mu, -10, 10, reflect=FALSE)
    p <- pnorm(mu)
    LL <- sum(dbinom(Data$y, Data$n, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(length(p), Data$n, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0,J)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Binomial Robit

Binomial robit regression uses a Student-t CDF link for grouped binary
data, providing robustness to outlying covariate values compared to the
binomial logit and probit models. The binomial robit model applies the
robust Student-t link of Liu [\[75\]](#ref75) to grouped binary data,
down-weighting aberrant group-level proportions. Robit links are useful
in epidemiology when a few study sites report anomalous disease
prevalence due to local reporting artefacts.

#### Form

\\\textbf{y} \sim \mathcal{BIN}(\textbf{p}, \textbf{n})\\ \\\textbf{p} =
\textbf{T}\_\nu(\mu)\\ \\\mu = \textbf{X} \beta \in \[-10,10\]\\
\\\beta_j \sim \mathcal{N}(0,1000), \quad j=1,\dots,J\\ \\\nu \sim
\mathcal{U}(5, 10)\\ where \\\textbf{T}\_\nu\\ is the CDF of the
standard t-distribution with \\\nu\\ degrees of freedom.

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Student-t CDF link (`pst`) is a custom distribution function
not available in the DSL.

#### Ground truth and data

``` r

set.seed(42122)
N <- 200
J <- 3
true.beta <- c(0.4, -0.8, 0.5)
true.nu <- 7
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu <- tcrossprod(X, t(true.beta))
mu <- pmin(pmax(mu, -10), 10)
p <- pt(mu, df = true.nu)
n <- sample(20:100, N, replace = TRUE)
y <- rbinom(N, n, p)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), nu=0))
pos.beta <- grep("beta", parm.names)
pos.nu <- grep("nu", parm.names)
PGF <- function(Data) return(c(rnorm(Data$J), runif(1,5,10)))
Data <- list(J=J, N=N, PGF=PGF, X=X, n=n, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.nu=pos.nu, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 5, 10)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    nu.prior <- dunif(nu, 5, 10, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    mu <- interval(mu, -10, 10, reflect=FALSE)
    p <- pst(mu, nu=nu)
    LL <- sum(dbinom(Data$y, Data$n, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + nu.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(length(p), Data$n, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 7)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.nu]   <- true.nu[seq_along(pos.nu)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Overdispersed Logistic Regression **\[NEW\]**

Overdispersed logistic regression extends standard binomial regression
by introducing observation-level random effects on the logit scale. Each
observation \\i\\ receives a latent logit-probability
\\\text{logit}(p_i)\\ drawn from a normal distribution centred on the
linear predictor \\\textbf{X}\_i\beta\\ with variance \\\sigma_p^2\\, so
that the effective probability varies more widely than the binomial
model alone would permit. This extra-binomial variation is common in
grouped binary data where unobserved heterogeneity inflates the variance
beyond what the binomial assumption supports. The approach was
introduced by Williams [\[143\]](#ref143) for dose-response assays and
provides a direct Bayesian alternative to quasi-likelihood corrections.
Overdispersed logistic models are applied in toxicology to analyse
dose-response experiments with extra-binomial variation among
experimental units.

#### Form

\\y_i \sim \mathcal{B}(K_i, p_i), \quad i=1,\dots,N\\
\\\text{logit}(p_i) \sim \mathcal{N}(\textbf{X}\_i\beta, \sigma_p^2)\\
\\\beta_j \sim \mathcal{N}(0, 100), \quad j=1,\dots,J\\ \\\sigma_p \sim
\mathcal{HC}(5)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the observation-level random effects introduce \\N\\ latent
logit-probability parameters whose hierarchical prior depends on the
linear predictor, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42225)
N <- 100
K <- 20
J <- 3
true.beta <- c(0.5, 1, -0.5)
true.sigma.p <- 0.5

#### Simulate overdispersed binomial data
x1 <- rnorm(N)
x2 <- rnorm(N)
X <- cbind(1, x1, x2)
eta <- as.vector(X %*% true.beta) + rnorm(N, 0, true.sigma.p)
p <- plogis(eta)
y <- rbinom(N, K, p)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0, J), sigma.p=0,
    logit.p=rep(0, N)))
pos.beta <- grep("beta", parm.names)
pos.sigma.p <- grep("sigma.p", parm.names)
pos.logit.p <- grep("logit.p", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma.p <- runif(1, 0.1, 2)
    logit.p <- rnorm(Data$N, Data$X %*% beta, sigma.p)
    return(c(beta, sigma.p, logit.p))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta,
    pos.logit.p=pos.logit.p, pos.sigma.p=pos.sigma.p, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma.p <- interval(parm[Data$pos.sigma.p], 1e-100, Inf)
    parm[Data$pos.sigma.p] <- sigma.p
    logit.p <- parm[Data$pos.logit.p]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 100, log=TRUE))
    sigma.p.prior <- dhalfcauchy(sigma.p, 5, log=TRUE)
    logit.p.prior <- sum(dnorm(logit.p, as.vector(Data$X %*% beta),
        sigma.p, log=TRUE))
    ### Log-Likelihood
    p <- plogis(logit.p)
    LL <- sum(dbinom(Data$y, Data$K, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.p.prior + logit.p.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(Data$N, Data$K, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1, rep(0, N))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Generalized linear models: count data

### Conway-Maxwell-Poisson Regression **\[NEW\]**

The Conway-Maxwell-Poisson (COM-Poisson) distribution generalises the
Poisson by adding a dispersion parameter \\\nu\\ that controls tail
weight: when \\\nu = 1\\ the model reduces to Poisson, \\\nu \< 1\\
produces overdispersion, and \\\nu \> 1\\ yields underdispersion. The
probability mass function includes a normalising constant \\Z(\lambda,
\nu) = \sum\_{s=0}^{\infty} \lambda^s / (s!)^\nu\\ that has no closed
form and must be evaluated by truncated summation. Shmueli et
al. [\[155\]](#ref155) developed the COM-Poisson regression framework
and demonstrated its flexibility for count data that violates the
Poisson equidispersion assumption. COM-Poisson regression is applied in
manufacturing quality control to model defect counts that exhibit
underdispersion due to process constraints.

#### Form

\\P(Y_i = y_i) = \frac{\lambda_i^{y_i}}{(y_i !)^\nu \\ Z(\lambda_i,
\nu)}\\ \\Z(\lambda, \nu) = \sum\_{s=0}^{S} \frac{\lambda^s}{(s!)^\nu}\\
\\\log(\lambda_i) = \textbf{x}\_i^\top \boldsymbol{\beta}\\ \\\beta_j
\sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\nu \sim
\mathcal{HC}(25), \quad \nu \in (0.01, 10)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the COM-Poisson likelihood requires computing a normalising
constant \\Z(\lambda, \nu)\\ by truncated summation over a series whose
length depends on the current parameter values. This custom normalising
constant falls outside the distributional registry of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42206)
N <- 150
J <- 2

#### True parameter values
true.beta <- c(1.5, 0.8)
true.nu <- 1.5

#### Covariates
x <- matrix(rnorm(N * J), nrow=N, ncol=J)

#### COM-Poisson sampler via PMF inversion
rcomp <- function(n, lambda, nu, S.max=200) {
    y <- integer(n)
    for (i in 1:n) {
        s <- 0:S.max
        log.pmf <- s * log(lambda[i]) - nu * lgamma(s + 1)
        log.pmf <- log.pmf - max(log.pmf)
        pmf <- exp(log.pmf)
        pmf <- pmf / sum(pmf)
        y[i] <- sample(s, 1, prob=pmf)
    }
    return(y)
    }

#### Generate response
lambda.true <- exp(x %*% true.beta)
y <- rcomp(N, lambda.true, true.nu)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), nu=0))
pos.beta <- grep("beta", parm.names)
pos.nu <- grep("nu", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    nu <- rhalfcauchy(1, 5)
    return(c(beta, nu))
    }
Data <- list(N=N, J=J, S.max=200, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.nu=pos.nu, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    nu <- interval(parm[Data$pos.nu], 1e-2, 10)
    parm[Data$pos.nu] <- nu
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    nu.prior <- dhalfcauchy(nu, 25, log=TRUE)
    ### Log-Likelihood
    log.lambda <- Data$x %*% beta
    lambda <- exp(log.lambda)
    ### Compute log(Z) for each observation using log-sum-exp
    s <- 0:Data$S.max
    LL <- 0
    for (i in 1:Data$N) {
        log.terms <- s * log.lambda[i] - nu * lgamma(s + 1)
        max.term <- max(log.terms)
        log.Z <- max.term + log(sum(exp(log.terms - max.term)))
        LL <- LL + Data$y[i] * log.lambda[i] -
            nu * lgamma(Data$y[i] + 1) - log.Z
    }
    ### Log-Posterior
    LP <- LL + beta.prior + nu.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(Data$N, lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
for (j in 1:J)
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(fit$Summary2[pos.beta[j], "Mean"], 3), "\n")
cat("nu -- true:", true.nu,
    " post. mean:", round(fit$Summary2[pos.nu, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.nu]   <- true.nu[seq_along(pos.nu)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Hurdle Model (Zero-Altered Poisson) **\[NEW\]**

The hurdle model, introduced by Mullahy [\[137\]](#ref137) and discussed
extensively by Zuur et al. [\[138\]](#ref138), separates the mechanism
that generates zeros from the mechanism that generates positive counts.
A logistic sub-model governs the binary outcome of whether a zero or a
positive count is observed, while a truncated Poisson sub-model
describes the distribution of positive counts conditional on crossing
the zero hurdle. This two-part structure is more flexible than a
standard Poisson because the zero probability is decoupled from the
count intensity; the zero mass can be larger or smaller than
\\\exp(-\lambda)\\. Six regression parameters link two covariates plus
intercepts to each sub-model. Hurdle models are applied in health
economics to model healthcare utilisation where many individuals have
zero visits.

#### Form

\\P(y_i = 0) = 1 - \pi_i\\ \\P(y_i = k \mid y_i \> 0) = \pi_i \cdot
\frac{\text{Pois}(k;\\ \lambda_i)}{1 - \exp(-\lambda_i)}, \quad
k=1,2,\dots\\ \\\text{logit}(\pi_i) = \textbf{X}\_i \alpha\\
\\\log(\lambda_i) = \textbf{X}\_i \beta\\ \\\alpha_j \sim \mathcal{N}(0,
100), \quad j=1,\dots,3\\ \\\beta_j \sim \mathcal{N}(0, 100), \quad
j=1,\dots,3\\

#### model_spec() notation

The hurdle model cannot be expressed in
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because its piecewise likelihood combines a logistic zero component with
a truncated Poisson component for positive counts. This two-part
structure with distinct link functions and a truncation normalisation
constant falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42216)
N <- 300
J <- 3

#### True parameter values
true.alpha <- c(-0.5, 1.0, 0.3)
true.beta  <- c(1.0, 0.5, -0.2)

#### Generate covariates and design matrix
x1 <- rnorm(N)
x2 <- rnorm(N)
X <- cbind(1, x1, x2)

#### Generate hurdle data
pi.true <- plogis(X %*% true.alpha)
lambda.true <- exp(X %*% true.beta)
y <- integer(N)
for (i in 1:N) {
    if (runif(1) < pi.true[i]) {
        ## Draw from truncated Poisson (reject zeros)
        yi <- 0L
        while (yi == 0L) yi <- rpois(1, lambda.true[i])
        y[i] <- yi
    } else {
        y[i] <- 0L
    }
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J), beta=rep(0,J)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    beta <- rnorm(Data$J)
    return(c(alpha, beta))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 100, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 100, log=TRUE))
    ### Log-Likelihood
    pi.i <- plogis(tcrossprod(Data$X, t(alpha)))
    lambda.i <- exp(tcrossprod(Data$X, t(beta)))
    LL <- 0
    idx0 <- which(Data$y == 0)
    idx1 <- which(Data$y > 0)
    ## Zero component
    if (length(idx0) > 0)
        LL <- LL + sum(log(1 - pi.i[idx0]))
    ## Positive-count component (truncated Poisson)
    if (length(idx1) > 0)
        LL <- LL + sum(log(pi.i[idx1]) +
            dpois(Data$y[idx1], lambda.i[idx1], log=TRUE) -
            log(1 - exp(-lambda.i[idx1])))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior
    ### yhat: approximate hurdle draw
    yhat <- ifelse(runif(Data$N) < pi.i, rpois(Data$N, lambda.i), 0L)
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0, 2 * J)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha, "\n")
cat("         post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("beta  -- true:", true.beta, "\n")
cat("         post. mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Negative Binomial Regression

The negative binomial regression extends Poisson regression to handle
overdispersed count data by introducing a dispersion parameter
\\\kappa\\ (also called `size`). When \\\kappa \to \infty\\ the model
reduces to Poisson regression; smaller values of \\\kappa\\ indicate
greater overdispersion. The mean-dispersion parameterization used here
passes `mu` and `size` directly to R’s `dnbinom`. This example was
originally contributed by Jim Robison-Cox. The negative binomial
regression model handles overdispersed count data and was systematised
in a Bayesian context by Hilbe [\[56\]](#ref56). Negative binomial
models are standard in ecology for modelling species abundance counts
where variance exceeds the mean due to spatial aggregation.

#### Form

\\\textbf{y} \sim \mathcal{NB}(\mu, \kappa)\\ \\p =
\frac{\kappa}{\kappa + \mu}\\ \\\mu = \exp(\textbf{X} \beta)\\ \\\beta_j
\sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\kappa \sim
\mathcal{HC}(25) \in (0,\infty\]\\

#### model_spec

The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry entry for `NegBinomial` uses the `(size, prob)`
parameterization mapping to `dnbinom(size, prob)`. Because this model
uses the `(size, mu)` alternative parameterization through
`dnbinom(y, size=kappa, mu=mu)`, the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL cannot express it directly without converting to the probability
form. The manual Model function is therefore the recommended approach
for this model.

#### Data

``` r

set.seed(42163)
N <- 300
J <- 4
true.beta <- c(1.5, 0.3, -0.5, 0.2)
true.kappa <- 3
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 0.5), N, J - 1))
mu.true <- exp(tcrossprod(X, t(true.beta)))
y <- rnbinom(N, size = true.kappa, mu = mu.true)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), kappa=0))
pos.beta <- grep("beta", parm.names)
pos.kappa <- grep("kappa", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    kappa <- runif(1)
    return(c(beta, kappa))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.kappa=pos.kappa, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$J + 1] <- kappa <- interval(parm[Data$pos.kappa],
        .Machine$double.xmin, Inf)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    kappa.prior <- dhalfcauchy(kappa, 25, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(exp(tcrossprod(Data$X, t(beta))))
    #p <- kappa / (kappa + mu)
    LL <- sum(dnbinom(Data$y, size=kappa, mu=mu, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + kappa.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnbinom(length(mu), size=kappa, mu=mu), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "kappa"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.kappa] <- true.kappa[seq_along(pos.kappa)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Panel, Autoregressive Poisson

The panel autoregressive Poisson model accommodates count data observed
across \\N\\ cross-sectional units over \\T\\ time periods.
Unit-specific intercepts \\\alpha_i\\ are drawn from a hierarchical
normal distribution, a common covariate effect \\\beta\\ shifts the
log-rate, and the autoregressive coefficient \\\rho\\ captures temporal
persistence through the lagged log-counts. The panel autoregressive
Poisson model extends count regression to longitudinal data with
temporal dependence, following Zeger [\[119\]](#ref119). Panel Poisson
models are used in criminology to analyse monthly crime counts across
police precincts over multiple years.

#### Form

\\\textbf{Y} \sim \mathcal{P}(\Lambda)\\ \\\Lambda\_{1:N,1} =
\exp(\alpha + \beta \textbf{x})\\ \\\Lambda\_{1:N,t} = \exp(\alpha +
\beta \textbf{x} + \rho \log(\textbf{Y}\_{1:N,t-1})), \quad
t=2,\dots,T\\ \\\alpha_i \sim \mathcal{N}(\alpha\_\mu,
\alpha^2\_\sigma), \quad i=1,\dots,N\\ \\\alpha\_\mu \sim \mathcal{N}(0,
1000)\\ \\\alpha\_\sigma \sim \mathcal{HC}(25)\\ \\\beta \sim
\mathcal{N}(0, 1000)\\ \\\rho \sim \mathcal{N}(0, 1000)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the autoregressive component \\\rho \log(Y\_{t-1})\\ introduces
recursive time dependence on observed lagged counts, which falls outside
the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42165)
N <- 10
T <- 10

#### True parameter values
true.alpha.mu    <- 2.0
true.alpha.sigma <- 0.5
true.beta        <- 0.5
true.rho         <- 0.5

#### Simulate panel data
alpha <- rnorm(N, true.alpha.mu, true.alpha.sigma)
x <- runif(N, 0, 1)
Y <- matrix(NA, N, T)
Y[,1] <- rpois(N, exp(alpha + true.beta * x))
for (t in 2:T) {
    Y[,t] <- rpois(N, exp(alpha + true.beta * x +
        true.rho * log(pmax(Y[,t-1], 1))))
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,N), alpha.mu=0,
    alpha.sigma=0, beta=0, rho=0))
pos.alpha <- 1:N
pos.alpha.mu <- grep("alpha.mu", parm.names)
pos.alpha.sigma <- grep("alpha.sigma", parm.names)
pos.beta <- grep("beta", parm.names)
pos.rho <- grep("rho", parm.names)
PGF <- function(Data) {
    alpha.mu <- rnorm(1)
    alpha.sigma <- runif(1)
    alpha <- rnorm(Data$N, alpha.mu, alpha.sigma)
    beta <- rnorm(1)
    rho <- rnorm(1)
    return(c(alpha, alpha.mu, alpha.sigma, beta, rho))
    }
Data <- list(N=N, PGF=PGF, T=T, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.alpha.mu=pos.alpha.mu,
    pos.alpha.sigma=pos.alpha.sigma, pos.beta=pos.beta, pos.rho=pos.rho,
    x=x)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    alpha.mu <- parm[Data$pos.alpha.mu]
    alpha.sigma <- interval(parm[Data$pos.alpha.sigma], 1e-100, Inf)
    parm[Data$pos.alpha.sigma] <- alpha.sigma
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    rho <- parm[Data$pos.rho]
    ### Log-Hyperprior
    alpha.mu.prior <- dnormv(alpha.mu, 0, 1000, log=TRUE)
    alpha.sigma.prior <- dhalfcauchy(alpha.sigma, 25, log=TRUE)
    ### Log-Prior
    alpha.prior <- sum(dnorm(alpha, alpha.mu, alpha.sigma, log=TRUE))
    beta.prior <- dnormv(beta, 0, 1000, log=TRUE)
    rho.prior <- dnormv(rho, 0, 1000, log=TRUE)
    ### Log-Likelihood
    Lambda <- Data$Y
    Lambda[,1] <- exp(alpha + beta*Data$x)
    Lambda[,2:Data$T] <- exp(alpha + beta*Data$x +
        rho*log(pmax(Data$Y[,1:(Data$T-1)], 1)))
    LL <- sum(dpois(Data$Y, Lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + alpha.mu.prior + alpha.sigma.prior +
        beta.prior + rho.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(prod(dim(Lambda)), Lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,N), 0, 1, 0, 0)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True alpha.mu:   ", true.alpha.mu, "\n")
cat("Posterior mean:   ", round(fit$Summary2[pos.alpha.mu, "Mean"], 3), "\n")
cat("True beta:        ", true.beta, "\n")
cat("Posterior mean:   ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True rho:         ", true.rho, "\n")
cat("Posterior mean:   ", round(fit$Summary2[pos.rho, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.rho]  <- true.rho[seq_along(pos.rho)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Poisson Regression

Poisson regression models count data through the log link, so that the
linear predictor \\\textbf{X}\beta\\ is exponentiated to produce the
rate parameter \\\lambda\\. This is the canonical generalized linear
model for non-negative integer responses when the variance is assumed
equal to the mean. When overdispersion is present, consider the negative
binomial alternative in section [link](#negbin.reg). Poisson regression
was established within the generalised linear model framework by Nelder
and Wedderburn [\[89\]](#ref89). Poisson models are standard in
epidemiology for analysing disease incidence rates as a function of risk
factors.

#### Form

\\\textbf{y} \sim \mathcal{P}(\lambda)\\ \\\lambda =
\exp(\textbf{X}\beta)\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\

#### model_spec

The Poisson regression is directly expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL.

``` r

spec <- model_spec("
  y ~ Poisson(lambda)
  lambda = exp(X %*% beta)
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
")
```

#### Data

``` r

set.seed(42161)
N <- 300
J <- 3
true.beta <- c(1.0, 0.3, -0.5)
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 0.5), N, J - 1))
lambda.true <- exp(tcrossprod(X, t(true.beta)))
y <- rpois(N, lambda.true)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    return(beta)
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    lambda <- exp(tcrossprod(Data$X, t(beta)))
    LL <- sum(dpois(Data$y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(length(lambda), lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0,J)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Zero-Inflated Poisson (ZIP)

The zero-inflated Poisson model handles count data with excess zeros by
mixing a point mass at zero (governed by the logistic sub-model on
\\\alpha\\) with a standard Poisson count process (governed by the
log-linear sub-model on \\\beta\\). The latent indicator \\z_i=1\\ when
observation \\i\\ is a structural zero, and \\z_i=0\\ when it comes from
the Poisson component. The zero-inflated Poisson model was introduced by
Lambert [\[69\]](#ref69) to handle excess zeros arising from a mixture
of structural and sampling zeros. ZIP models are applied in ecology for
modelling species count data where many sites have zero abundance
because the species is truly absent, not merely undetected.

#### Form

\\\textbf{y} \sim \mathcal{P}(\Lambda\_{1:N,2})\\ \\\textbf{z} \sim
\mathcal{BERN}(\Lambda\_{1:N,1})\\ \\\textbf{z}\_i = \left\\
\begin{array}{l l} 1 & \quad \mbox{if \$\textbf{y}\_i = 0\$}\\ 0 \\
\end{array} \right. \\ \\\Lambda\_{i,2} = \left\\ \begin{array}{l l} 0 &
\quad \mbox{if \$\Lambda\_{i,1} \ge 0.5\$}\\ \Lambda\_{i,2} \\
\end{array} \right. \\ \\\Lambda\_{1:N,1} = \frac{1}{1 +
\exp(-\textbf{X}\_1 \alpha)}\\ \\\Lambda\_{1:N,2} = \exp(\textbf{X}\_2
\beta)\\ \\\alpha_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J_1\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J_2\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the zero-inflation mechanism requires a conditional zeroing of
the Poisson rate based on the logistic component, mixing two sub-models
with distinct link functions and a deterministic switch, which falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42214)
N <- 500
J1 <- 4
J2 <- 3

#### True parameter values
true.alpha <- c(-0.5, 0.3, -0.2, 0.4)
true.beta  <- c(0.8, -0.3, 0.5)

#### Generate design matrices
X1 <- matrix(runif(N*J1,-2,2), N, J1); X1[,1] <- 1
X2 <- matrix(runif(N*J2,-2,2), N, J2); X2[,1] <- 1

#### Generate zero-inflated Poisson data
p <- invlogit(tcrossprod(X1, t(true.alpha)))
mu <- exp(tcrossprod(X2, t(true.beta)))
y <- ifelse(p > 0.5, 0, rpois(N, mu))
z <- ifelse(y == 0, 1, 0)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J1), beta=rep(0,J2)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J1)
    beta <- rnorm(Data$J2)
    return(c(alpha, beta))
    }
Data <- list(J1=J1, J2=J2, N=N, PGF=PGF, X1=X1, X2=X2,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.beta=pos.beta, y=y, z=z)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[Data$pos.alpha] <- alpha <- interval(parm[Data$pos.alpha], -5, 5)
    parm[Data$pos.beta] <- beta <- interval(parm[Data$pos.beta], -5, 5)
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 5, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 5, log=TRUE))
    ### Log-Likelihood
    Lambda <- matrix(NA, Data$N, 2)
    Lambda[,1] <- invlogit(tcrossprod(Data$X1, t(alpha)))
    Lambda[,2] <- exp(tcrossprod(Data$X2, t(beta))) + 1e-100
    Lambda[which(Lambda[,1] >= 0.5),2] <- 0
    LL <- sum(dbern(Data$z, Lambda[,1], log=TRUE),
        dpois(Data$y, Lambda[,2], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(nrow(Lambda), Lambda[,2]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- GIV(Model, Data, n=10000)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
post.alpha <- fit$Summary2[pos.alpha, "Mean"]
post.beta  <- fit$Summary2[pos.beta, "Mean"]
cat("True alpha:      ", true.alpha, "\n")
cat("Posterior mean:   ", round(post.alpha, 3), "\n")
cat("True beta:       ", true.beta, "\n")
cat("Posterior mean:   ", round(post.beta, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Generalized linear models: continuous, survival, and censored

### Beta Regression

Beta regression models a continuous response bounded in \\(0,1)\\ using
the Beta distribution reparameterized in terms of its mean \\\mu\\ and
precision \\\phi\\. The mean is linked to a linear predictor through the
probit (normal CDF) link, ensuring that \\\mu\\ lies strictly within the
unit interval. This approach is standard for modelling rates,
proportions, and other fractional data. Beta regression was introduced
by Ferrari and Cribari-Neto [\[40\]](#ref40) for modelling continuous
proportions bounded in (0,1) with heteroscedastic variance. In ecology,
beta regression models proportional vegetation cover as a function of
elevation, soil moisture, and disturbance history.

#### Form

\\\textbf{y} \sim \mathcal{BETA}(a,b)\\ \\a = \mu \phi\\ \\b = (1 - \mu)
\phi\\ \\\mu = \Phi(\beta_1 + \beta_2 \textbf{x}), \quad \mu \in (0,
1)\\ \\\beta_j \sim \mathcal{N}(0, 10), \quad j=1,\dots,J\\ \\\phi \sim
\mathcal{HC}(25)\\ where \\\Phi\\ is the normal CDF.

#### model_spec

The Beta distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry with `shape1` and `shape2` arguments. The mean-precision
reparameterization is expressed through deterministic statements, with
the probit link written as
[`pnorm()`](https://rdrr.io/r/stats/Normal.html).

``` r

spec <- model_spec("
  y ~ Beta(a, b)
  a = mu * phi
  b = (1 - mu) * phi
  mu = pnorm(beta1 + beta2 * x)
  beta1 ~ Normal(0, sqrt(10))
  beta2 ~ Normal(0, sqrt(10))
  phi ~ HalfCauchy(25)
")
```

Note that
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
uses `Beta(shape1, shape2)`, which maps to R’s `dbeta(shape1, shape2)`.
The deterministic statements for `a` and `b` perform the mean-precision
to shape reparameterization. Because this model uses `dnormv` (variance
parameterization) in the manual Model function with variance 10, the
equivalent
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
prior uses `Normal(0, sqrt(10))` (standard deviation parameterization).

#### Data

``` r

set.seed(42114)
N <- 200
J <- 2
true.beta <- c(0.3, 0.6)
true.phi <- 15
x <- rnorm(N)
mu.true <- pnorm(true.beta[1] + true.beta[2] * x)
a.true <- mu.true * true.phi
b.true <- (1 - mu.true) * true.phi
y <- rbeta(N, a.true, b.true)
mon.names <- "LP"
parm.names <- c("beta[1]","beta[2]","phi")
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
PGF <- function(Data) return(c(rnormv(2,0,10), rhalfcauchy(1,5)))
Data <- list(PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.beta=pos.beta, pos.phi=pos.phi, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 10, log=TRUE))
    phi.prior <- dhalfcauchy(phi, 25, log=TRUE)
    ### Log-Likelihood
    mu <- interval(pnorm(beta[1] + beta[2]*Data$x), 0.001, 0.999,
        reflect=FALSE)
    a <- mu * phi
    b <- (1 - mu) * phi
    LL <- sum(dbeta(Data$y, a, b, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + phi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbeta(length(mu), a, b), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,2), 0.01)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "phi"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]  <- true.phi[seq_along(pos.phi)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Beta-Binomial

The Beta-Binomial model places independent Beta priors on
per-observation success probabilities \\\pi_i\\, each governing a
Binomial likelihood for the observed count \\y_i\\ out of \\n_i\\
trials. This is the simplest hierarchical structure for overdispersed
binomial data, where each unit gets its own probability rather than
sharing a single global rate. The beta-binomial model accounts for
overdispersion in binomial data through a conjugate beta mixing
distribution, as described by Williams [\[115\]](#ref115) and Skellam
[\[99\]](#ref99). Beta-binomial models are applied in toxicology to
analyse litter-level variation in teratogenicity studies, where foetal
outcomes within litters are correlated.

#### Form

\\\textbf{y}\_i \sim \mathcal{BIN}(\textbf{n}\_i, \pi_i), \quad
i=1,\dots,N\\ \\\pi_i \sim \mathcal{BETA}(\alpha, \beta) \in
\[0.001,0.999\]\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires \\N\\ separate per-observation probability
parameters \\\pi_i\\, each with its own Beta prior. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
framework is designed for a fixed, low-dimensional parameter vector with
vectorized likelihood evaluation, not for observation-indexed latent
parameters. The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42115)
N <- 20
n <- round(runif(N, 50, 100))

#### True parameter values
true.pi <- rbeta(N, 2, 5)

#### Generate response from the true model
y <- rbinom(N, n, true.pi)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(pi=rep(0,N)))
PGF <- function(Data) {
    pi <- rbeta(Data$N, 1, 1)
    return(pi)
    }
Data <- list(N=N, PGF=PGF, mon.names=mon.names, n=n,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[1:Data$N] <- pi <- interval(parm[1:Data$N], 0.001, 0.999)
    ### Log-Prior
    pi.prior <- sum(dbeta(pi, 1, 1, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dbinom(Data$y, Data$n, pi, log=TRUE))
    ### Log-Posterior
    LP <- LL + pi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(Data$N, Data$n, pi), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0.5, N))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: per-observation probabilities
post.pi <- fit$Summary2[1:N, "Mean"]
cat("True pi (first 5):      ", round(true.pi[1:5], 3), "\n")
cat("Posterior mean (first 5):", round(post.pi[1:5], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("pi", parm.names)] <- true.pi[seq_along(grep("pi", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Exponential Survival Regression **\[NEW\]**

Exponential survival regression models time-to-event data under the
memoryless assumption by relating the hazard rate \\\lambda_i\\ to
covariates through a log-linear predictor \\\textbf{X}\_i \beta\\, so
that the survival time follows an exponential distribution whose rate
varies with subject characteristics. The exponential model is the
simplest parametric survival model and serves as a building block for
more flexible specifications such as Weibull or piecewise-exponential
formulations. Parametric survival regression was developed by Lawless
[\[128\]](#ref128), and the Bayesian treatment of survival models was
systematised by Ibrahim et al. [\[60\]](#ref60). Exponential survival
models are used in reliability engineering to estimate the failure rate
of electronic components as a function of operating conditions.

#### Form

\\t_i \sim \text{Exp}(\lambda_i), \quad i=1,\dots,N\\ \\\log(\lambda_i)
= \textbf{X}\_i \beta\\ \\\beta_j \sim \mathcal{N}(0, 100), \quad
j=1,\dots,J\\

#### model_spec() notation

The Exponential distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry with a `rate` argument. The log link maps the linear predictor
to the hazard rate through a deterministic
[`exp()`](https://rdrr.io/r/base/Log.html) node.

``` r

spec <- model_spec("
  y ~ Exponential(lambda)
  lambda = exp(X %*% beta)
  beta[j] ~ NormalV(0, 100), j = 1,...,J
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42212)
N <- 200
J <- 3   # Intercept + 2 covariates

#### True parameter values
true.beta <- c(-1.0, 0.5, -0.3)

#### Generate design matrix and survival times
x1 <- rnorm(N)
x2 <- rnorm(N)
X <- cbind(1, x1, x2)
lambda.true <- exp(X %*% true.beta)
y <- rexp(N, lambda.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J)))
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J, 0, 0.5)
    return(beta)
    }
Data <- list(N=N, J=J, X=X, y=y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 100, log=TRUE))
    ### Log-Likelihood
    lambda <- exp(tcrossprod(Data$X, t(beta)))
    LL <- sum(dexp(Data$y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rexp(Data$N, lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0, J)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Gamma Regression

Gamma regression models a strictly positive, right-skewed continuous
response through a log link relating the mean to a linear predictor
\\\textbf{X}\beta\\; the shape parameter \\\alpha\\ controls the
dispersion. It is appropriate for response variables such as durations,
costs, or concentrations where the variance increases with the mean.
Gamma regression is a generalised linear model for positive continuous
outcomes with constant coefficient of variation, formalised by Nelder
and Wedderburn [\[89\]](#ref89) and McCullagh and Nelder
[\[80\]](#ref80). Gamma regression is used in insurance for modelling
claim amounts, where larger claims exhibit proportionally larger
variance.

#### Form

\\\textbf{y} \sim \mathcal{G}(\lambda \tau, \tau)\\ \\\lambda =
\exp(\textbf{X} \beta)\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\tau \sim \mathcal{HC}(25)\\

#### model_spec() notation

The Gamma distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry with `shape` and `rate` arguments. The log-link
reparameterization maps the linear predictor to the mean through a
deterministic [`exp()`](https://rdrr.io/r/base/Log.html) node, and the
shape is expressed as `tau * lambda`.

``` r

spec <- model_spec("
  y ~ Gamma(shape, rate)
  shape = tau * lambda
  rate = tau
  lambda = exp(X %*% beta)
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  tau ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42131)
N <- 200
J <- 3

#### True parameter values
true.beta <- c(1.0, 0.5, -0.3)  # log-link coefficients (intercept + 2 predictors)
true.tau <- 4.0                   # precision (inverse dispersion)

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
lambda.true <- exp(X %*% true.beta)
y <- rgamma(N, shape = true.tau * lambda.true, rate = true.tau)

#### Assemble Data list
mon.names <- c("LP", "sigma2")
parm.names <- as.parm.names(list(beta=rep(0,J), tau=0))
pos.beta <- grep("beta", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    tau <- runif(1)
    return(c(beta, tau))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.tau=pos.tau, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    parm[Data$pos.tau] <- tau
    sigma2 <- 1/tau
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    ### Log-Likelihood
    lambda <- exp(tcrossprod(Data$X, t(beta)))
    LL <- sum(dgamma(Data$y, tau*lambda, tau, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma2),
        yhat=rgamma(nrow(lambda), tau*lambda, tau), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True tau:      ", true.tau, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.tau, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.tau]  <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Inverse Gaussian Regression

The inverse Gaussian (or Wald) distribution models positive,
right-skewed data and arises naturally as the first-passage time
distribution for Brownian motion with drift. This regression specifies
the mean through an exponential link \\\mu = \exp(\textbf{X}\beta) + C\\
with shape parameter \\\lambda\\ and half-Cauchy prior, making it
suitable for response times and reaction-time data where both the mean
and variance are positive. The inverse Gaussian distribution was
characterised by Tweedie [\[110\]](#ref110) and its regression
applications developed by Chhikara and Folks [\[23\]](#ref23). Inverse
Gaussian regression is applied in industrial engineering to model
failure times of components under accelerated life testing.

#### Form

\\\textbf{y} \sim \mathcal{IG}(\mu, \lambda)\\ \\\mu =
\exp(\textbf{X}\beta) + C\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\lambda \sim \mathcal{HC}(25)\\ where \\C\\ is a small
constant, such as 1.0E-10.

#### model_spec() notation

The inverse Gaussian distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry as `InvGaussian(mu, lambda)`. The exponential link with offset
constant is expressed through a deterministic node.

``` r

spec <- model_spec("
  y ~ InvGaussian(mu, lambda)
  mu = exp(X %*% beta) + 1e-10
  beta[j] ~ NormalV(0, 1000), j = 1,...,J
  lambda ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42133)
N <- 200
J <- 3

#### True parameter values
true.beta <- c(0.8, 0.4, -0.3)   # log-link coefficients
true.lambda <- 5.0                 # shape parameter

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- exp(X %*% true.beta) + 1e-10
y <- rinvgaussian(N, mu.true, true.lambda)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), lambda=0))
pos.beta <- grep("beta", parm.names)
pos.lambda <- grep("lambda", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    lambda <- runif(1)
    return(c(beta, lambda))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.lambda=pos.lambda, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    lambda.prior <- dhalfcauchy(lambda, 25, log=TRUE)
    ### Log-Likelihood
    mu <- exp(tcrossprod(Data$X, t(beta))) + 1.0E-10
    LL <- sum(dinvgaussian(Data$y, mu, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + lambda.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rinvgaussian(length(mu), mu, lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True lambda:   ", true.lambda, "\n")
cat("Posterior mean:", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Proportional Hazards Regression, Weibull

The Weibull proportional hazards model is a parametric survival model
that assumes the baseline hazard follows a Weibull distribution with
shape \\\alpha\\ and scale determined by the covariates. The likelihood
factorizes into contributions from censored and uncensored observations,
with the log-hazard linear in \\\textbf{X}\beta\\. Although the
dependent variable is usually denoted as \\\textbf{t}\\ in survival
analysis, it is denoted here as \\\textbf{y}\\ so Lucifer recognizes it
as a dependent variable for posterior predictive checks. The Weibull
proportional hazards model was placed in a Bayesian framework by Ibrahim
et al. [\[60\]](#ref60), extending the semiparametric model of Cox
[\[27\]](#ref27). Weibull hazard models are applied in engineering
reliability to estimate component lifetimes under different operating
conditions.

#### Form

\\\textbf{y}\_i \sim \mathcal{WEIB}(\gamma, \mu_i), \quad i=1,\dots,N\\
\\\mu = \exp(\textbf{X} \beta)\\ \\\beta_j \sim \mathcal{N}(0, 1000),
\quad j=1,\dots,J\\ \\\gamma \sim \mathcal{G}(1, 0.001)\\

#### model_spec

The Weibull proportional hazards model involves censoring-indicator
logic and a parameterization that does not map directly to a single
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
distribution entry. The manual Model function is therefore the
recommended approach.

#### Ground truth and data

``` r

set.seed(42161)
N <- 200
J <- 4
true.beta  <- c(1.5, 0.5, -0.3, 0.2)
true.gamma <- 1.8  # Weibull shape

#### Design matrix with intercept
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 0.5), N, J - 1))

#### Generate Weibull survival times
scale.true <- exp(as.vector(tcrossprod(X, t(true.beta))))
y <- rweibull(N, shape = true.gamma, scale = scale.true)

mon.names <- "LP"
parm.names <- as.parm.names(list(beta = rep(0, J), gamma = 0))
pos.beta  <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    beta  <- rnorm(Data$J)
    gamma <- rgamma(1, 1e-3)
    return(c(beta, gamma))
    }
Data <- list(J = J, N = N, PGF = PGF, X = X, mon.names = mon.names,
    parm.names = parm.names, pos.beta = pos.beta, pos.gamma = pos.gamma,
    y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.gamma] <- gamma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- dgamma(gamma, 1, 1.0E-3, log=TRUE)
    ### Log-Likelihood
    mu <- exp(tcrossprod(Data$X, t(beta)))
    LL <- sum(dweibull(Data$y, gamma, mu, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rweibull(length(mu), gamma, mu), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "gamma"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Tobit Regression (Censored) **\[NEW\]**

The Tobit model handles left-censored continuous responses by positing a
latent variable \\y^\*\_i = \textbf{X}\_i \beta + \epsilon_i\\ that is
observed only when positive; otherwise the response is recorded as zero.
The log-likelihood therefore has two components: a normal density for
uncensored observations and a normal CDF evaluated at the censoring
threshold for the censored ones. This piecewise structure makes Tobit a
standard tool for limited dependent variables. The Tobit model was
introduced by Tobin [\[141\]](#ref141) for analysing household
expenditure with a mass point at zero, and Bayesian estimation via data
augmentation was developed by Chib [\[142\]](#ref142). Tobit regression
is widely applied in consumer economics to model household spending
categories where zero expenditure is common.

#### Form

\\y^\*\_i = \textbf{X}\_i \beta + \epsilon_i, \quad \epsilon_i \sim
\mathcal{N}(0, \sigma^2)\\ \\y_i = \max(0, y^\*\_i)\\ \\\mathcal{L}(y_i
\mid \beta, \sigma) = \begin{cases} \Phi\left(\frac{-\textbf{X}\_i
\beta}{\sigma}\right) & \text{if } y_i = 0 \\ \phi\left(\frac{y_i -
\textbf{X}\_i \beta}{\sigma}\right) / \sigma & \text{if } y_i \> 0
\end{cases}\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\
\\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The censored likelihood with a piecewise log-probability depending on
whether \\y_i = 0\\ or \\y_i \> 0\\ falls outside the declarative scope
of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42231)
N <- 200
J <- 3

#### True parameter values
true.beta <- c(1.0, 1.5, -0.8)
true.sigma <- 1.0

#### Generate design matrix and latent response
X <- cbind(1, rnorm(N), rnorm(N))
y.star <- as.vector(X %*% true.beta) + rnorm(N, 0, true.sigma)
y <- pmax(0, y.star)
cat("Proportion censored:", round(mean(y == 0), 3), "\n")

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(N=N, J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma,
    y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(Data$X %*% beta)
    censored <- Data$y == 0
    LL <- sum(ifelse(censored,
        pnorm(0, mu, sigma, log.p=TRUE),
        dnorm(Data$y, mu, sigma, log=TRUE)))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=pmax(0, rnorm(Data$N, mu, sigma)), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:     ", true.beta, "\n")
cat("Post. mean:    ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Post. mean:    ", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Tweedie Regression **\[NEW\]**

The Tweedie distribution with power parameter \\p \in (1, 2)\\ is a
compound Poisson-gamma model that places a discrete mass at zero while
producing continuous positive values elsewhere, making it the natural
likelihood for zero-inflated positive-continuous responses. A Tweedie
GLM links the mean \\\mu_i\\ to the linear predictor through a log link
and estimates the dispersion \\\phi\\ and the power \\p\\ jointly with
the regression coefficients. The Tweedie family was characterised by
Jorgensen [\[145\]](#ref145) within the exponential dispersion model
framework, and practical density evaluation was developed by Dunn and
Smyth [\[146\]](#ref146) using a series expansion. Tweedie regression is
the standard model in actuarial science for aggregate insurance claims
where many policies produce zero claims and the positive claims are
right-skewed.

#### Form

\\y_i \sim \text{Tweedie}(\mu_i, \phi, p), \quad i=1,\dots,N\\
\\\log(\mu_i) = \textbf{X}\_i \beta\\ \\\beta_j \sim \mathcal{N}(0,
100), \quad j=1,\dots,J\\ \\\phi \sim \mathcal{HC}(25)\\ \\p \sim
\mathcal{U}(1.01, 1.99)\\

#### model_spec() notation

The Tweedie distribution is registered in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL and can be expressed directly.

``` r

spec <- model_spec("
  y ~ Tweedie(mu, phi, p)
  log(mu) = beta1 + beta2 * x1 + beta3 * x2
  beta1 ~ NormalV(0, 100)
  beta2 ~ NormalV(0, 100)
  beta3 ~ NormalV(0, 100)
  phi ~ HalfCauchy(25)
  p ~ Uniform(1.01, 1.99)
")
```

#### Ground truth and data

``` r

set.seed(42232)
N <- 200
J <- 3

#### True parameter values
true.beta <- c(1.0, 0.5, -0.3)
true.phi <- 1.5
true.p <- 1.5

#### Generate design matrix and response
X <- cbind(1, rnorm(N), rnorm(N))
mu.true <- exp(X %*% true.beta)
y <- rtweedie(N, as.vector(mu.true), true.phi, true.p)
cat("Proportion of zeros:", round(mean(y == 0), 3), "\n")

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), phi=0, p=0))
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.p <- grep("^p$", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J, 0, 0.5)
    phi <- runif(1, 0.5, 3)
    p <- runif(1, 1.01, 1.99)
    return(c(beta, phi, p))
    }
Data <- list(N=N, J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.phi=pos.phi,
    pos.p=pos.p, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
    parm[Data$pos.phi] <- phi
    p <- interval(parm[Data$pos.p], 1.01, 1.99)
    parm[Data$pos.p] <- p
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 100, log=TRUE))
    phi.prior <- dhalfcauchy(phi, 25, log=TRUE)
    p.prior <- dunif(p, 1.01, 1.99, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(exp(Data$X %*% beta))
    LL <- sum(dtweedie(Data$y, mu, phi, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + phi.prior + p.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rtweedie(Data$N, mu, phi, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1, 1.5)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:    ", true.beta, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True phi:     ", true.phi, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")
cat("True p:       ", true.p, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.p, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.p]    <- true.p[seq_along(pos.p)]
ground_truth[pos.phi]  <- true.phi[seq_along(pos.phi)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Ordinal and categorical models

### Contingency Table

The two-way contingency table, matrix \\\textbf{Y}\\, can easily be
extended to more dimensions. Contingency table \\\textbf{Y}\\ has J rows
and K columns. The cell counts are fit with Poisson regression,
according to intercept \\\alpha\\, main effects \\\beta_j\\ for each
row, main effects \\\gamma_k\\ for each column, and interaction effects
\\\delta\_{j,k}\\ for dependence effects. An omnibus (all cells) test of
independence is done by estimating two models (one with \\\delta\\, and
one without), and a large enough Bayes Factor indicates a violation of
independence when the model with \\\delta\\ fits better than the model
without \\\delta\\. In an ANOVA-like style, main effects contrasts can
be used to distinguish rows or groups of rows from each other, as well
as with columns. Likewise, interaction effects contrasts can be used to
test independence in groups of \\\delta\_{j,k}\\ elements. Finally,
single-cell interactions can be used to indicate violations of
independence for a given cell, such as when zero is not within its 95%
probability interval. Bayesian analysis of contingency tables was
developed by Albert [\[1\]](#ref1) using Dirichlet priors on cell
probabilities. Contingency table models are applied in social science to
test the association between categorical variables such as educational
attainment and voting preference.

#### Form

\\\textbf{Y}\_{j,k} \sim \mathcal{P}(\lambda\_{j,k}), \quad j=1,\dots,J,
\quad k=1,\dots,K\\ \\\lambda\_{j,k} = \exp(\alpha + \beta_j +
\gamma_k + \delta\_{j,k}), \quad j=1,\dots,J, \quad k=1,\dots,K\\
\\\alpha \sim \mathcal{N}(0, 1000)\\ \\\beta_j \sim \mathcal{N}(0,
\beta^2\_\sigma), \quad j=1,\dots,J\\ \\\beta_J = -
\displaystyle\sum^{J-1}\_{j=1} \beta_j\\ \\\beta\_\sigma \sim
\mathcal{HC}(25)\\ \\\gamma_k \sim \mathcal{N}(0, \gamma^2\_\sigma),
\quad k=1,\dots,K\\ \\\gamma_K = - \displaystyle\sum^{K-1}\_{k=1}
\gamma_k\\ \\\gamma\_\sigma \sim \mathcal{HC}(25)\\ \\\delta\_{j,k} \sim
\mathcal{N}(0, \delta^2\_\sigma)\\ \\\delta\_{J,K} = - \displaystyle\sum
\delta\_{-J,-K}\\ \\\delta\_\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires sum-to-zero constraints on multiple parameter
vectors (\\\beta_J = -\sum \beta\_{-J}\\, \\\gamma_K = -\sum
\gamma\_{-K}\\, \\\delta\_{J,K} = -\sum \delta\_{-J,-K}\\). These
constraints derive the last element of each effect vector from the
sampled parameters through imperative operations that fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42121)
J <- 4  # Rows
K <- 4  # Columns

#### True parameter values
true.alpha <- 3.0
true.beta.free <- c(0.5, -0.3, 0.1)
true.beta <- c(true.beta.free, -sum(true.beta.free))
true.gamma.free <- c(-0.2, 0.4, -0.1)
true.gamma <- c(true.gamma.free, -sum(true.gamma.free))
true.delta.free <- c(0.3, -0.2, 0.1, -0.4, 0.2, -0.1,
    0.0, 0.3, -0.2, 0.1, -0.3, 0.2, -0.1, 0.0, 0.1)
true.delta <- c(true.delta.free, -sum(true.delta.free))
true.delta.mat <- matrix(true.delta, J, K)

#### Generate contingency table from Poisson model
beta.mat <- matrix(true.beta, J, K)
gamma.mat <- matrix(true.gamma, J, K, byrow=TRUE)
lambda.true <- exp(true.alpha + beta.mat + gamma.mat + true.delta.mat)
Y <- matrix(rpois(J*K, lambda.true), J, K)
rownames(Y) <- c("Black", "Blond", "Brunette", "Red")
colnames(Y) <- c("Blue", "Brown", "Green", "Hazel")

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,J-1),
    gamma=rep(0,K-1), delta=rep(0,J*K-1), b.sigma=0, g.sigma=0,
    d.sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.b.sigma <- grep("b.sigma", parm.names)
pos.g.sigma <- grep("g.sigma", parm.names)
pos.d.sigma <- grep("d.sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1, log(mean(Y)), 1)
    beta <- rnorm(Data$J-1)
    gamma <- rnorm(Data$K-1)
    delta <- rnorm(Data$J*Data$K-1)
    sigma <- runif(3)
    return(c(alpha, beta, gamma, delta, sigma))
    }
Data <- list(J=J, K=K, PGF=PGF, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.gamma=pos.gamma, pos.delta=pos.delta, pos.b.sigma=pos.b.sigma,
    pos.g.sigma=pos.g.sigma, pos.d.sigma=pos.d.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    beta.sigma <- interval(parm[Data$pos.b.sigma], 1e-100, Inf)
    parm[Data$pos.b.sigma] <- beta.sigma
    gamma.sigma <- interval(parm[Data$pos.g.sigma], 1e-100, Inf)
    parm[Data$pos.g.sigma] <- gamma.sigma
    delta.sigma <- interval(parm[Data$pos.d.sigma], 1e-100, Inf)
    parm[Data$pos.d.sigma] <- delta.sigma
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    beta <- c(beta, -sum(beta))
    gamma <- parm[Data$pos.gamma]
    gamma <- c(gamma, -sum(gamma))
    delta <- parm[Data$pos.delta]
    delta <- c(delta, -sum(delta))
    delta <- matrix(delta, Data$J, Data$K)
    ### Log-Hyperprior
    beta.sigma.prior <- dhalfcauchy(beta.sigma, 25, log=TRUE)
    gamma.sigma.prior <- dhalfcauchy(gamma.sigma, 25, log=TRUE)
    delta.sigma.prior <- dhalfcauchy(delta.sigma, 25, log=TRUE)
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnorm(beta, 0, beta.sigma, log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, gamma.sigma, log=TRUE))
    delta.prior <- sum(dnorm(delta, 0, delta.sigma, log=TRUE))
    ### Log-Likelihood
    beta <- matrix(beta, Data$J, Data$K)
    gamma <- matrix(gamma, Data$J, Data$K, byrow=TRUE)
    lambda <- exp(alpha + beta + gamma + delta)
    LL <- sum(dpois(Data$Y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + beta.sigma.prior +
        gamma.prior + gamma.sigma.prior + delta.prior +
        delta.sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(length(lambda), lambda),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(log(mean(Y)), rep(0, J-1), rep(0, K-1), rep(0, J*K-1),
    rep(1, 3))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: intercept
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")

### Parameter recovery: row effects
post.beta <- fit$Summary2[pos.beta, "Mean"]
cat("True beta (free):", round(true.beta.free, 3), "\n")
cat("Posterior mean:  ", round(post.beta, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multinomial Logit

The multinomial logit is the standard baseline-category model for
unordered polytomous response variables. Each category \\j \< J\\ has
its own coefficient vector \\\beta_j\\, with the last category \\J\\
serving as the reference (utility fixed at zero). Choice probabilities
are computed via the softmax transformation of the linear predictors.
The multinomial logit model was introduced by McFadden [\[82\]](#ref82)
within the random utility maximisation framework. Multinomial logit is
used in land cover classification to predict vegetation types from
satellite-derived spectral indices.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires per-category coefficient vectors \\\beta_j\\
with a structured softmax likelihood, where the reference category is
fixed at zero utility. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
framework does not support this kind of multi-equation categorical
choice specification. The manual `Model` function below handles this
directly.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(\textbf{p}\_{i,1:J}), \quad
i=1,\dots,N\\ \\\textbf{p}\_{i,j} = \frac{\phi\_{i,j}}{\sum^J\_{j=1}
\phi\_{i,j}}, \quad \sum^J\_{j=1} \textbf{p}\_{i,j} = 1\\ \\\phi =
\exp(\mu)\\ \\\mu\_{i,J} = 0, \quad i=1,\dots,N\\ \\\mu\_{i,j} =
\textbf{X}\_{i,1:K} \beta\_{j,1:K} \in \[-700,700\], \quad
j=1,\dots,(J-1)\\ \\\beta\_{j,k} \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,(J-1), \quad k=1,\dots,K\\

#### Ground truth and data

``` r

set.seed(42178)
N <- 200
J <- 3   # Number of categories
K <- 3   # Number of predictors (including intercept)

#### True parameter values
true.beta <- matrix(c(0.5, -0.3, 0.8, -0.4, 0.6, -0.2), J-1, K)

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
for (j in 2:K) X[, j] <- CenterScale(X[, j])
mu.true <- cbind(tcrossprod(X, true.beta), 0)
phi.true <- exp(mu.true)
p.true <- phi.true / rowSums(phi.true)
y <- rcat(N, p.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,J-1,K)))
PGF <- function(Data) {
    beta <- rnorm((Data$J-1)*Data$K)
    return(beta)
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm, Data$J-1, Data$K)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(0, Data$N, Data$J)
    mu[,-Data$J] <- tcrossprod(Data$X, beta)
    mu <- interval(mu, -700, 700, reflect=FALSE)
    phi <- exp(mu)
    p <- phi / rowSums(phi)
    LL <- sum(dcat(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(p), p),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, (J-1)*K))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
post.beta <- matrix(colMeans(fit$Posterior2), J-1, K)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("beta", parm.names)] <- true.beta[seq_along(grep("beta", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multinomial Logit, Nested

The nested logit relaxes the independence of irrelevant alternatives
(IIA) property of the standard multinomial logit by grouping
alternatives into nests. The nesting parameter \\\alpha\\ (bounded
between 0 and 2) controls the degree of correlation among alternatives
within the nest; when \\\alpha = 1\\, the model collapses to the
standard multinomial logit. Here the first alternative is isolated
(non-nested), while the second and third are nested together. The nested
logit model was developed by Ben-Akiva [\[7\]](#ref7) and McFadden
[\[83\]](#ref83) to relax the IIA property of standard multinomial
logit. Nested logit is applied in transportation planning where
alternatives naturally group into nests, such as public versus private
transport modes.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the likelihood involves a nested logit structure with an
inclusive value term \\\textbf{I} = \log(\textbf{V})\\, the nesting
parameter \\\alpha\\ modulating both the isolated and nested branch
probabilities, and a product \\\iota = \alpha \beta_1\\ linking the
nesting parameter to regression coefficients. This hierarchical choice
structure falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(\textbf{P}\_{i,1:J}), \quad
i=1,\dots,N\\ \\\textbf{P}\_{1:N,1} = \frac{\textbf{R}}{\textbf{R} +
\exp(\alpha \textbf{I})}\\ \\\textbf{P}\_{1:N,2} = \frac{(1 -
\textbf{P}\_{1:N,1}) \textbf{S}\_{1:N,1}}{\textbf{V}}\\
\\\textbf{P}\_{1:N,3} = \frac{(1 - \textbf{P}\_{1:N,1})
\textbf{S}\_{1:N,2}}{\textbf{V}}\\ \\\textbf{R}\_{1:N} =
\exp(\mu\_{1:N,1})\\ \\\textbf{S}\_{1:N,1:2} = \exp(\mu\_{1:N,2:3})\\
\\\textbf{I} = \log(\textbf{V})\\ \\\textbf{V}\_i =
\displaystyle\sum^K\_{k=1} \textbf{S}\_{i,k}, \quad i=1,\dots,N\\
\\\mu\_{1:N,1} = \textbf{X} \iota \in \[-700,700\]\\ \\\mu\_{1:N,2} =
\textbf{X} \beta\_{2,1:K} \in \[-700,700\]\\ \\\iota = \alpha
\beta\_{1,1:K}\\ \\\alpha \sim \mathcal{EXP}(1) \in \[0,2\]\\
\\\beta\_{j,k} \sim \mathcal{N}(0, 1000), \quad j=1,\dots,(J-1) \quad
k=1,\dots,K\\ where there are \\J=3\\ categories of \\\textbf{y}\\,
\\K=3\\ predictors, \\\textbf{R}\\ is the non-nested alternative,
\\\textbf{S}\\ is the nested alternative, \\\textbf{V}\\ is the observed
utility in the nest, \\\alpha\\ is effectively 1 - correlation and has a
truncated exponential distribution, and \\\iota\\ is a vector of
regression effects for the isolated alternative after \\\alpha\\ is
taken into account. The third alternative is the reference category.

#### Ground truth and data

``` r

set.seed(42179)
N <- 200
J <- 3   # Number of categories
K <- 3   # Number of predictors (including intercept)

#### True parameter values
true.alpha <- 0.7
true.beta <- matrix(c(0.4, -0.3, 0.5, -0.2, 0.6, -0.4), J-1, K)

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
for (j in 2:K) X[, j] <- CenterScale(X[, j])

#### Generate response from the true nested logit model
iota.true <- true.alpha * true.beta[1, ]
mu.true <- matrix(0, N, J)
mu.true[, 1] <- tcrossprod(X, t(iota.true))
mu.true[, 2] <- tcrossprod(X, t(true.beta[2, ]))
R <- exp(mu.true[, 1])
S <- exp(mu.true[, -1])
V <- rowSums(S)
I <- log(V)
P <- matrix(0, N, J)
P[, 1] <- R / (R + exp(true.alpha * I))
P[, 2] <- (1 - P[, 1]) * S[, 1] / V
P[, 3] <- (1 - P[, 1]) * S[, 2] / V
y <- rcat(N, P)

#### Assemble Data list
mon.names <- c("LP", as.parm.names(list(iota=rep(0,K))))
parm.names <- as.parm.names(list(alpha=0, beta=matrix(0,J-1,K)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    alpha <- rtrunc(1, "exp", a=0, b=2, rate=1)
    beta <- rnorm((Data$J-1)*Data$K)
    return(c(alpha, beta))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    alpha.rate <- 1
    ### Parameters
    parm[Data$pos.alpha] <- alpha <- interval(parm[Data$pos.alpha],0,2)
    beta <- matrix(parm[Data$pos.beta], Data$J-1, Data$K)
    ### Log-Prior
    alpha.prior <- dtrunc(alpha, "exp", a=0, b=2, rate=alpha.rate,
        log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- P <- matrix(0, Data$N, Data$J)
    iota <- alpha * beta[1,]
    mu[,1] <- tcrossprod(Data$X, t(iota))
    mu[,2] <- tcrossprod(Data$X, t(beta[2,]))
    mu <- interval(mu, -700, 700, reflect=FALSE)
    R <- exp(mu[,1])
    S <- exp(mu[,-1])
    V <- rowSums(S)
    I <- log(V)
    P[,1] <- R / (R + exp(alpha*I))
    P[,2] <- (1 - P[,1]) * S[,1] / V
    P[,3] <- (1 - P[,1]) * S[,2] / V
    LL <- sum(dcat(Data$y, P, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,iota),
        yhat=rcat(nrow(P), P), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0.5, rep(0.1, (J-1)*K))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], J-1, K)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multinomial Probit

The multinomial probit model handles unordered categorical choice data
by introducing latent utility variables \\W\_{i,j}\\ for each
alternative, with multivariate normal errors whose covariance \\\Sigma\\
captures substitution patterns between alternatives. The chosen category
corresponds to the alternative with the highest latent utility,
implemented through sign constraints on \\W\\. The multinomial probit
model was developed for Bayesian estimation by McCulloch and Rossi
[\[81\]](#ref81) using data augmentation. MNP is used in marketing for
analysing consumer choice among multiple brands where the IIA assumption
of logit is unrealistic.

#### Form

\\\textbf{W}\_{i,1:(J-1)} \sim \mathcal{N}\_{J-1}(\mu\_{i,1:(J-1)},
\Sigma), \quad i=1,\dots,N\\ \\\textbf{W}\_{i,j} \in \left\\
\begin{array}{l l} \$\[0,10\]\$ & \quad \mbox{if \$\textbf{y}\_i =
j\$}\\ \$\[-10,0\]\$ \\ \end{array} \right. \\ \\\mu\_{1:N,j} =
\textbf{X} \beta\_{j,1:K}\\ \\\Sigma = \textbf{U}^T \textbf{U}\\
\\\beta\_{j,k} \sim \mathcal{N}(0, 10), \quad j=1,\dots,(J-1), \quad
k=1,\dots,K\\ \\\textbf{U}\_{j,k} \sim \mathcal{N}(0,1), \quad
j=1,\dots,(J-1), \quad k=1,\dots,(J-1), \quad j \ge k, \quad j \ne k =
1\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it requires latent utility variables \\W\\ with sign constraints
that depend on the observed choice, a Cholesky-parameterized covariance
matrix, and data augmentation logic, all of which fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42180)
N <- 200
J <- 3
K <- 3

#### Simulate multinomial probit data
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
for (j in 2:K) X[,j] <- CenterScale(X[,j])
true.beta <- matrix(c(0.5, 0.8, -0.3,
                       -0.4, 0.2, 0.6), J - 1, K)
mu <- tcrossprod(X, true.beta)
Sigma.true <- matrix(c(1.0, 0.3, 0.3, 1.0), J - 1, J - 1)
W <- mu + rmvn(N, rep(0, J - 1), Sigma.true)
y <- max.col(cbind(W, 0))  # J-th category is reference

#### Assemble Data list
S <- diag(J-1)
U <- matrix(NA,J-1,J-1)
U[upper.tri(U, diag=TRUE)] <- 0
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,(J-1),K),
    U=U, W=matrix(0,N,J-1)))
parm.names <- parm.names[-which(parm.names == "U[1,1]")]
pos.beta <- grep("beta", parm.names)
pos.U <- grep("U", parm.names)
pos.W <- grep("W", parm.names)
PGF <- function(Data) {
    beta <- rnorm((Data$J-1)*Data$K)
    U <- rnorm((Data$J-2) + (factorial(Data$J-1) /
        (factorial(Data$J-1-2)*factorial(2))))
    W <- matrix(runif(Data$N*(Data$J-1),-10,0), Data$N, Data$J-1)
    Y <- as.indicator.matrix(Data$y)
    W <- ifelse(Y[,-Data$J] == 1, abs(W), W)
    return(c(beta, U, as.vector(W)))}
Data <- list(J=J, K=K, N=N, PGF=PGF, S=S, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.U=pos.U, pos.W=pos.W,
    y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$J-1, Data$K)
    u <- c(0, parm[Data$pos.U])
    U <- diag(Data$J-1)
    U[upper.tri(U, diag=TRUE)] <- u
    diag(U) <- exp(diag(U))
    Sigma <- t(U) %*% U
    Sigma[1,] <- Sigma[,1] <- U[1,]
    W <- matrix(parm[Data$pos.W], Data$N, Data$J-1)
    Y <- as.indicator.matrix(Data$y)
    temp <- which(Y[,-c(Data$J)] == 1)
    W[temp] <- interval(W[temp], 0, 10)
    temp <- which(Y[,-c(Data$J)] == 0)
    W[temp] <- interval(W[temp], -10, 0)
    parm[Data$pos.W] <- as.vector(W)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 10, log=TRUE))
    U.prior <- sum(dnorm(u[-1], 0, 1, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, beta)
    #eta <- exp(cbind(mu,0))
    #p <- eta / rowSums(eta)
    LL <- sum(dmvn(W, mu, Sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=max.col(cbind(rmvn(nrow(mu), mu, Sigma),0)), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- GIV(Model, Data, PGF=TRUE)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multinomial-Dirichlet **\[NEW\]**

The Dirichlet-Multinomial model handles count compositions observed
across multiple groups by integrating out the group-specific category
probabilities analytically. Given \\K\\ categories with counts summing
to a fixed trial size \\n_i\\ within each of \\N\\ groups, the Dirichlet
prior on the probability simplex is conjugate to the multinomial
likelihood, yielding the Polya (Dirichlet-Multinomial) compound
distribution. The concentration parameters \\\alpha_k\\ govern how much
probability mass concentrates on each category across groups. The
Dirichlet-Multinomial distribution was developed by Mosimann
[\[144\]](#ref144) as a model for compositional data with
overdispersion. Dirichlet-Multinomial models are applied in species
abundance analysis to estimate community composition from replicated
quadrat counts.

#### Form

\\\textbf{y}\_i \sim \text{DirMult}(n_i, \boldsymbol{\alpha}), \quad
i=1,\dots,N\\ \\\log p(\textbf{y}\_i \mid \boldsymbol{\alpha}) =
\log\Gamma(\sum_k \alpha_k) - \log\Gamma(n_i + \sum_k \alpha_k) +
\sum\_{k=1}^{K} \left\[\log\Gamma(y\_{ik} + \alpha_k) -
\log\Gamma(\alpha_k)\right\]\\ \\\alpha_k \sim \mathcal{G}(1, 1), \quad
k=1,\dots,K\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Dirichlet-Multinomial compound distribution, which
integrates out the group-level probability simplexes analytically, is
not available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry. The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42221)
N <- 20
K <- 4
n.trials <- 50
true.alpha <- c(5, 3, 2, 1)

#### Simulate compositional count data
Y <- matrix(0, N, K)
for (i in 1:N) {
    p.i <- rdirichlet(1, true.alpha)
    Y[i, ] <- rmultinom(1, n.trials, p.i)
    }

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0, K)))
pos.alpha <- grep("alpha", parm.names)
n.vec <- rep(n.trials, N)
PGF <- function(Data) {
    return(rgamma(Data$K, 1, 1))
    }
Data <- list(K=K, N=N, PGF=PGF, Y=Y, mon.names=mon.names,
    n.vec=n.vec, parm.names=parm.names, pos.alpha=pos.alpha)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- interval(parm[Data$pos.alpha], 1e-100, Inf)
    parm[Data$pos.alpha] <- alpha
    ### Log-Prior
    alpha.prior <- sum(dgamma(alpha, 1, 1, log=TRUE))
    ### Log-Likelihood (Dirichlet-Multinomial / Polya)
    alpha.sum <- sum(alpha)
    LL <- 0
    for (i in 1:Data$N) {
        LL <- LL + lgamma(alpha.sum) - lgamma(Data$n.vec[i] + alpha.sum) +
            sum(lgamma(Data$Y[i, ] + alpha) - lgamma(alpha))
        }
    ### Log-Posterior
    LP <- LL + alpha.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=apply(matrix(rdirichlet(Data$N, alpha), Data$N, Data$K),
            1, function(p) rmultinom(1, Data$n.vec[1], p)[1, ]),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(1, K)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Ordinal Logit

Ordinal logistic regression models an ordered categorical response \\y
\in \\1,\dots,J\\\\ through cumulative probabilities linked to a linear
predictor via the logistic function. The cutpoint parameters
\\\delta_j\\ partition the latent continuous scale into ordered
categories. The cumulative logit model for ordinal data was introduced
by McCullagh [\[79\]](#ref79) using proportional odds, with Bayesian
treatment following Albert and Chib [\[2\]](#ref2). Ordinal logit models
are used in pain research to analyse patient-reported pain severity on
Likert scales (none, mild, moderate, severe).

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
due to its complex cumulative link structure with sorted cutpoints and
matrix-level probability construction.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(P\_{i,1:J})\\ \\P\_{,J} = 1 -
Q\_{,(J-1)}\\ \\P\_{,j} = \|Q\_{,j} - Q\_{,(j-1)}\|, \quad
j=2,\dots,(J-1)\\ \\P\_{,1} = Q\_{,1}\\ \\Q = \frac{1}{1 + \exp(\mu)}\\
\\\mu\_{,j} = \delta_j - \textbf{X} \beta \in \[-5,5\]\\ \\\beta_k \sim
\mathcal{N}(0, 1000), \quad k=1,\dots,K\\ \\\delta_j \sim
\mathcal{N}(0, 1) \in \[(j-1),j\], \quad j=1,\dots,(J-1)\\

#### Ground truth and data

``` r

set.seed(42155)
N <- 400
J <- 3   # Number of ordinal categories
K <- 3   # Number of predictors
true.beta  <- c(0.8, -0.5, 0.6)
true.delta <- c(-0.5, 0.8)  # Sorted cutpoints (J-1 = 2)
X <- matrix(rnorm(N * K), N, K)
for (k in 1:K) X[,k] <- (X[,k] - mean(X[,k])) / sd(X[,k])
linpred <- tcrossprod(X, t(true.beta))
# Cumulative probabilities via logistic CDF
Q <- matrix(NA, N, J - 1)
for (j in 1:(J - 1)) {
    Q[,j] <- 1 / (1 + exp(-(true.delta[j] - linpred)))
}
# Category probabilities
P <- matrix(NA, N, J)
P[,1] <- Q[,1]
for (j in 2:(J - 1)) P[,j] <- Q[,j] - Q[,j-1]
P[,J] <- 1 - Q[,J-1]
P <- pmax(P, 1e-10)
P <- P / rowSums(P)
# Simulate ordinal responses
y <- integer(N)
for (i in 1:N) y[i] <- sample(1:J, 1, prob = P[i,])
table(y)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,K), delta=rep(0,J-1)))
pos.beta <- grep("beta", parm.names)
pos.delta <- grep("delta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K)
    delta <- sort(rnorm(Data$J-1))
    return(c(beta, delta))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.delta=pos.delta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    delta <- interval(parm[Data$pos.delta], -5, 5)
    delta <- sort(delta)
    parm[Data$pos.delta] <- delta
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    delta.prior <- sum(dtrunc(delta, "norm", a=-5, b=5, log=TRUE,
        mean=0, sd=1))
    ### Log-Likelihood
    mu <- matrix(delta, Data$N, Data$J-1, byrow=TRUE) -
        matrix(tcrossprod(Data$X, t(beta)), Data$N, Data$J-1)
    P <- Q <- invlogit(mu)
    P[,-1] <- abs(Q[,-1] - Q[,-(Data$J-1)])
    P <- cbind(P, 1 - Q[,(Data$J-1)])
    LL <- sum(dcat(Data$y, P, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + delta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(P), P),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K), seq(from=-1, to=1, len=(J-1)))
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Ordinal Probit

Ordinal probit regression replaces the logistic CDF of the ordinal logit
with the standard normal CDF \\\Phi\\, which is the natural choice when
the latent variable generating the ordered categories is assumed
normally distributed. The ordinal probit model uses threshold parameters
on a latent normal variable, with Bayesian estimation via data
augmentation following Albert and Chib [\[2\]](#ref2). Ordinal probit is
applied in credit rating analysis to model the ordered default risk
categories assigned by rating agencies.

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
due to its complex cumulative link structure with sorted cutpoints and
matrix-level probability construction.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(P\_{i,1:J})\\ \\P\_{,J} = 1 -
Q\_{,(J-1)}\\ \\P\_{,j} = \|Q\_{,j} - Q\_{,(j-1)}\|, \quad
j=2,\dots,(J-1)\\ \\P\_{,1} = Q\_{,1}\\ \\Q = \phi(\mu)\\ \\\mu\_{,j} =
\delta_j - \textbf{X} \beta \in \[-5,5\]\\ \\\beta_k \sim \mathcal{N}(0,
1000), \quad k=1,\dots,K\\ \\\delta_j \sim \mathcal{N}(0, 1) \in
\[(j-1),j\], \quad j=1,\dots,(J-1)\\

#### Ground truth and data

``` r

set.seed(42156)
N <- 400
J <- 3   # Number of ordinal categories
K <- 3   # Number of predictors
true.beta  <- c(0.6, -0.4, 0.7)
true.delta <- c(-0.3, 0.9)  # Sorted cutpoints (J-1 = 2)
X <- matrix(rnorm(N * K), N, K)
for (k in 1:K) X[,k] <- (X[,k] - mean(X[,k])) / sd(X[,k])
linpred <- tcrossprod(X, t(true.beta))
# Cumulative probabilities via normal CDF
Q <- matrix(NA, N, J - 1)
for (j in 1:(J - 1)) {
    Q[,j] <- pnorm(true.delta[j] - linpred)
}
# Category probabilities
P <- matrix(NA, N, J)
P[,1] <- Q[,1]
for (j in 2:(J - 1)) P[,j] <- Q[,j] - Q[,j-1]
P[,J] <- 1 - Q[,J-1]
P <- pmax(P, 1e-10)
P <- P / rowSums(P)
# Simulate ordinal responses
y <- integer(N)
for (i in 1:N) y[i] <- sample(1:J, 1, prob = P[i,])
table(y)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,K), delta=rep(0,J-1)))
pos.beta <- grep("beta", parm.names)
pos.delta <- grep("delta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K)
    delta <- sort(rnorm(Data$J-1))
    return(c(beta, delta))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.delta=pos.delta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    delta <- interval(parm[Data$pos.delta], -5, 5)
    delta <- sort(delta)
    parm[Data$pos.delta] <- delta
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    delta.prior <- sum(dtrunc(delta, "norm", a=-5, b=5, log=TRUE,
        mean=0, sd=1))
    ### Log-Likelihood
    mu <- matrix(delta, Data$N, Data$J-1, byrow=TRUE) -
        matrix(tcrossprod(Data$X, t(beta)), Data$N, Data$J-1)
    P <- Q <- pnorm(mu)
    P[,-1] <- abs(Q[,-1] - Q[,-(Data$J-1)])
    P <- cbind(P, 1 - Q[,(Data$J-1)])
    LL <- sum(dcat(Data$y, P, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + delta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(P), P),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K), seq(from=-1, to=1, len=(J-1)))
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = "beta")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Discrete choice models

### Discrete Choice, Conditional Logit

The conditional logit model extends multinomial logit by incorporating
both individual-specific attributes \\\textbf{X}\\ (with
category-varying coefficients \\\beta\\) and choice-specific attributes
\\\textbf{Z}\\ (with generic coefficients \\\gamma\\). The last category
\\J\\ serves as the reference level, so \\\beta\\ coefficients are
estimated for categories \\1,\dots,J-1\\ only. The conditional logit
model for discrete choice was introduced by McFadden [\[82\]](#ref82),
who received the Nobel Prize in Economics for this work. Conditional
logit is the standard model in transportation research for predicting
mode choice (car, bus, train) as a function of travel time and cost.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the utility function combines individual attributes with
category-varying coefficients and choice-specific attributes with
generic coefficients, requiring per-category linear predictors with
different coefficient structures. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
framework does not support this kind of structured multinomial utility
specification. The manual `Model` function below handles this directly.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(\textbf{p}\_{i,1:J}), \quad
i=1,\dots,N, \quad j=1,\dots,J\\ \\\textbf{p}\_{i,j} =
\frac{\phi\_{i,j}}{\sum^J\_{j=1} \phi\_{i,j}}\\ \\\phi = \exp(\mu)\\
\\\mu\_{i,j} = \beta\_{j,1:K} \textbf{X}\_{i,1:K} + \gamma
\textbf{Z}\_{i,1:C} \in \[-700,700\], \quad j=1,\dots,(J-1)\\
\\\mu\_{i,J} = \gamma \textbf{Z}\_{i,1:C}\\ \\\beta\_{j,k} \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,(J-1)\\ \\\gamma_c \sim
\mathcal{N}(0, 1000)\\

#### Ground truth and data

``` r

set.seed(42122)
N <- 200
J <- 3   # Number of categories
K <- 3   # Number of individual attributes (including intercept)
C <- 4   # Number of choice-based attributes

#### True parameter values
true.beta <- matrix(c(0.5, -0.3, 0.8, -0.2, 0.4, -0.6), J-1, K)
true.gamma <- c(0.3, -0.5, 0.2, -0.1)

#### Generate design matrices
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
Z <- matrix(rnorm(N * C), N, C)
for (j in 2:K) X[, j] <- CenterScale(X[, j])
for (j in 1:C) Z[, j] <- CenterScale(Z[, j])

#### Generate response from the true model
mu.true <- matrix(tcrossprod(true.gamma, Z), N, J)
mu.true[, -J] <- mu.true[, -J] + tcrossprod(X, true.beta)
phi.true <- exp(mu.true)
p.true <- phi.true / rowSums(phi.true)
y <- rcat(N, p.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,J-1,K), gamma=rep(0,C)))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    beta <- rnorm((Data$J-1)*Data$K)
    gamma <- rnorm(Data$C)
    return(c(beta, gamma))
    }
Data <- list(C=C, J=J, K=K, N=N, PGF=PGF, X=X, Z=Z, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$J-1, Data$K)
    gamma <- parm[Data$pos.gamma]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(tcrossprod(gamma, Data$Z), Data$N, Data$J)
    mu[,-Data$J] <- mu[,-Data$J] + tcrossprod(Data$X, beta)
    mu <- interval(mu, -700, 700, reflect=FALSE)
    phi <- exp(mu)
    p <- phi / rowSums(phi)
    LL <- sum(dcat(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(p), p),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, (J-1)*K), rep(0, C))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], J-1, K)
post.gamma <- fit$Summary2[pos.gamma, "Mean"]
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))
cat("True gamma:     ", round(true.gamma, 3), "\n")
cat("Posterior mean: ", round(post.gamma, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Discrete Choice, Mixed Logit

The mixed logit extends the conditional logit by allowing
individual-level random coefficients \\\beta\_{j,k,i}\\ that vary across
observations according to normal hyperpriors with unknown means
\\\zeta^\mu\\ and standard deviations \\\zeta^\sigma\\. This captures
preference heterogeneity across decision-makers, which the
fixed-coefficient conditional logit cannot accommodate. The mixed logit
model generalises the conditional logit by allowing random taste
variation, as described by Train [\[109\]](#ref109). Mixed logit models
are used in health economics to estimate willingness-to-pay for
different attributes of health care programmes.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a three-dimensional array of random coefficients
\\\beta\_{j,k,i}\\ indexed by category, predictor, and observation, with
hierarchical normal hyperpriors on the mean and variance. This
high-dimensional random-coefficient structure falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(\textbf{p}\_{i,1:J}), \quad
i=1,\dots,N\\ \\\textbf{p}\_{i,j} = \frac{\phi\_{i,j}}{\sum^J\_{j=1}
\phi\_{i,j}}\\ \\\phi = \exp(\mu)\\ \\\mu\_{i,j} = \beta\_{j,1:K,i}
\textbf{X}\_{i,1:K} + \gamma \textbf{Z}\_{i,1:C} \in \[-700,700\], \quad
i=1,\dots,N, \quad j=1,\dots,(J-1)\\ \\\mu\_{i,J} = \gamma
\textbf{Z}\_{i,1:C}\\ \\\beta\_{j,k,i} \sim
\mathcal{N}(\zeta^\mu\_{j,k}, \zeta^\sigma2\_{j,k}), \quad i=1,\dots,N,
\quad j=1,\dots,(J-1), \quad k=1,\dots,K\\ \\\gamma_c \sim
\mathcal{N}(0, 1000), \quad c=1,\dots,C\\ \\\zeta^\mu\_{j,k} \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,(J-1), \quad k=1,\dots,K\\
\\\zeta^\sigma\_{j,k} \sim \mathcal{HC}(25), \quad j=1,\dots,(J-1),
\quad k=1,\dots,K\\

#### Ground truth and data

``` r

set.seed(42123)
N <- 100
J <- 3   # Number of categories
K <- 2   # Number of individual predictors (including intercept)
C <- 3   # Number of choice-based attributes

#### True parameter values
true.zeta.mu <- matrix(c(0.5, -0.3, 0.2, 0.4), J-1, K)
true.zeta.sigma <- matrix(c(0.3, 0.2, 0.25, 0.15), J-1, K)
true.gamma <- c(0.3, -0.2, 0.1)

#### Generate design matrices
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
Z <- matrix(rnorm(N * C), N, C)
for (j in 2:K) X[, j] <- CenterScale(X[, j])
for (j in 1:C) Z[, j] <- CenterScale(Z[, j])

#### Generate individual-level random coefficients
true.beta <- array(0, dim=c(J-1, K, N))
for (i in 1:N)
    true.beta[,,i] <- matrix(rnorm((J-1)*K, true.zeta.mu, true.zeta.sigma),
        J-1, K)

#### Generate response from the true model
mu.true <- matrix(tcrossprod(Z, t(true.gamma)), N, J)
for (j in 1:(J-1)) mu.true[,j] <- mu.true[,j] + rowSums(X * t(true.beta[j,,]))
phi.true <- exp(mu.true)
p.true <- phi.true / rowSums(phi.true)
y <- rcat(N, p.true)

#### Assemble Data list
S <- diag(J-1)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=array(0, dim=c(J-1,K,N)),
    gamma=rep(0,C), zeta.mu=matrix(0,J-1,K), zeta.sigma=matrix(0,J-1,K)))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.zeta.mu <- grep("zeta.mu", parm.names)
pos.zeta.sigma <- grep("zeta.sigma", parm.names)
PGF <- function(Data) {
    zeta.mu <- matrix(rnorm((Data$J-1)*Data$K), Data$J-1, Data$K)
    zeta.sigma <- matrix(runif((Data$J-1)*Data$K), Data$J-1, Data$K)
    beta <- array(rnorm((Data$J-1)*Data$K*Data$N),
        dim=c(Data$J-1, Data$K, Data$N))
    gamma <- rnorm(Data$C)
    return(c(beta, gamma, as.vector(zeta.mu), as.vector(zeta.sigma)))
    }
Data <- list(C=C, J=J, K=K, N=N, PGF=PGF, S=S, X=X, Z=Z,
    mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta,
    pos.gamma=pos.gamma, pos.zeta.mu=pos.zeta.mu,
    pos.zeta.sigma=pos.zeta.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- array(parm[Data$pos.beta], dim=c(Data$J-1, Data$K, Data$N))
    gamma <- parm[Data$pos.gamma]
    zeta.mu <- matrix(parm[Data$pos.zeta.mu], Data$J-1, Data$K)
    zeta.sigma <- matrix(interval(parm[Data$pos.zeta.sigma], 1e-100, Inf),
        Data$J-1, Data$K)
    parm[Data$pos.zeta.sigma] <- as.vector(zeta.sigma)
    ### Log-Hyperprior
    zeta.mu.prior <- sum(dnormv(zeta.mu, 0, 1000, log=TRUE))
    zeta.sigma.prior <- sum(dhalfcauchy(zeta.sigma, 25, log=TRUE))
    ### Log-Prior
    beta.prior <- sum(dnorm(beta, zeta.mu, zeta.sigma, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(tcrossprod(Data$Z, t(gamma)), Data$N, Data$J)
    for (j in 1:(Data$J-1)) mu[,j] <- rowSums(Data$X * t(beta[j, , ]))
    mu <- interval(mu, -700, 700, reflect=FALSE)
    phi <- exp(mu)
    p <- phi / rowSums(phi)
    LL <- sum(dcat(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + zeta.mu.prior + zeta.sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(p), p),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, (J-1)*K*N), rep(0, C), rep(0, (J-1)*K),
    rep(1, (J-1)*K))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: hyperparameters
post.zeta.mu <- matrix(fit$Summary2[pos.zeta.mu, "Mean"], J-1, K)
post.zeta.sigma <- matrix(fit$Summary2[pos.zeta.sigma, "Mean"], J-1, K)
cat("True zeta.mu:\n"); print(round(true.zeta.mu, 3))
cat("Posterior mean zeta.mu:\n"); print(round(post.zeta.mu, 3))
cat("True zeta.sigma:\n"); print(round(true.zeta.sigma, 3))
cat("Posterior mean zeta.sigma:\n"); print(round(post.zeta.sigma, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Discrete Choice, Multinomial Probit

The multinomial probit for discrete choice uses latent utility
differences \\\textbf{W}\\ with a multivariate normal structure,
allowing unrestricted substitution patterns through the covariance
matrix \\\Sigma\\. Individual attributes enter via category-varying
coefficients \\\beta\\, and choice-specific attributes enter via generic
coefficients \\\gamma\\. The latent utilities are constrained to be
positive when the corresponding alternative is chosen and negative
otherwise. The multinomial probit model with full covariance was
developed for Bayesian inference by McCulloch and Rossi
[\[81\]](#ref81). MNP models are applied in marketing to analyse brand
choice where the independence-of-irrelevant-alternatives assumption of
logit is violated.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires latent utility variables \\\textbf{W}\\ that are
constrained to lie in different intervals depending on the observed
choice, with a multivariate normal likelihood parameterized through a
Cholesky factor \\\textbf{U}\\. The conditional truncation of latent
variables and the Cholesky decomposition of the covariance matrix fall
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Form

\\\textbf{W}\_{i,1:(J-1)} \sim \mathcal{N}\_{J-1}(\mu\_{i,1:(J-1)},
\Sigma), \quad i=1,\dots,N\\ \\\textbf{W}\_{i,j} \in \left\\
\begin{array}{l l} \$\[0,10\]\$ & \quad \mbox{if \$\textbf{y}\_i =
j\$}\\ \$\[-10,0\]\$ \\ \end{array} \right. \\ \\\mu\_{1:N,j} =
\textbf{X} \beta\_{j,1:K} + \textbf{Z} \gamma\\ \\\Sigma = \textbf{U}^T
\textbf{U}\\ \\\beta\_{j,k} \sim \mathcal{N}(0, 10), \quad
j=1,\dots,(J-1), \quad k=1,\dots,K\\ \\\gamma_c \sim \mathcal{N}(0, 10),
\quad c=1,\dots,C\\ \\\textbf{U}\_{j,k} \sim \mathcal{N}(0,1), \quad
j=1,\dots,(J-1), \quad k=1,\dots,(J-1), \quad j \ge k, \quad j \ne k =
1\\

#### Ground truth and data

``` r

set.seed(42124)
N <- 150
J <- 3   # Number of categories
K <- 3   # Number of individual predictors (including intercept)
C <- 4   # Number of choice-based attributes

#### True parameter values
true.beta <- matrix(c(0.4, -0.3, 0.6, -0.2, 0.3, -0.5), J-1, K)
true.gamma <- c(0.2, -0.3, 0.1, -0.2)

#### Generate design matrices
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
Z <- matrix(rnorm(N * C), N, C)
for (j in 2:K) X[, j] <- CenterScale(X[, j])
for (j in 1:C) Z[, j] <- CenterScale(Z[, j])

#### Generate response from the true model (via MNP latent utilities)
mu.true <- tcrossprod(X, true.beta) +
    as.vector(tcrossprod(Z, t(true.gamma)))
W.true <- mu.true + rmvn(N, rep(0, J-1), diag(J-1))
y <- max.col(cbind(W.true, 0))

#### Assemble Data list
S <- diag(J-1)
U <- matrix(NA, J-1, J-1)
U[upper.tri(U, diag=TRUE)] <- 0
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,(J-1),K), gamma=rep(0,C),
    U=U, W=matrix(0,N,J-1)))
parm.names <- parm.names[-which(parm.names == "U[1,1]")]
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.U <- grep("U", parm.names)
pos.W <- grep("W", parm.names)
PGF <- function(Data) {
    beta <- rnorm((Data$J-1)*Data$K)
    gamma <- rnorm(Data$C)
    U <- rnorm((Data$J-2) + (factorial(Data$J-1) /
        (factorial(Data$J-1-2)*factorial(2))), 0, 1)
    W <- matrix(runif(Data$N*(Data$J-1), -10, 0), Data$N, Data$J-1)
    Y <- as.indicator.matrix(Data$y)
    W <- ifelse(Y[,-Data$J] == 1, abs(W), W)
    return(c(beta, gamma, U, as.vector(W)))}
Data <- list(C=C, J=J, K=K, N=N, PGF=PGF, S=S, X=X, Z=Z,
    mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta,
    pos.gamma=pos.gamma, pos.U=pos.U, pos.W=pos.W, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$J-1, Data$K)
    gamma <- parm[Data$pos.gamma]
    u <- c(0, parm[Data$pos.U])
    U <- diag(Data$J-1)
    U[upper.tri(U, diag=TRUE)] <- u
    diag(U) <- exp(diag(U))
    Sigma <- t(U) %*% U
    Sigma[1,] <- Sigma[,1] <- U[1,]
    W <- matrix(parm[Data$pos.W], Data$N, Data$J-1)
    Y <- as.indicator.matrix(Data$y)
    temp <- which(Y[,-c(Data$J)] == 1)
    W[temp] <- interval(W[temp], 0, 10)
    temp <- which(Y[,-c(Data$J)] == 0)
    W[temp] <- interval(W[temp], -10, 0)
    parm[Data$pos.W] <- as.vector(W)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 10, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 10, log=TRUE))
    U.prior <- sum(dnorm(u[-1], 0, 1, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, beta) +
        as.vector(tcrossprod(Data$Z, t(gamma)))
    #eta <- exp(cbind(mu,0))
    #p <- eta / rowSums(eta)
    LL <- sum(dmvn(W, mu, Sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=max.col(cbind(rmvn(nrow(mu), mu, Sigma),0)), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- GIV(Model, Data, PGF=TRUE)
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], J-1, K)
post.gamma <- fit$Summary2[pos.gamma, "Mean"]
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))
cat("True gamma:     ", round(true.gamma, 3), "\n")
cat("Posterior mean: ", round(post.gamma, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Mixed Logit

The mixed logit (without choice-specific attributes) models categorical
choice with individual-level random coefficients \\\beta\_{j,k,i}\\
drawn from normal hyperpriors. Unlike the DC mixed logit variant, this
formulation includes only individual attributes \\\textbf{X}\\ with
random slopes, making it appropriate when the heterogeneity of interest
lies in how individuals weight the same set of predictors differently.
The mixed logit model allows random coefficients in discrete choice, as
developed by McFadden and Train [\[84\]](#ref84) and synthesised by
Train [\[109\]](#ref109). Mixed logit is applied in environmental
economics to estimate heterogeneous willingness-to-pay for conservation
programmes through choice experiments.

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a three-dimensional array of random coefficients
\\\beta\_{j,k,i}\\ indexed by category, predictor, and observation, with
hierarchical normal hyperpriors. This high-dimensional
random-coefficient structure falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Form

\\\textbf{y}\_i \sim \mathcal{CAT}(\textbf{p}\_{i,1:J}), \quad
i=1,\dots,N\\ \\\textbf{p}\_{i,j} = \frac{\phi\_{i,j}}{\sum^J\_{j=1}
\phi\_{i,j}}\\ \\\phi = \exp(\mu)\\ \\\mu\_{i,j} = \beta\_{j,1:K,i}
\textbf{X}\_{i,1:K} \in \[-700,700\], \quad i=1,\dots,N, \quad
j=1,\dots,(J-1)\\ \\\mu\_{i,J} = 0\\ \\\beta\_{j,k,i} \sim
\mathcal{N}(\zeta^\mu\_{j,k}, \zeta^\sigma2\_{j,k}), \quad i=1,\dots,N,
\quad j=1,\dots,(J-1), \quad k=1,\dots,K\\ \\\zeta^\mu\_{j,k} \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,(J-1), \quad k=1,\dots,K\\
\\\zeta^\sigma\_{j,k} \sim \mathcal{HC}(25), \quad j=1,\dots,(J-1),
\quad k=1,\dots,K\\

#### Ground truth and data

``` r

set.seed(42175)
N <- 100
J <- 3   # Number of categories
K <- 2   # Number of predictors (including intercept)

#### True parameter values
true.zeta.mu <- matrix(c(0.6, -0.4, 0.3, 0.5), J-1, K)
true.zeta.sigma <- matrix(c(0.25, 0.20, 0.30, 0.15), J-1, K)

#### Generate design matrix
X <- cbind(1, matrix(rnorm(N * (K - 1)), N, K - 1))
for (j in 2:K) X[, j] <- CenterScale(X[, j])

#### Generate individual-level random coefficients
true.beta <- array(0, dim=c(J-1, K, N))
for (i in 1:N)
    true.beta[,,i] <- matrix(rnorm((J-1)*K, true.zeta.mu, true.zeta.sigma),
        J-1, K)

#### Generate response from the true model
mu.true <- matrix(0, N, J)
for (j in 1:(J-1)) mu.true[,j] <- rowSums(X * t(true.beta[j,,]))
phi.true <- exp(mu.true)
p.true <- phi.true / rowSums(phi.true)
y <- rcat(N, p.true)

#### Assemble Data list
S <- diag(J-1)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=array(0, dim=c(J-1,K,N)),
    zeta.mu=matrix(0,J-1,K), zeta.sigma=matrix(0,J-1,K)))
pos.beta <- grep("beta", parm.names)
pos.zeta.mu <- grep("zeta.mu", parm.names)
pos.zeta.sigma <- grep("zeta.sigma", parm.names)
PGF <- function(Data) {
    zeta.mu <- matrix(rnorm((Data$J-1)*Data$K), Data$J-1, Data$K)
    zeta.sigma <- matrix(runif((Data$J-1)*Data$K), Data$J-1, Data$K)
    beta <- array(rnorm((Data$J-1)*Data$K*Data$N),
        dim=c(Data$J-1, Data$K, Data$N))
    return(c(beta, as.vector(zeta.mu), as.vector(zeta.sigma)))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, S=S, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.zeta.mu=pos.zeta.mu,
    pos.zeta.sigma=pos.zeta.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- array(parm[Data$pos.beta], dim=c(Data$J-1, Data$K, Data$N))
    zeta.mu <- matrix(parm[Data$pos.zeta.mu], Data$J-1, Data$K)
    zeta.sigma <- matrix(interval(parm[Data$pos.zeta.sigma], 1e-100, Inf),
        Data$J-1, Data$K)
    parm[Data$pos.zeta.sigma] <- as.vector(zeta.sigma)
    ### Log-Hyperprior
    zeta.mu.prior <- sum(dnormv(zeta.mu, 0, 1000, log=TRUE))
    zeta.sigma.prior <- sum(dhalfcauchy(zeta.sigma, 25, log=TRUE))
    ### Log-Prior
    beta.prior <- sum(dnorm(beta, zeta.mu, zeta.sigma, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(0, Data$N, Data$J)
    for (j in 1:(Data$J-1)) mu[,j] <- rowSums(Data$X * t(beta[j, , ]))
    mu <- interval(mu, -700, 700, reflect=FALSE)
    phi <- exp(mu)
    p <- phi / rowSums(phi)
    LL <- sum(dcat(Data$y, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + zeta.mu.prior + zeta.sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rcat(nrow(p), p),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, (J-1)*K*N), rep(0, (J-1)*K), rep(1, (J-1)*K))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: hyperparameters
post.zeta.mu <- matrix(fit$Summary2[pos.zeta.mu, "Mean"], J-1, K)
post.zeta.sigma <- matrix(fit$Summary2[pos.zeta.sigma, "Mean"], J-1, K)
cat("True zeta.mu:\n"); print(round(true.zeta.mu, 3))
cat("Posterior mean zeta.mu:\n"); print(round(post.zeta.mu, 3))
cat("True zeta.sigma:\n"); print(round(true.zeta.sigma, 3))
cat("Posterior mean zeta.sigma:\n"); print(round(post.zeta.sigma, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multiple Discrete-Continuous Choice

This form of a multivariate discrete-continuous choice model was
introduced in Kim et al. [\[62\]](#ref62) and referred to as a variety
model. The original version is presented with log-normally distributed
errors, but a gamma regression form is used here instead, which has
always mixed better in testing. The \\\gamma\\ parameters are fixed at 1
for identifiability, as recommended by the original authors. The model
accommodates simultaneous purchase decisions across \\J\\ products by
\\G\\ decision-makers, with multilevel random effects on the baseline
preferences \\\psi_1\\. The multiple discrete-continuous choice model
was introduced by Bhat [\[11\]](#ref11) for situations where consumers
choose both which alternatives to consume and how much. MDCC models are
applied in time-use studies to analyse how individuals allocate their
leisure time across multiple activities simultaneously.

#### Form

\\\textbf{Y} \sim \mathcal{G}(\lambda\tau, \tau)\\ \\\lambda\_{i,j} =
\exp(\textbf{Z}\_{i,j}\log(\psi1\_{m\[i\],j}) +
\textbf{X1}\_{i,1:K}\log(\beta) +
\textbf{X2}\_{i,1:L}\log(\delta))(\textbf{Y}\_{i,j} +
\gamma_j)^\alpha_j), \quad i=1,\dots,N, j=1,\dots,J\\ \\\alpha_j \sim
\mathcal{U}(0,1), \quad j=1,\dots,J\\ \\\log(\beta_k) \sim
\mathcal{N}(0,1000), \quad k=1,\dots,K\\ \\\gamma_j = 1, \quad
j=1,\dots,J\\ \\\log(\delta\_{j,l}) \sim \mathcal{N}(0,1000), \quad
j=1,\dots,(J-1), \quad l=1,\dots,L\\ \\\log(\psi0_j) \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\ \\\log(\psi1\_{g,j}) \sim
\mathcal{N}\_{J}(\log(\psi0), \Omega^{-1}), \quad g=1,\dots,G, \quad
j1=,\dots,J\\ \\\Omega \sim \mathcal{W}\_{J+1}(\textbf{S}), \quad
\textbf{S} = \textbf{I}\_J\\ \\\tau_j \sim \mathcal{HC}(25), \quad
j=1,\dots,J\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it combines a gamma likelihood with log-space multilevel random
effects, Wishart-distributed precision, satiation parameters
\\\alpha_j\\, and a complex product-attribute interaction structure, all
of which fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42181)
G <- 6   # decision-makers
J <- 3   # products
K <- 4   # product attributes
L <- 5   # decision-maker attributes
N <- 30  # records

#### Simulate MDCC data
X1 <- matrix(rnorm(N*K), N, K)
X2 <- matrix(rnorm(N*L), N, L)
Sigma <- matrix(runif((J-1)*(J-1),-1,1),J-1,J-1)
diag(Sigma) <- runif(J-1,1,5)
Sigma <- as.positive.definite(Sigma) / 100
alpha <- runif(J)
log.beta <- rnorm(K,0,0.1)
log.delta <- matrix(rnorm((J-1)*L,0,0.1), J-1, L)
log.psi0 <- rnorm(J)
log.psi1 <- rmvn(G, log.psi0, Sigma)
m <- rcat(N, rep(1/G,G))
Z <- as.indicator.matrix(m)
Y <- matrix(0, N, J)
Y <- round(exp(tcrossprod(Z, t(cbind(log.psi1,0))) +
    matrix(tcrossprod(X1, t(log.beta)), N, J) +
    tcrossprod(X2, rbind(log.delta, colSums(log.delta)*-1))) *
    (Y + 1)^ matrix(alpha,N,J,byrow=TRUE) +
    matrix(rnorm(N*J,0,0.1),N,J))

#### Assemble Data list
S <- diag(J)
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J), log.beta=rep(0,K),
    log.delta=matrix(0,J-1,L), log.psi0=rep(0,J),
    log.psi1=matrix(0,G,J), tau=rep(0,J), U=S),
    uppertri=c(0,0,0,0,0,0,1))
pos.alpha <- grep("alpha", parm.names)
pos.log.beta <- grep("log.beta", parm.names)
pos.log.delta <- grep("delta", parm.names)
pos.log.psi0 <- grep("log.psi0", parm.names)
pos.log.psi1 <- grep("log.psi1", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    alpha <- runif(Data$J,0.9,1)
    log.beta <- rnorm(Data$K,0,0.1)
    log.delta <- rnorm((Data$J-1)*Data$L,0,0.1)
    log.psi0 <- rnorm(Data$J)
    U <- rwishartc(Data$J+1, Data$S)
    log.psi1 <- as.vector(rmvnpc(Data$G, log.psi0, U))
    tau <- runif(Data$J)
    return(c(alpha, log.beta, log.delta, log.psi0, log.psi1, tau,
        U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(G=G, J=J, K=K, L=L, N=N, PGF=PGF, S=S, X1=X1, X2=X2, Y=Y,
    Z=Z, m=m, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.log.beta=pos.log.beta,
    pos.log.delta=pos.log.delta, pos.log.psi0=pos.log.psi0,
    pos.log.psi1=pos.log.psi1, pos.tau=pos.tau)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[Data$pos.alpha] <- alpha <- interval(parm[Data$pos.alpha], 0, 1)
    log.beta <- parm[Data$pos.log.beta]
    log.delta <- matrix(parm[Data$pos.log.delta], Data$J-1, Data$L)
    log.psi0 <- parm[Data$pos.log.psi0]
    log.psi1 <- matrix(parm[Data$pos.log.psi1], Data$G, Data$J)
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    U <- as.parm.matrix(U, Data$J, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    lambda <- tcrossprod(Data$Z, t(log.psi1)) +
        matrix(tcrossprod(Data$X1, t(log.beta)), Data$N, Data$J) +
        tcrossprod(Data$X2, rbind(log.delta, colSums(log.delta)*-1))
    ### Log-Prior
    U.prior <- dwishartc(U, Data$J+1, Data$S, log=TRUE)
    alpha.prior <- sum(dunif(alpha, 0, 1, log=TRUE))
    log.beta.prior <- sum(dnormv(log.beta, 0, 1000, log=TRUE))
    log.delta.prior <- sum(dnormv(log.delta, 0, 1000, log=TRUE))
    log.psi0.prior <- sum(dnormv(log.psi0, 0, 1000, log=TRUE))
    log.psi1.prior <- sum(dmvnpc(lambda,
        matrix(log.psi0, Data$N, Data$J, byrow=TRUE), U, log=TRUE))
    tau.prior <- sum(dhalfcauchy(tau, 25, log=TRUE))
    ### Log-Likelihood
    alpha <- matrix(alpha, Data$N, Data$J, byrow=TRUE)
    lambda <- exp(lambda)*(Data$Y + 1)^ alpha
    tau <- matrix(tau, Data$N, Data$J, byrow=TRUE)
    LL <- sum(dgamma(Data$Y+1, lambda*tau, tau, log=TRUE))
    ### Log-Posterior
    LP <- LL + U.prior + alpha.prior + log.beta.prior + log.delta.prior +
        log.psi0.prior + log.psi1.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rgamma(prod(dim(lambda)), lambda*tau, tau)-1,
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(runif(J,0.9,1), rnorm(K,0,0.1),
    rnorm((J-1)*L,0,0.1), rnorm(J,0,0.1),
    rmvnpc(G, rnorm(J,0,0.1), rwishartc(J+1,S)), runif(J),
    upper.triangle(rwishartc(J+1,S), diag=TRUE))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Multivariate and compositional models

### Additive Main Effects and Multiplicative Interaction (AMMI) **\[NEW\]**

The AMMI model decomposes genotype-by-environment interaction tables
into additive main effects (genotype and environment means) plus a
multiplicative interaction term obtained by singular value decomposition
of the residual matrix. The multiplicative component captures structured
interaction patterns that purely additive models miss, making AMMI a
cornerstone of multi-environment trial analysis. AMMI was systematised
by Gauch [\[159\]](#ref159), who showed that retaining the first few
multiplicative components often explains most of the
genotype-environment interaction variance while filtering noise. AMMI
models are standard in plant breeding programmes for identifying stable
crop varieties across diverse agro-ecological zones.

#### Form

\\Y\_{ij} \sim \mathcal{N}(\mu\_{ij}, \sigma^2)\\ \\\mu\_{ij} = \mu +
\alpha_i + \beta_j + \sum\_{q=1}^{Q} \lambda_q \gamma\_{iq}
\delta\_{jq}\\ \\\mu \sim \mathcal{N}(0, 1000)\\ \\\alpha_i \sim
\mathcal{N}(0, 1000)\\ \\\beta_j \sim \mathcal{N}(0, 1000)\\ \\\lambda
\sim \mathcal{HC}(25), \quad \lambda \> 0\\ \\\gamma\_{iq} \sim
\mathcal{N}(0, 1), \quad \delta\_{jq} \sim \mathcal{N}(0, 1)\\ \\\sigma
\sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean function contains a bilinear interaction term
\\\lambda_q \gamma\_{iq} \delta\_{jq}\\ that is a product of three
estimated parameters. The DSL supports linear combinations and standard
link functions but not arbitrary multiplicative interactions among
latent parameters. The manual `Model` function below handles this
directly.

#### Ground truth and data

``` r

set.seed(42201)
I <- 5   # genotypes
J <- 9   # environments
Q <- 1   # multiplicative components
N <- I * J

#### True parameter values
true.mu <- 100
true.alpha <- c(-1, -1, 0, 1, 1)
true.beta <- -4:4
true.lambda <- 12
true.gamma <- c(0.63, -0.32, 0, 0.32, -0.63)
true.delta <- c(0.5, 0.5, 0, 0, 0, 0, 0, -0.5, -0.5)
true.sigma <- 1.5

#### Generate data on genotype-environment grid
grid <- expand.grid(geno = 1:I, env = 1:J)
mu.vec <- true.mu + true.alpha[grid$geno] + true.beta[grid$env] +
    true.lambda * true.gamma[grid$geno] * true.delta[grid$env]
y <- rnorm(N, mu.vec, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, alpha=rep(0,I), beta=rep(0,J),
    lambda=0, gamma=rep(0,I), delta=rep(0,J), sigma=0))
pos.mu <- grep("^mu$", parm.names)
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    mu <- rnorm(1, 100, 10)
    alpha <- rnorm(Data$I)
    beta <- rnorm(Data$J)
    lambda <- abs(rnorm(1, 0, 5))
    gamma <- rnorm(Data$I)
    delta <- rnorm(Data$J)
    sigma <- rhalfcauchy(1, 5)
    return(c(mu, alpha, beta, lambda, gamma, delta, sigma))
    }
Data <- list(I=I, J=J, Q=Q, N=N, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.alpha=pos.alpha,
    pos.beta=pos.beta, pos.lambda=pos.lambda, pos.gamma=pos.gamma,
    pos.delta=pos.delta, pos.sigma=pos.sigma, geno=grid$geno,
    env=grid$env, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    gamma <- parm[Data$pos.gamma]
    gamma[1] <- interval(gamma[1], 1e-100, Inf)
    parm[Data$pos.gamma[1]] <- gamma[1]
    delta <- parm[Data$pos.delta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 1000, log=TRUE)
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    lambda.prior <- dhalfcauchy(lambda, 25, log=TRUE)
    gamma.prior <- sum(dnorm(gamma, 0, 1, log=TRUE))
    delta.prior <- sum(dnorm(delta, 0, 1, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu.vec <- mu + alpha[Data$geno] + beta[Data$env] +
        lambda * gamma[Data$geno] * delta[Data$env]
    LL <- sum(dnorm(Data$y, mu.vec, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + alpha.prior + beta.prior + lambda.prior +
        gamma.prior + delta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu.vec, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(100, rep(0, I), rep(0, J), 5, rep(0.1, I),
    rep(0.1, J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
cat("mu -- true:", true.mu,
    " post. mean:", round(fit$Summary2[pos.mu, "Mean"], 3), "\n")
cat("lambda -- true:", true.lambda,
    " post. mean:", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha]  <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta]  <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma]  <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.mu]     <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Compositional Covariates Regression **\[NEW\]**

Compositional data lie on the simplex (parts sum to a constant,
typically one) and cannot enter standard regression directly because of
the resulting perfect multicollinearity. The isometric log-ratio (ILR)
transform maps the \\D\\-part simplex to \\(D-1)\\-dimensional Euclidean
space using an orthonormal basis, after which standard linear regression
applies without constraint. The ILR coordinates are interpretable as
sequential log-contrasts between groups of parts defined by a chosen
partition. The compositional data framework was established by Aitchison
[\[160\]](#ref160), and Egozcue et al. [\[161\]](#ref161) introduced the
ILR transform as the theoretically optimal coordinate system for the
simplex. Compositional regression with ILR coordinates is standard in
geochemistry for predicting rock formation age from mineral oxide
proportions.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \beta_1 + \beta_2
\\ \textbf{x}\_{\text{ilr},1} + \beta_3 \\ \textbf{x}\_{\text{ilr},2}\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,3\\ \\\sigma \sim
\mathcal{HC}(25)\\

where \\\textbf{x}\_{\text{ilr}}\\ is the ILR transform of a \\D=3\\
compositional covariate.

#### model_spec() notation

Once the ILR coordinates are precomputed, this model is standard linear
regression and is expressible in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = beta1 + beta2 * x_ilr1 + beta3 * x_ilr2
  beta1 ~ Normal(0, sqrt(1000))
  beta2 ~ Normal(0, sqrt(1000))
  beta3 ~ Normal(0, sqrt(1000))
  sigma ~ HalfCauchy(25)
")
```

#### Ground truth and data

``` r

set.seed(42207)
N <- 200
D <- 3   # compositional parts

#### True parameter values
true.beta <- c(2, 1.5, -0.8)
true.sigma <- 0.5

#### Generate compositional covariates from Dirichlet(5,5,5)
x.comp <- matrix(rgamma(N * D, shape=5), nrow=N, ncol=D)
x.comp <- x.comp / rowSums(x.comp)

#### ILR transform (Helmert sub-composition for D=3)
x.ilr1 <- sqrt(1/2) * log(x.comp[,1] / x.comp[,2])
x.ilr2 <- sqrt(2/3) * log(sqrt(x.comp[,1] * x.comp[,2]) / x.comp[,3])

#### Design matrix with intercept
X <- cbind(1, x.ilr1, x.ilr2)

#### Generate response
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
J <- ncol(X)
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- rhalfcauchy(1, 5)
    return(c(beta, sigma))
    }
Data <- list(N=N, J=J, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma,
    X=X, x_ilr1=x.ilr1, x_ilr2=x.ilr2, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- Data$X %*% beta
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
for (j in 1:J)
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(fit$Summary2[pos.beta[j], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Dirichlet Regression **\[NEW\]**

Dirichlet regression models a compositional response \\\textbf{y}\_i\\
on the simplex by assigning each of \\K\\ components its own linear
predictor through a log link, so that covariates can shift the expected
proportion of each component in different directions while the Dirichlet
distribution captures the inherent negative correlation among parts.
Each component has a separate intercept and set of slope coefficients,
yielding a flexible model for multivariate proportions bounded between
zero and one. The Dirichlet regression framework was developed by Hijazi
and Jernigan [\[139\]](#ref139) and implemented in software by Maier
[\[140\]](#ref140). Dirichlet regression is used in political science to
model vote shares across parties as a function of economic indicators
and demographic covariates.

#### Form

\\\textbf{y}\_i \sim \mathcal{D}(\alpha_i), \quad i=1,\dots,N\\
\\\log(\alpha\_{i,k}) = \textbf{X}\_i \beta_k, \quad k=1,\dots,K\\
\\\beta\_{k,j} \sim \mathcal{N}(0, 100), \quad j=1,\dots,J, \quad
k=1,\dots,K\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the response is a vector on the simplex with
component-specific linear predictors, each producing a separate
Dirichlet concentration parameter, requiring a multivariate
distributional structure that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42210)
N <- 200
K <- 3    # Compositional parts
J <- 3    # Predictors (intercept + 2 covariates)

#### True parameter values (rows = predictors, cols = components)
true.beta <- matrix(c(2.0,  1.0,  0.5,
                       0.5, -0.3,  0.2,
                      -0.8,  0.4,  0.1), nrow=J, byrow=TRUE)

#### Generate design matrix
x1 <- rnorm(N)
x2 <- rnorm(N)
X <- cbind(1, x1, x2)

#### Generate response
log.alpha.true <- X %*% true.beta
alpha.true <- exp(log.alpha.true)
y <- matrix(NA, N, K)
for (i in 1:N) y[i,] <- rdirichlet(1, alpha.true[i,])

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0, J*K)))
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J * Data$K, 0, 0.5)
    return(beta)
    }
Data <- list(N=N, J=J, K=K, X=X, y=y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$J, Data$K)
    ### Log-Prior
    beta.prior <- sum(dnormv(as.vector(beta), 0, 100, log=TRUE))
    ### Log-Likelihood
    log.alpha <- Data$X %*% beta
    alpha <- exp(log.alpha)
    LL <- 0
    for (i in 1:Data$N) LL <- LL + ddirichlet(Data$y[i,], alpha[i,], log=TRUE)
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rdirichlet(Data$N, alpha[Data$N,])[,1], parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- rep(0, J*K)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], J, K)
cat("True beta matrix:\n"); print(round(true.beta, 3))
cat("Posterior mean beta matrix:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### MANCOVA

The multivariate analysis of covariance extends the univariate ANCOVA
(see section [link](#ancova)) to multiple response variables
simultaneously, modelling their joint distribution through a
multivariate normal likelihood with a Wishart prior on the precision
matrix. Each response has its own intercept, factor effects with
sum-to-zero constraints, and covariate regression coefficients, while
the cross-response dependence is captured by the shared precision matrix
\\\Omega\\. Multivariate analysis of covariance extends Wilks’ (1932)
MANOVA framework by including continuous covariates alongside
categorical factors. MANCOVA is applied in neuroimaging to compare
multiple brain region volumes across clinical groups while adjusting for
age and total intracranial volume.

#### Form

\\\textbf{Y}\_{i,1:J} \sim \mathcal{N}\_K(\mu\_{i,1:J}, \Sigma), \quad
i=1,\dots,N\\ \\\mu\_{i,k} = \alpha_k + \beta\_{k,\textbf{X}\[i,1\]} +
\gamma\_{k,\textbf{X}\[i,1\]} + \textbf{X}\_{1:N,3:(C+J)}
\delta\_{k,1:C}\\ \\\epsilon\_{i,k} = \textbf{Y}\_{i,k} - \mu\_{i,k}\\
\\\alpha_k \sim \mathcal{N}(0, 1000), \quad k=1,\dots,K\\ \\\beta\_{k,l}
\sim \mathcal{N}(0, \sigma^2_1), \quad l=1,\dots,(L-1)\\
\\\beta\_{1:K,L} = - \sum^{L-1}\_{l=1} \beta\_{1:K,l}\\ \\\gamma\_{k,m}
\sim \mathcal{N}(0, \sigma^2_2), \quad m=1,\dots,(M-1)\\
\\\gamma\_{1:K,M} = - \sum^{M-1}\_{m=1} \gamma\_{1:K,m}\\
\\\delta\_{k,c} \sim \mathcal{N}(0, 1000)\\ \\\Omega \sim
\mathcal{W}\_{K+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_K\\
\\\Sigma = \Omega^{-1}\\ \\\sigma\_{1:J} \sim \mathcal{HC}(25)\\

#### model_spec() notation

The MANCOVA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a Wishart prior on the Cholesky factor of the
precision matrix, sum-to-zero constraints on the factor effects, and
hierarchical variance components on the group-level standard deviations.
These constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42193)
C <- 2   # Number of covariates
J <- 2   # Number of factors (treatments)
K <- 3   # Number of endogenous (dependent) variables
L <- 4   # Number of levels in factor (treatment) 1
M <- 5   # Number of levels in factor (treatment) 2
N <- 100

#### True parameter values
true.alpha <- c(0.8, -0.5, 1.2)
true.beta <- matrix(c(-1.0, 0.5, 0.3, 0.2,
                       0.7, -0.3, 0.1, -0.5,
                      -0.4, 0.6, -0.8, 0.6), K, L, byrow = TRUE)
true.beta[,L] <- -rowSums(true.beta[,-L])
true.gamma <- matrix(c(0.4, -0.2, 0.5, -0.3, -0.4,
                       -0.6, 0.3, 0.1, 0.4, -0.2,
                        0.2, -0.5, 0.3, 0.2, -0.2), K, M, byrow = TRUE)
true.gamma[,M] <- -rowSums(true.gamma[,-M])
true.delta <- matrix(c(0.9, -0.4,
                       0.3,  0.7,
                      -0.5,  0.2), K, C, byrow = TRUE)

#### Generate design matrix and response
X <- matrix(c(rcat(N, rep(1/L, L)), rcat(N, rep(1/M, M)),
    runif(N * C, 0, 1)), N, J + C)
Y <- matrix(NA, N, K)
for (k in 1:K) {
    Y[,k] <- true.alpha[k] + true.beta[k, X[,1]] +
        true.gamma[k, X[,2]] +
        tcrossprod(true.delta[k,], X[,-c(1,2)]) + rnorm(N, 0, 0.3)
}
S <- diag(K)

#### Assemble Data list
mon.names <- c("LP", "s.o.beta", "s.o.gamma", "s.o.epsilon",
    as.parm.names(list(s.beta=rep(0,K), s.gamma=rep(0,K),
    s.epsilon=rep(0,K))))
parm.names <- as.parm.names(list(alpha=rep(0,K), beta=matrix(0,K,(L-1)),
    gamma=matrix(0,K,(M-1)), delta=matrix(0,K,C), U=diag(K),
    sigma=rep(0,2)), uppertri=c(0,0,0,0,1,0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$K)
    sigma <- runif(2)
    beta <- rnorm(Data$K*(Data$L-1), 0, sigma[1])
    gamma <- rnorm(Data$K*(Data$M-1), 0, sigma[2])
    delta <- rnorm(Data$K*Data$C)
    U <- rwishartc(Data$K+1, Data$S)
    return(c(alpha, beta, gamma, delta, U[upper.tri(U, diag=TRUE)],
        sigma))
    }
Data <- list(C=C, J=J, K=K, L=L, M=M, N=N, PGF=PGF, S=S, X=X, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.beta=pos.beta, pos.gamma=pos.gamma, pos.delta=pos.delta,
    pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- matrix(c(parm[Data$pos.beta], rep(0,Data$K)), Data$K, Data$L)
    beta[,Data$L] <- -rowSums(beta[,-Data$L])
    gamma <- matrix(c(parm[Data$pos.gamma],
        rep(0,Data$K)), Data$K, Data$M)
    gamma[,Data$M] <- -rowSums(gamma[,-Data$M])
    delta <- matrix(parm[Data$pos.delta], Data$K, Data$C)
    U <- as.parm.matrix(U, Data$K, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnorm(beta, 0, sigma[1], log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma[2], log=TRUE))
    delta.prior <- sum(dnormv(delta, 0, 1000, log=TRUE))
    U.prior <- dwishartc(U, Data$K+1, Data$S, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(0,Data$N,Data$K)
    for (k in 1:Data$K) {
        mu[,k] <- alpha[k] + beta[k,Data$X[,1]] + gamma[k,Data$X[,2]] +
        tcrossprod(Data$X[,-c(1,2)], t(delta[k,]))}
    LL <- sum(dmvnpc(Data$Y, mu, U, log=TRUE))
    ### Variance Components, Omnibus
    s.o.beta <- sd(as.vector(beta))
    s.o.gamma <- sd(as.vector(gamma))
    s.o.epsilon <- sd(as.vector(Data$Y - mu))
    ### Variance Components, Univariate
    s.beta <- sqrt(apply(beta, 1, var))
    s.gamma <- sqrt(apply(gamma, 1, var))
    s.epsilon <- sqrt(apply(Data$Y - mu, 2, var))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + delta.prior +
        U.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, s.o.beta, s.o.gamma,
        s.o.epsilon, s.beta, s.gamma, s.epsilon),
        yhat=rmvnpc(nrow(mu), mu, U), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K), rep(0,K*(L-1)), rep(0,K*(M-1)),
    rep(0,C*K), upper.triangle(S, diag=TRUE), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: intercepts
for (k in 1:K)
    cat("alpha[", k, "] -- true:", true.alpha[k],
        " post. mean:", round(fit$Summary2[pos.alpha[k], "Mean"], 3), "\n")

#### Parameter recovery: covariate effects
post.delta <- matrix(fit$Summary2[pos.delta, "Mean"], K, C)
for (k in 1:K)
    for (cc in 1:C)
        cat("delta[", k, ",", cc, "] -- true:", true.delta[k, cc],
            " post. mean:", round(post.delta[k, cc], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### MANOVA

The multivariate analysis of variance extends the two-way ANOVA (see
section [link](#anova.two.way)) to \\K\\ response variables jointly,
using a multivariate normal likelihood with a Wishart prior on the
precision matrix. Two crossed factors with sum-to-zero constraints
partition the mean structure, while hierarchical half-Cauchy priors
govern the group-level standard deviations. Multivariate analysis of
variance was formalised by Wilks [\[114\]](#ref114) and its Bayesian
treatment developed in the context of multivariate normal theory. MANOVA
is standard in pharmaceutical research for simultaneously comparing
multiple efficacy endpoints across treatment arms.

#### Form

\\\textbf{Y}\_{i,1:J} \sim \mathcal{N}\_K(\mu\_{i,1:J}, \Omega^{-1}),
\quad i=1,\dots,N\\ \\\mu\_{i,k} = \alpha_k +
\beta\_{k,\textbf{X}\[i,1\]} + \gamma\_{k,\textbf{X}\[i,1\]}\\
\\\epsilon\_{i,k} = \textbf{Y}\_{i,k} - \mu\_{i,k}\\ \\\alpha_k \sim
\mathcal{N}(0, 1000), \quad k=1,\dots,K\\ \\\beta\_{k,l} \sim
\mathcal{N}(0, \sigma^2_1), \quad l=1,\dots,(L-1)\\ \\\beta\_{1:K,L} = -
\sum^{L-1}\_{l=1} \beta\_{1:K,l}\\ \\\gamma\_{k,m} \sim \mathcal{N}(0,
\sigma^2_2), \quad m=1,\dots,(M-1)\\ \\\gamma\_{1:K,M} = -
\sum^{M-1}\_{m=1} \gamma\_{1:K,m}\\ \\\Omega \sim
\mathcal{W}\_{K+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_K\\
\\\sigma\_{1:J} \sim \mathcal{HC}(25)\\

#### model_spec() notation

The MANOVA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a Wishart prior on the Cholesky factor of the
precision matrix and sum-to-zero constraints on the factor effects.
These constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42194)
J <- 2   # Number of factors (treatments)
K <- 3   # Number of endogenous (dependent) variables
L <- 4   # Number of levels in factor (treatment) 1
M <- 5   # Number of levels in factor (treatment) 2
N <- 100

#### True parameter values
true.alpha <- c(1.0, -0.3, 0.6)
true.beta <- matrix(c(-0.8, 0.4, 0.5, -0.1,
                       0.6, -0.5, 0.2, -0.3,
                      -0.3, 0.7, -0.6, 0.2), K, L, byrow = TRUE)
true.beta[,L] <- -rowSums(true.beta[,-L])
true.gamma <- matrix(c(0.3, -0.4, 0.6, -0.2, -0.3,
                       -0.5, 0.2, 0.4, 0.1, -0.2,
                        0.1, -0.3, 0.2, 0.3, -0.3), K, M, byrow = TRUE)
true.gamma[,M] <- -rowSums(true.gamma[,-M])

#### Generate design matrix and response
X <- cbind(rcat(N, rep(1/L, L)), rcat(N, rep(1/M, M)))
Y <- matrix(NA, N, K)
for (k in 1:K) {
    Y[,k] <- true.alpha[k] + true.beta[k, X[,1]] +
        true.gamma[k, X[,2]] + rnorm(N, 0, 0.3)
}
S <- diag(K)

#### Assemble Data list
mon.names <- c("LP", "s.o.beta", "s.o.gamma", "s.o.epsilon",
    as.parm.names(list(s.beta=rep(0,K), s.gamma=rep(0,K),
    s.epsilon=rep(0,K))))
parm.names <- as.parm.names(list(alpha=rep(0,K), beta=matrix(0,K,(L-1)),
    gamma=matrix(0,K,(M-1)), U=diag(K), sigma=rep(0,2)),
    uppertri=c(0,0,0,1,0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$K)
    sigma <- runif(2)
    beta <- rnorm(Data$K*(Data$L-1), 0, sigma[1])
    gamma <- rnorm(Data$K*(Data$M-1), 0, sigma[2])
    U <- rwishartc(Data$K+1, Data$S)
    return(c(alpha, beta, gamma, U[upper.tri(U, diag=TRUE)], sigma))
    }
Data <- list(J=J, K=K, L=L, M=M, N=N, PGF=PGF, S=S, X=X, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.beta=pos.beta, pos.gamma=pos.gamma, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- matrix(c(parm[Data$pos.beta], rep(0,Data$K)), Data$K, Data$L)
    beta[,Data$L] <- -rowSums(beta[,-Data$L])
    gamma <- matrix(c(parm[Data$pos.gamma],
        rep(0,Data$K)), Data$K, Data$M)
    gamma[,Data$M] <- -rowSums(gamma[,-Data$M])
    U <- as.parm.matrix(U, Data$K, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnorm(beta, 0, sigma[1], log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma[2], log=TRUE))
    U.prior <- dwishartc(U, Data$K+1, Data$S, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(0,Data$N,Data$K)
    for (k in 1:Data$K) {
        mu[,k] <- alpha[k] + beta[k,Data$X[,1]] + gamma[k,Data$X[,2]]}
    LL <- sum(dmvnpc(Data$Y, mu, U, log=TRUE))
    ### Variance Components, Omnibus
    s.o.beta <- sd(as.vector(beta))
    s.o.gamma <- sd(as.vector(gamma))
    s.o.epsilon <- sd(as.vector(Data$Y - mu))
    ### Variance Components, Univariate
    s.beta <- sqrt(apply(beta, 1, var))
    s.gamma <- sqrt(apply(gamma, 1, var))
    s.epsilon <- sqrt(apply(Data$Y - mu, 2, var))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + U.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, s.o.beta, s.o.gamma,
        s.o.epsilon, s.beta, s.gamma, s.epsilon),
        yhat=rmvnpc(nrow(mu), mu, U), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K), rep(0,K*(L-1)), rep(0,K*(M-1)),
    upper.triangle(S, diag=TRUE), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: intercepts
for (k in 1:K)
    cat("alpha[", k, "] -- true:", true.alpha[k],
        " post. mean:", round(fit$Summary2[pos.alpha[k], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multivariate Binary Probit

The multivariate binary probit extends the univariate binary probit (see
section [link](#binary.probit)) to \\J\\ correlated binary outcomes
using a latent multivariate normal representation. The latent variables
\\\textbf{W}\\ are constrained to match the observed binary pattern, and
the correlation structure \\\rho\\ captures the dependence between
outcomes. The multivariate binary probit model with correlated latent
utilities was developed by Chib and Greenberg [\[25\]](#ref25) using
Gibbs sampling. Multivariate probit models are applied in health
economics to jointly model the adoption of multiple preventive health
behaviours.

#### Form

\\\textbf{W}\_{i,1:J} \sim \mathcal{N}\_J(\mu\_{i,1:J}, \Omega^{-1}),
\quad i=1,\dots,N\\ \\\textbf{W}\_{i,j} \in \left\\ \begin{array}{l l}
\$\[0,10\]\$ & \quad \mbox{if \$\textbf{y}\_i = j\$}\\ \$\[-10,0\]\$ \\
\end{array} \right. \\ \\\mu\_{1:N,j} = \textbf{X} \beta\_{j,1:K}\\
\\\Omega = \rho^{-1}\\ \\\beta\_{j,k} \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,(J-1), \quad k=1,\dots,K\\ \\\beta\_{J,k} = -
\sum^{J-1}\_{j=1} \beta\_{j,k}\\ \\\rho \sim \mathcal{U}(-1, 1)\\

#### model_spec() notation

The multivariate binary probit cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires latent variable augmentation with truncation
constraints tied to the observed binary outcomes, plus a correlation
matrix prior. These constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42195)
N <- 30
J <- 2   # Number of binary dependent variables
K <- 3   # Number of columns in design matrix X

#### True parameter values
true.beta <- matrix(c(0.5, -0.8, 0.3,
                      -0.4, 0.6, -0.2), J, K, byrow = TRUE)
true.rho <- 0.4

#### Generate design matrix and latent response
X <- cbind(1, matrix(rnorm(N * (K - 1), 0, 1), N, K - 1))
mu.true <- tcrossprod(X, true.beta)
rho.mat <- diag(J)
rho.mat[1, 2] <- rho.mat[2, 1] <- true.rho
Omega.true <- solve(rho.mat)
U.true <- chol(Omega.true)
W <- interval(rmvnpc(N, mu.true, U.true), -10, 10)
Y <- 1 * (W >= 0)
apply(Y, 2, table)

#### Assemble Data list
mon.names <- "LP"
rho <- matrix(NA, J, J)
rho[upper.tri(rho)] <- 0
parm.names <- as.parm.names(list(beta=matrix(0,J,K), rho=rho,
    W=matrix(0,N,J)))
pos.beta <- grep("beta", parm.names)
pos.rho <- grep("rho", parm.names)
pos.W <- grep("W", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J*Data$K)
    rho <- rep(0, length(which(upper.tri(diag(Data$J)))))
    W <- matrix(runif(Data$N*Data$J,-10,0), Data$N, Data$J)
    W <- ifelse(Y == 1, abs(W), W)
    return(c(beta, rho, as.vector(W)))}
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta,
    pos.rho=pos.rho, pos.W=pos.W)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$J, Data$K)
    u <- interval(parm[Data$pos.rho], -1, 1)
    rho <- diag(Data$J)
    rho[upper.tri(rho)] <- u
    rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
    if(is.positive.semidefinite(rho) == FALSE)
        rho <- as.positive.semidefinite(rho)
    parm[Data$pos.rho] <- upper.triangle(rho)
    Omega <- as.inverse(rho)
    U <- chol(Omega)
    W <- matrix(parm[Data$pos.W], Data$N, Data$J)
    W[Data$Y == 0] <- interval(W[Data$Y == 0], -10, 0)
    W[Data$Y == 1] <- interval(W[Data$Y == 1], 0, 10)
    parm[Data$pos.W] <- as.vector(W)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    rho.prior <- sum(dunif(u, -1, 1, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, beta)
    LL <- sum(dmvnpc(W, mu, U, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + rho.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=1*(rmvnpc(nrow(mu), mu, U) >= 0), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- GIV(Model, Data, PGF=TRUE)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: regression coefficients
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], J, K)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Parameter recovery: correlation
post.rho <- fit$Summary2[pos.rho, "Mean"]
cat("rho -- true:", true.rho,
    " post. mean:", round(post.rho, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
ground_truth[pos.rho]  <- true.rho[seq_along(pos.rho)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multivariate Laplace Regression

The multivariate Laplace regression models \\K\\ response variables
jointly through a multivariate Laplace likelihood, which provides
heavier tails than the multivariate normal and therefore greater
robustness to outliers. The precision matrix receives a Wishart prior,
and each response column is linearly related to the design matrix
through its own coefficient vector. Multivariate Laplace regression
extends the asymmetric Laplace likelihood to vector-valued responses,
following Kozumi and Kobayashi [\[67\]](#ref67). Multivariate Laplace
regression is used in finance for joint quantile modelling of correlated
asset returns during tail events.

#### Form

\\\textbf{Y}\_{i,k} \sim \mathcal{L}\_K(\mu\_{i,k}, \Sigma), \quad
i=1,\dots,N; \quad k=1,\dots,K\\ \\\mu\_{i,k} = \textbf{X}\_{1:N,k}
\beta\_{k,1:J}\\ \\\Sigma = \Omega^{-1}\\ \\\Omega \sim
\mathcal{W}\_{K+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_K\\
\\\beta\_{k,j} \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\

#### model_spec() notation

The multivariate Laplace distribution is available in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
registry as `MVLaplaceC(mu, U)`, but this model requires the Wishart
prior on the Cholesky factor of the precision matrix and matrix-valued
beta parameters indexed per response, which cannot be expressed in the
declarative DSL. The manual `Model` function below is recommended.

#### Ground truth and data

``` r

set.seed(42196)
N <- 50    # Number of records
J <- 4     # Number of columns in design matrix
K <- 2     # Number of DVs

#### True parameter values
true.beta <- matrix(c(15.0, 0.02, -0.03, -3.0,
                        5.0, 0.01, -0.01, -1.5), K, J, byrow = TRUE)

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 50), N, J - 1))
mu.true <- tcrossprod(X, true.beta)
Y <- mu.true + matrix(rnorm(N * K, 0, 2), N, K)
S <- diag(K)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,K,J), U=diag(K)),
    uppertri=c(0,1))
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K*Data$J)
    U <- rwishartc(Data$K+1, Data$S)
    return(c(beta, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, S=S, X=X, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$K, Data$J)
    U <- as.parm.matrix(U, Data$K, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    U.prior <- dwishartc(U, Data$K+1, Data$S, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, beta)
    LL <- sum(dmvlc(Data$Y, mu, U, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rmvlc(nrow(mu), mu, U), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J*K), upper.triangle(S, diag=TRUE))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], K, J)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multivariate Poisson Regression

The multivariate Poisson regression models \\K\\ correlated count
responses by introducing observation-level multivariate normal random
effects \\\gamma\\ that induce dependence across the Poisson margins.
Each response is log-linearly related to the design matrix, and the
random effects receive a Wishart prior on their precision matrix.
Multivariate Poisson regression models correlated count data through
shared latent factors, following the approach of Chib and Winkelmann
[\[26\]](#ref26). Multivariate Poisson models are applied in sports
analytics to jointly model home and away goal counts in football
matches.

#### Form

\\\textbf{Y}\_{i,k} \sim \mathcal{P}(\lambda\_{i,k}), \quad i=1,\dots,N
\quad k=1,\dots,K\\ \\\lambda\_{i,k} =
\exp(\textbf{X}\_{i,k}\beta\_{k,1:J} + \gamma\_{i,k}), \quad
i=1,\dots,N, \quad k=1,\dots,K\\ \\\beta\_{k,j} \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J, \quad k=1,\dots,K\\ \\\gamma\_{i,1:K} \sim
\mathcal{N}\_K(0, \Omega^{-1}), \quad i=1,\dots,N\\ \\\Omega \sim
\mathcal{W}\_{K+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_K\\

#### model_spec() notation

The multivariate Poisson regression cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires observation-level multivariate normal random
effects linked through the log rate, combined with a Wishart prior on
the precision matrix. These constructs fall outside the declarative
scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42197)
N <- 20   # Number of records
J <- 4    # Number of columns in design matrix
K <- 3    # Number of DVs

#### True parameter values
true.beta <- matrix(c(0.5, 0.3, -0.2, 0.1,
                      0.8, -0.1, 0.4, -0.3,
                      0.3, 0.2, 0.1, -0.1), K, J, byrow = TRUE)

#### Generate design matrix and response
X <- matrix(runif(N * J), N, J); X[,1] <- 1
Omega <- matrix(c(2.0, 0.5, 0.3,
                  0.5, 1.5, 0.2,
                  0.3, 0.2, 1.8), K, K)
gamma <- rmvnp(N, 0, Omega)
Y <- matrix(rpois(N * K, exp(tcrossprod(X, true.beta) + gamma)),
    N, K)
S <- diag(K)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,K,J), gamma=matrix(0,N,K),
    U=S), uppertri=c(0,0,1))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K*Data$J)
    gamma <- rnorm(Data$N*Data$K)
    U <- rwishartc(Data$K+1, Data$S)
    return(c(beta, gamma, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, S=S, X=X, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$K, Data$J)
    gamma <- matrix(parm[Data$pos.gamma], Data$N, Data$K)
    U <- as.parm.matrix(U, Data$K, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dmvnpc(gamma, 0, U, log=TRUE))
    U.prior <- dwishartc(U, Data$K+1, Data$S, log=TRUE)
    ### Log-Likelihood
    lambda <- exp(tcrossprod(Data$X, beta) + gamma)
    LL <- sum(dpois(Data$Y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(prod(dim(lambda)), lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K*J), rep(0,N*K), rep(0,K*(K+1)/2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], K, J)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Multivariate Regression

The multivariate regression models \\K\\ continuous response variables
jointly through a multivariate normal likelihood, using the Huang-Wand
prior on the covariance matrix rather than a Wishart prior. The
Huang-Wand prior decomposes the covariance into a correlation matrix and
marginal variances, providing more interpretable shrinkage on the
correlation structure. Multivariate normal regression with an
inverse-Wishart prior on the covariance was formalised by Zellner
[\[121\]](#ref121) in the seemingly unrelated regression framework.
Multivariate regression is used in nutritional epidemiology to jointly
model multiple dietary biomarkers as functions of food intake variables.

#### Form

\\\textbf{Y}\_{i,k} \sim \mathcal{N}\_K(\mu\_{i,k}, \Sigma), \quad
i=1,\dots,N; \quad k=1,\dots,K\\ \\\mu\_{i,k} = \textbf{X}\_{1:N,k}
\beta\_{k,1:J}\\ \\\Sigma \sim \mathcal{HW}\_{2}(\gamma, 1e6)\\
\\\beta\_{k,j} \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J, \quad
k=1,\dots,K\\

#### model_spec() notation

The multivariate regression uses a Huang-Wand prior on the Cholesky
factor of the covariance matrix, which cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The manual `Model` function below is recommended.

#### Ground truth and data

``` r

set.seed(42198)
N <- 50    # Number of records
J <- 4     # Number of columns in design matrix
K <- 2     # Number of DVs

#### True parameter values
true.beta <- matrix(c(18.0, 0.01, -0.04, -2.5,
                        6.0, 0.005, -0.01, -1.0), K, J, byrow = TRUE)

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1), 0, 50), N, J - 1))
mu.true <- tcrossprod(X, true.beta)
Y <- mu.true + matrix(rnorm(N * K, 0, 2), N, K)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,K,J), gamma=rep(0,K),
    U=diag(K)), uppertri=c(0,0,1))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K*Data$J)
    gamma <- runif(Data$K)
    U <- rhuangwandc(2, gamma, rep(1,Data$K))
    return(c(beta, gamma, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$K, Data$J)
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.gamma] <- gamma
    U <- as.parm.matrix(U, Data$K, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    HW.prior <- dhuangwandc(U, 2, gamma, rep(1e6,Data$K), log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, beta)
    LL <- sum(dmvnc(Data$Y, mu, U, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + HW.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rmvnc(nrow(mu), mu, U), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J*K), rep(1,K), rep(0,K*(K+1)/2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], K, J)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Seemingly Unrelated Regression (SUR)

The seemingly unrelated regression models multiple equations that appear
unrelated but share correlated errors, introduced by Zellner
[\[120\]](#ref120). Each equation has its own regressors and
coefficients, but the cross-equation error correlation is captured by a
joint multivariate normal likelihood. This implementation uses the
Yang-Berger prior for the precision matrix, which is a reference prior
for multivariate normal models that provides less informative shrinkage
than the Wishart. The data are Grunfeld’s investment equations for
General Electric and Westinghouse. Seemingly unrelated regression was
introduced by Zellner [\[120\]](#ref120) to exploit cross-equation error
correlations for efficiency gains. SUR models are used in economics to
jointly estimate demand equations for multiple goods, where common
unobserved factors create correlated residuals.

#### Form

\\\textbf{Y}\_{t,k} \sim \mathcal{N}\_K(\mu\_{t,k}, \Omega^{-1}), \quad
t=1,\dots,T, \quad k=1,\dots,K\\ \\\mu\_{1,t} = \alpha_1 + \alpha_2
\textbf{X}\_{t-1,1} + \alpha_3 \textbf{X}\_{t-1,2}, \quad t=2,\dots,T\\
\\\mu\_{2,t} = \beta_1 + \beta_2 \textbf{X}\_{t-1,3} + \beta_3
\textbf{X}\_{t-1,4}, \quad t=2,\dots,T\\ \\\alpha_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ where \\J=3\\, \\K=2\\, and \\T=20\\.

#### model_spec() notation

The SUR model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires equation-specific regressors with a shared
multivariate normal error structure and the Yang-Berger reference prior
on the precision matrix. These constructs fall outside the declarative
scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

This example uses the historical Grunfeld investment data as originally
presented by Zellner [\[120\]](#ref120), so no simulated ground truth is
provided. The dependent variables are gross investment for General
Electric and Westinghouse, and the regressors are lagged values of firm
value and capital stock.

``` r

set.seed(42199)
T <- 20 #Time-periods
year <- c(1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,1945,1946,
    1947,1948,1949,1950,1951,1952,1953,1954)
IG <- c(33.1,45.0,77.2,44.6,48.1,74.4,113.0,91.9,61.3,56.8,93.6,159.9,
    147.2,146.3,98.3,93.5,135.2,157.3,179.5,189.6)
VG <- c(1170.6,2015.8,2803.3,2039.7,2256.2,2132.2,1834.1,1588.0,1749.4,
    1687.2,2007.7,2208.3,1656.7,1604.4,1431.8,1610.5,1819.4,2079.7,
    2371.6,2759.9)
CG <- c(97.8,104.4,118.0,156.2,172.6,186.6,220.9,287.8,319.9,321.3,319.6,
    346.0,456.4,543.4,618.3,647.4,671.3,726.1,800.3,888.9)
IW <- c(12.93,25.90,35.05,22.89,18.84,28.57,48.51,43.34,37.02,37.81,
    39.27,53.46,55.56,49.56,32.04,32.24,54.38,71.78,90.08,68.60)
VW <- c(191.5,516.0,729.0,560.4,519.9,628.5,537.1,561.2,617.2,626.7,
    737.2,760.5,581.4,662.3,583.8,635.2,723.8,864.1,1193.5,1188.9)
CW <- c(1.8,0.8,7.4,18.1,23.5,26.5,36.2,60.8,84.4,91.2,92.4,86.0,111.1,
    130.6,141.8,136.7,129.7,145.5,174.8,213.5)
J <- 2 #Number of dependent variables
Y <- matrix(c(IG,IW), T, J)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,3), beta=rep(0,3),
    U=diag(J)), uppertri=c(0,0,1))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(3)
    beta <- rnorm(3)
    U <- runif(Data$J*(Data$J+1)/2)
    return(c(alpha, beta, U))
    }
Data <- list(J=J, PGF=PGF, T=T, Y=Y, CG=CG, CW=CW, IG=IG, IW=IW,
    VG=VG, VW=VW, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.beta=pos.beta)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    U <- as.parm.matrix(U, Data$J, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    U.prior <- dyangbergerc(U, log=TRUE)
    ### Log-Likelihood
    mu <- Data$Y
    mu[-1,1] <- alpha[1] + alpha[2]*Data$CG[-Data$T] +
        alpha[3]*Data$VG[-Data$T]
    mu[-1,2] <- beta[1] + beta[2]*Data$CW[-Data$T] +
        beta[3]*Data$VW[-Data$T]
    LL <- sum(dmvnpc(Data$Y[-1,], mu[-1,], U, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rmvnpc(nrow(mu), mu, U), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,3), rep(0,3), rep(0,J*(J+1)/2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter estimates (no ground truth for historical data)
post.alpha <- fit$Summary2[pos.alpha, "Mean"]
post.beta  <- fit$Summary2[pos.beta, "Mean"]
cat("GE equation (alpha):", round(post.alpha, 3), "\n")
cat("WH equation (beta): ", round(post.beta, 3), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Simultaneous Equations

Simultaneous equation models, pioneered by Haavelmo [\[54\]](#ref54)
(Nobel Prize 1989) for handling endogeneity through structural
identification, are fundamental in macroeconomics for modelling the
joint determination of output, interest rates, and inflation. This
example uses Klein’s Model I [\[64\]](#ref64) regarding economic
fluctuations in the United States in 1920-1941 (\\N\\=22). Usually, this
example is modeled with 3-stage least squares (3SLS), excluding the
uncertainty from multiple stages. By constraining each element in the
instrumental variables matrix \\\nu \in \[-10,10\]\\, this example
estimates the model without resorting to stages. The dependent variable
is matrix \\\textbf{Y}\\, in which \\\textbf{Y}\_{1,1:N}\\ is
\\\textbf{C}\\ (consumption), \\\textbf{Y}\_{2,1:N}\\ is \\\textbf{I}\\
(investment), and \\\textbf{Y}\_{3,1:N}\\ is \\\textbf{Wp}\\ (private
wages). Here is a data dictionary:

``` r
    A = Time Trend measured as years from 1931
    C = Consumption
    G = Government Nonwage Spending
    I = Investment
    K = Capital Stock
    P = Private (Corporate) Profits
    T = Indirect Business Taxes Plus Neg Exports
    Wg = Government Wage Bill
    Wp = Private Wages
    X = Equilibrium Demand (GNP)
```

See Klein [\[64\]](#ref64) for more information. Simultaneous equation
models were pioneered by Haavelmo [\[54\]](#ref54) (Nobel Prize 1989)
for handling endogeneity through structural identification. Simultaneous
equations are fundamental in macroeconomics for modelling the joint
determination of output, interest rates, and inflation.

#### Form

\\\textbf{Y} \sim \mathcal{N}\_3(\mu, \Omega^{-1})\\ \\ \mu\_{1,1} =
\alpha_1 + \alpha_2 \nu\_{1,1} + \alpha_4 \nu\_{2,1}\\ \\ \mu\_{1,i} =
\alpha_1 + \alpha_2 \nu\_{1,i} + \alpha_3 \textbf{P}\_{i-1} + \alpha_4
\nu\_{2,i}, \quad i=2,\dots,N\\ \\ \mu\_{2,1} = \beta_1 + \beta_2
\nu\_{1,1} + \beta_4 \textbf{K}\_1\\ \\ \mu\_{2,i} = \beta_1 + \beta_2
\nu\_{1,i} + \beta_3 \textbf{P}\_{i-1} + \beta_4 \textbf{K}\_i, \quad
i=2,\dots,N\\ \\\mu\_{3,1} = \gamma_1 + \gamma_2 \nu\_{3,1} + \gamma_4
\textbf{A}\_1\\ \\\mu\_{3,i} = \gamma_1 + \gamma_2 \nu\_{3,i} + \gamma_3
\textbf{X}\_{i-1} + \gamma_4 \textbf{A}\_i, \quad i=2,\dots,N\\
\\\textbf{Z}\_{j,i} \sim \mathcal{N}(\nu\_{j,i}, \sigma^2_j), \quad
j=1,\dots,3\\ \\\nu\_{j,1} = \pi\_{j,1} + \pi\_{j,3} \textbf{K}\_1 +
\pi\_{j,5} \textbf{A}\_1 + \pi\_{j,6} \textbf{T}\_1 + \pi\_{j,7}
\textbf{G}\_1, \quad j=1,\dots,3\\ \\\nu\_{j,i} = \pi\_{j,1} +
\pi\_{j,2} \textbf{P}\_{i-1} + \pi\_{j,3} \textbf{K}\_i + \pi\_{j,4}
\textbf{X}\_{i-1} + \pi\_{j,5} \textbf{A}\_i + \pi\_{j,6}
\textbf{T}\_i + \pi\_{j,7} \textbf{G}\_i, \quad i=1,\dots,N, \quad
j=1,\dots,3\\ \\\alpha_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,4\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,4\\ \\\gamma_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,4\\ \\\pi\_{j,i} \sim
\mathcal{N}(0, 1000) \in \[-10,10\], \quad j=1,\dots,3, \quad
i=1,\dots,N\\ \\\sigma_j \sim \mathcal{HC}(25), \quad j=1,\dots,3\\
\\\Omega \sim \mathcal{W}\_4(\textbf{S}), \quad \textbf{S} =
\textbf{I}\_3\\

#### model_spec() notation

The simultaneous equations model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires instrumental variable equations, lagged
endogenous variables, bounded latent parameters, and a multivariate
normal structural equation system with a Wishart prior on the precision
matrix. These constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

This example uses Klein’s historical data for the US economy 1920-1941,
so no simulated ground truth is provided.

``` r

set.seed(42200)
N <- 22
A <- c(-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
C <- c(39.8,41.9,45,49.2,50.6,52.6,55.1,56.2,57.3,57.8,55,50.9,45.6,46.5,
    48.7,51.3,57.7,58.7,57.5,61.6,65,69.7)
G <- c(2.4,3.9,3.2,2.8,3.5,3.3,3.3,4,4.2,4.1,5.2,5.9,4.9,3.7,4,4.4,2.9,4.3,
    5.3,6.6,7.4,13.8)
I <- c(2.7,-0.2,1.9,5.2,3,5.1,5.6,4.2,3,5.1,1,-3.4,-6.2,-5.1,-3,-1.3,2.1,2,
    -1.9,1.3,3.3,4.9)
K <- c(180.1,182.8,182.6,184.5,189.7,192.7,197.8,203.4,207.6,210.6,215.7,
    216.7,213.3,207.1,202,199,197.7,199.8,201.8,199.9,201.2,204.5)
P <- c(12.7,12.4,16.9,18.4,19.4,20.1,19.6,19.8,21.1,21.7,15.6,11.4,7,11.2,
    12.3,14,17.6,17.3,15.3,19,21.1,23.5)
T <- c(3.4,7.7,3.9,4.7,3.8,5.5,7,6.7,4.2,4,7.7,7.5,8.3,5.4,6.8,7.2,8.3,6.7,
    7.4,8.9,9.6,11.6)
Wg <- c(2.2,2.7,2.9,2.9,3.1,3.2,3.3,3.6,3.7,4,4.2,4.8,5.3,5.6,6,6.1,7.4,
    6.7,7.7,7.8,8,8.5)
Wp <- c(28.8,25.5,29.3,34.1,33.9,35.4,37.4,37.9,39.2,41.3,37.9,34.5,29,28.5,
    30.6,33.2,36.8,41,38.2,41.6,45,53.3)
X <- c(44.9,45.6,50.1,57.2,57.1,61,64,64.4,64.5,67,61.2,53.4,44.3,45.1,
    49.7,54.4,62.7,65,60.9,69.5,75.7,88.4)
year <- c(1920,1921,1922,1923,1924,1925,1926,1927,1928,1929,1930,1931,1932,
    1933,1934,1935,1936,1937,1938,1939,1940,1941)
Y <- matrix(c(C,I,Wp), 3, N, byrow=TRUE)
Z <- matrix(c(P, Wp+Wg, X), 3, N, byrow=TRUE)
S <- diag(nrow(Y))

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,4), beta=rep(0,4),
    gamma=rep(0,4), pi=matrix(0,3,7), sigma=rep(0,3),
    U=diag(3)), uppertri=c(0,0,0,0,0,1))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.pi <- grep("pi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(4)
    beta <- rnorm(4)
    gamma <- rnorm(4)
    pi <- rnorm(3*7)
    sigma <- runif(3)
    U <- rwishartc(ncol(Data$Y)+1, Data$S)
    return(c(alpha, beta, gamma, pi, sigma, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(A=A, C=C, G=G, I=I, K=K, N=N, P=P, PGF=PGF, S=S, T=T, Wg=Wg,
    Wp=Wp, X=X, Y=Y, Z=Z, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.pi=pos.pi, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    parm[Data$pos.pi] <- pi <- interval(parm[Data$pos.pi], -10, 10)
    pi <- matrix(pi, 3, 7)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    U <- as.parm.matrix(U, nrow(Data$S), parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    pi.prior <- sum(dnormv(pi, 0, 1000, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    U.prior <- dwishartc(U, nrow(Data$S)+1, Data$S, log=TRUE)
    ### Log-Likelihood
    mu <- nu <- matrix(0,3,Data$N)
    for (i in 1:3) {
        nu[i,1] <- pi[i,1] + pi[i,3]*Data$K[1] + pi[i,5]*Data$A[1] +
    pi[i,6]*Data$T[1] + pi[i,7]*Data$G[1]
        nu[i,-1] <- pi[i,1] + pi[i,2]*Data$P[-Data$N] +
    pi[i,3]*Data$K[-1] + pi[i,4]*Data$X[-Data$N] +
    pi[i,5]*Data$A[-1] + pi[i,6]*Data$T[-1] +
    pi[i,7]*Data$G[-1]}
    LL <- sum(dnorm(Data$Z, nu, matrix(sigma, 3, Data$N), log=TRUE))
    mu[1,1] <- alpha[1] + alpha[2]*nu[1,1] + alpha[4]*nu[2,1]
    mu[1,-1] <- alpha[1] + alpha[2]*nu[1,-1] +
        alpha[3]*Data$P[-Data$N] + alpha[4]*nu[2,-1]
    mu[2,1] <- beta[1] + beta[2]*nu[1,1] + beta[4]*Data$K[1]
    mu[2,-1] <- beta[1] + beta[2]*nu[1,-1] +
        beta[3]*Data$P[-Data$N] + beta[4]*Data$K[-1]
    mu[3,1] <- gamma[1] + gamma[2]*nu[3,1] + gamma[4]*Data$A[1]
    mu[3,-1] <- gamma[1] + gamma[2]*nu[3,-1] +
        gamma[3]*Data$X[-Data$N] + gamma[4]*Data$A[-1]
    LL <- LL + sum(dmvnpc(t(Data$Y), t(mu), U, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + pi.prior +
        sigma.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=t(rmvnp(ncol(mu), t(mu), U)), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,4), rep(0,4), rep(0,4), rep(0,3*7), rep(1,3),
    upper.triangle(S, diag=TRUE))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter estimates (no ground truth for historical data)
post.alpha <- fit$Summary2[pos.alpha, "Mean"]
post.beta  <- fit$Summary2[pos.beta, "Mean"]
post.gamma <- fit$Summary2[pos.gamma, "Mean"]
cat("Consumption eq (alpha):", round(post.alpha, 3), "\n")
cat("Investment eq (beta):  ", round(post.beta, 3), "\n")
cat("Wages eq (gamma):      ", round(post.gamma, 3), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Hierarchical and multilevel models

### Linear Regression, Hierarchical Bayesian

Hierarchical (or multilevel) priors allow the regression coefficients
\\\beta_j\\ to share information through a common normal distribution
whose mean \\\gamma\\ and variance \\\delta\\ are themselves given
hyperpriors. An additional layer of hierarchy places a half-Cauchy
hyperprior on the scale \\\tau\\ of the residual standard deviation
prior. This structure induces partial pooling across coefficients,
shrinking noisy estimates toward the group mean while still allowing
each \\\beta_j\\ to adapt to the data. Hierarchical Bayesian regression
was developed by Lindley and Smith [\[72\]](#ref72) and extended by
Gelman and Hill [\[47\]](#ref47) for multilevel data structures.
Hierarchical models are standard in education research for estimating
school-level effects while pooling information across schools.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_j \sim \mathcal{N}(\gamma, \delta), \quad j=1,\dots,J\\ \\\gamma
\sim \mathcal{N}(0, 1000)\\ \\\delta \sim \mathcal{HC}(25)\\ \\\sigma
\sim \mathcal{HC}(\tau)\\ \\\tau \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the DSL does not currently support hierarchical hyperpriors
where prior parameters (e.g., \\\gamma\\, \\\delta\\, \\\tau\\) are
themselves given distributions. The nested prior structure \\\beta_j
\sim \mathcal{N}(\gamma, \delta)\\ with \\\gamma \sim
\mathcal{N}(0,1000)\\ and \\\delta \sim \mathcal{HC}(25)\\ requires the
manual `Model` function.

#### Ground truth and data

We simulate data with four predictors (including intercept) where the
true coefficients are drawn from a common normal distribution with mean
\\\gamma = 1.0\\ and standard deviation \\\delta = 0.8\\, reflecting the
hierarchical structure. The residual standard deviation is 0.6 and the
hyperprior scale is \\\tau = 2.0\\.

``` r

set.seed(42150)
N <- 200
J <- 4
true.gamma <- 1.0
true.delta <- 0.8
true.sigma <- 0.6
true.tau <- 2.0
true.beta <- rnorm(J, true.gamma, true.delta)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y <- rnorm(N, X %*% true.beta, true.sigma)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), gamma=0, delta=0, sigma=0,
    tau=0))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    gamma <- rnorm(1)
    delta <- runif(1)
    sigma <- runif(1)
    tau <- runif(1)
    return(c(beta, gamma, delta, sigma, tau))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.delta=pos.delta, pos.sigma=pos.sigma, pos.tau=pos.tau, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    gamma <- parm[Data$pos.gamma]
    delta <- interval(parm[Data$pos.delta], 1e-100, Inf)
    parm[Data$pos.delta] <- delta
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperprior
    gamma.prior <- dnormv(gamma, 0, 1000, log=TRUE)
    delta.prior <- dhalfcauchy(delta, 25, log=TRUE)
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, gamma, delta, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, tau, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + delta.prior + sigma.prior +
        tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), 0, rep(1, 3))
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare estimates to ground truth
cat("True beta:     ", round(true.beta, 3), "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True gamma:    ", true.gamma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.gamma]), 3), "\n")
cat("True delta:    ", true.delta, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.delta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.tau]   <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression, Multilevel

The multilevel (or random-coefficients) linear regression allows each
group \\g\\ to have its own coefficient vector \\\beta_g\\, drawn from a
common multivariate normal distribution with mean \\\gamma\\ and
precision matrix \\\Omega\\. This induces partial pooling across groups:
groups with sparse data borrow strength from the population-level mean,
while groups with abundant data retain their own estimates. The
precision matrix receives a Wishart prior, and the population-level
means receive vague normal priors. Multilevel (mixed effects) models
were unified in a Bayesian framework by Gelman and Hill
[\[47\]](#ref47), allowing partial pooling of group-level parameters.
Multilevel regression is applied in political science to estimate
state-level opinion from national survey data via multilevel regression
and poststratification (MRP).

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_i = \textbf{X}
\beta\_{\textbf{m}\[i\],1:J}\\ \\\beta\_{g,1:J} \sim
\mathcal{N}\_J(\gamma, \Omega^{-1}), \quad g=1,\dots,G\\ \\\Omega \sim
\mathcal{W}\_{J+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_J\\
\\\gamma_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\ where \\\textbf{m}\\ is a vector of length \\N\\, and
each element indicates the multilevel group (\\g=1,\dots,G\\) for the
associated record.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it requires group-indexed coefficient matrices \\\beta\_{g,j}\\,
Wishart-distributed precision matrices parameterized via a Cholesky
factor, and row-wise indexing through a group membership vector. These
structural features exceed the current capabilities of the declarative
DSL.

#### Ground truth and data

We simulate data for \\G = 3\\ groups with an intercept and one
predictor (\\J = 2\\). The population-level coefficients are \\\gamma =
(1.5, -0.8)\\, and the group-specific coefficients are drawn from a
multivariate normal with a known covariance matrix. The residual
standard deviation is 0.3.

``` r

set.seed(42151)
N <- 90
J <- 2
G <- 3
true.gamma <- c(1.5, -0.8)
true.Sigma <- matrix(c(0.5, 0.1, 0.1, 0.3), J, J)
true.Sigma <- as.positive.definite(true.Sigma)
true.sigma <- 0.3
true.beta <- matrix(NA, G, J)
for (g in 1:G) true.beta[g,] <- rmvn(1, true.gamma, true.Sigma)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
m <- rcat(N, rep(1/G, G))
y <- rowSums(true.beta[m, ] * X) + rnorm(N, 0, true.sigma)
S <- diag(J)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=matrix(0,G,J), gamma=rep(0,J),
    sigma=0, U=S), uppertri=c(0,0,0,1))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    U <- rwishartc(Data$J+1, Data$S)
    gamma <- rnorm(Data$J)
    beta <- as.vector(rmvnpc(Data$G, gamma, U))
    sigma <- runif(1)
    return(c(beta, gamma, sigma, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(G=G, J=J, N=N, PGF=PGF, S=S, X=X, m=m, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$G, Data$J)
    gamma <- parm[Data$pos.gamma]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    U <- as.parm.matrix(U, Data$J, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    U.prior <- dwishartc(U, Data$J+1, Data$S, log=TRUE)
    beta.prior <- sum(dmvnpc(beta, gamma, U, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 100, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- rowSums(beta[Data$m,] * Data$X)
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + U.prior + beta.prior + gamma.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, G * J), rep(0, J), 1,
    upper.triangle(S, diag=TRUE))
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare group-level coefficients to ground truth
cat("True beta (group coefficients):\n")
print(round(true.beta, 3))
cat("\nPosterior mean beta:\n")
print(matrix(round(colMeans(fit$Posterior1[, pos.beta]), 3), G, J))
cat("\nTrue gamma:    ", round(true.gamma, 3), "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.gamma]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[grep("Sigma", parm.names)] <- true.Sigma[seq_along(grep("Sigma", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Nested Hierarchical (3-Level) **\[NEW\]**

The three-level nested hierarchical model accommodates observations
clustered within sub-groups that are themselves nested within
higher-level groups. A grand mean \\\mu\\ anchors the population,
group-level deviations \\b_j\\ capture between-group heterogeneity, and
sub-group deviations \\c\_{jk}\\ capture within-group variation, with
residual noise \\\sigma^2\\ at the observation level. The three variance
components \\\sigma_b\\, \\\sigma_c\\, and \\\sigma\\ are estimated
simultaneously, partitioning total variation across the hierarchy. The
nested multilevel model follows the framework of Gelman and Hill
[\[47\]](#ref47), with comprehensive treatment by Raudenbush and Bryk
[\[135\]](#ref135) in the frequentist tradition. Three-level models are
applied in education research to analyse student outcomes nested within
classrooms within schools.

#### Form

\\y\_{ijk} \sim \mathcal{N}(\mu + b_j + c\_{jk}, \sigma^2)\\ \\b_j \sim
\mathcal{N}(0, \sigma_b^2), \quad j=1,\dots,M_1\\ \\c\_{jk} \sim
\mathcal{N}(0, \sigma_c^2), \quad k=1,\dots,M_2\\ \\\mu \sim
\mathcal{N}(0, 1000)\\ \\\sigma \sim \mathcal{HC}(25)\\ \\\sigma_b \sim
\mathcal{HC}(25)\\ \\\sigma_c \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because nested random effects with multiple variance components at
different hierarchical levels require explicit index mappings and
separate distributional assumptions for each grouping factor, which
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42222)
M1 <- 5
M2 <- 3
n.per <- 10
N <- M1 * M2 * n.per
true.mu <- 10
true.sigma <- 1
true.sigma.b <- 2
true.sigma.c <- 0.5

#### Simulate hierarchical data
true.b <- rnorm(M1, 0, true.sigma.b)
true.c <- rnorm(M1 * M2, 0, true.sigma.c)

#### Build index vectors
group.idx <- rep(rep(1:M1, each=M2 * n.per))
subgroup.idx <- rep(rep(1:(M1 * M2), each=n.per))
y <- true.mu + true.b[group.idx] + true.c[subgroup.idx] +
    rnorm(N, 0, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, b=rep(0, M1),
    c=rep(0, M1 * M2), sigma=0, sigma.b=0, sigma.c=0))
pos.mu <- grep("^mu$", parm.names)
pos.b <- grep("^b\\[", parm.names)
pos.c <- grep("^c\\[", parm.names)
pos.sigma <- grep("^sigma$", parm.names)
pos.sigma.b <- grep("^sigma\\.b$", parm.names)
pos.sigma.c <- grep("^sigma\\.c$", parm.names)
PGF <- function(Data) {
    sigma.b <- runif(1, 0.5, 3)
    sigma.c <- runif(1, 0.1, 1)
    return(c(rnorm(1), rnorm(Data$M1, 0, sigma.b),
        rnorm(Data$M1 * Data$M2, 0, sigma.c),
        runif(1), runif(1), runif(1)))
    }
Data <- list(M1=M1, M2=M2, N=N, PGF=PGF, group.idx=group.idx,
    mon.names=mon.names, parm.names=parm.names, pos.b=pos.b,
    pos.c=pos.c, pos.mu=pos.mu, pos.sigma=pos.sigma,
    pos.sigma.b=pos.sigma.b, pos.sigma.c=pos.sigma.c,
    subgroup.idx=subgroup.idx, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    b <- parm[Data$pos.b]
    cc <- parm[Data$pos.c]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    sigma.b <- interval(parm[Data$pos.sigma.b], 1e-100, Inf)
    parm[Data$pos.sigma.b] <- sigma.b
    sigma.c <- interval(parm[Data$pos.sigma.c], 1e-100, Inf)
    parm[Data$pos.sigma.c] <- sigma.c
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 1000, log=TRUE)
    b.prior <- sum(dnorm(b, 0, sigma.b, log=TRUE))
    c.prior <- sum(dnorm(cc, 0, sigma.c, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    sigma.b.prior <- dhalfcauchy(sigma.b, 25, log=TRUE)
    sigma.c.prior <- dhalfcauchy(sigma.c, 25, log=TRUE)
    ### Log-Likelihood
    mu.vec <- mu + b[Data$group.idx] + cc[Data$subgroup.idx]
    LL <- sum(dnorm(Data$y, mu.vec, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + b.prior + c.prior + sigma.prior +
        sigma.b.prior + sigma.c.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu.vec, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(mean(y), rep(0, M1), rep(0, M1 * M2), 1, 1, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.b]     <- true.b[seq_along(pos.b)]
ground_truth[pos.c]     <- true.c[seq_along(pos.c)]
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Normal, Multilevel

This is Gelman’s eight schools example [\[48\]](#ref48), a canonical
hierarchical model where \\J=8\\ school-level treatment effects share a
common population distribution. The group-level standard deviations are
known, and the model estimates school-specific effects \\\theta_j\\, a
population mean \\\theta\_\mu\\, and the between-school standard
deviation \\\theta\_\sigma\\. Uniform priors are used instead of the
default half-Cauchy priors, because the latter overwhelm the likelihood
in this small-sample setting. The Bayesian multilevel normal model
follows the hierarchical framework of Gelman and Hill [\[47\]](#ref47),
with partial pooling of group means. Multilevel normal models are
applied in meta-analysis to estimate an overall treatment effect while
accounting for between-study heterogeneity.

#### Form

\\\textbf{y}\_j \sim \mathcal{N}(\theta_j, \sigma^2_j), \quad
j=1,\dots,J\\ \\\theta_j \sim \mathcal{N}(\theta\_{\mu},
\theta\_\sigma^2)\\ \\\theta\_{\mu} \sim \mathcal{N}(0, 1000000)\\
\\\theta\_{\sigma\[j\]} \sim \mathcal{N}(0, 1000)\\ \\\sigma \sim
\mathcal{U}(0, 1000)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the known, observation-specific standard deviations \\\sigma_j\\
enter the likelihood as fixed constants rather than parameters, and the
hierarchical centring uses a parametric form with mixed prior families
(normal precision prior on \\\theta\_\mu\\, uniform on
\\\theta\_\sigma\\) that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42164)

#### Eight schools data (Gelman et al., 2004)
J <- 8
y <- c(28.4, 7.9, -2.8, 6.8, -0.6, 0.6, 18.0, 12.2)
sd <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(theta=rep(0,J), theta.mu=0,
    theta.sigma=0))
pos.theta <- 1:J
pos.theta.mu <- J+1
pos.theta.sigma <- J+2
PGF <- function(Data) {
    theta.mu <- rnorm(1)
    theta.sigma <- runif(1)
    theta <- rnorm(Data$J, theta.mu, theta.sigma)
    return(c(theta, theta.mu, theta.sigma))
    }
Data <- list(J=J, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.theta=pos.theta, pos.theta.mu=pos.theta.mu,
    pos.theta.sigma=pos.theta.sigma, sd=sd, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    theta.mu <- parm[Data$pos.theta.mu]
    theta.sigma <- interval(parm[Data$pos.theta.sigma], 1e-100, Inf)
    parm[Data$pos.theta.sigma] <- theta.sigma
    ### Parameters
    theta <- parm[Data$pos.theta]
    ### Log-Hyperprior
    theta.mu.prior <- dnormp(theta.mu, 0, 1.0E-6, log=TRUE)
    theta.sigma.prior <- dunif(theta.sigma, 0, 1000, log=TRUE)
    ### Log-Prior
    theta.prior <- sum(dnorm(theta, theta.mu, theta.sigma, log=TRUE))
    sigma.prior <- sum(dunif(Data$sd, 0, 1000, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnorm(Data$y, theta, Data$sd, log=TRUE))
    ### Log-Posterior
    LP <- LL + theta.prior + theta.mu.prior + theta.sigma.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(theta), theta, Data$sd), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(mean(y),J), mean(y), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Missing data models

### Linear Regression with Full Missingness

In the full-missingness setting, both the response \\\textbf{y}\\ and
one or more columns of the design matrix \\\textbf{X}\\ contain missing
values. The full-likelihood approach augments the state with auxiliary
parameters \\\gamma\\ (imputed covariate values) and \\\delta\\ (imputed
response values), which are sampled jointly with the regression
coefficients. Matrix \\\alpha\\ captures inter-predictor regression
effects used to model the distribution of the covariates themselves.
This approach is sound as long as the model remains identifiable; when
it is not, a preliminary imputation stage via the `MISS` function is
preferable. Bayesian treatment of missing data as latent variables
follows the framework of Little and Rubin [\[74\]](#ref74) and Gelman et
al. [\[48\]](#ref48). Full-missingness models are applied in
longitudinal clinical studies where both covariates and responses may be
missing due to patient dropout.

#### Form

\\\textbf{y}^{imp} \sim \mathcal{N}(\nu, \sigma^2_J)\\
\\\textbf{X}^{imp} \sim \mathcal{N}(\mu, \sigma^2\_{-J})\\ \\\nu =
\textbf{X}^{imp} \beta\\ \\\mu = \textbf{X}^{imp} \alpha\\
\\\textbf{y}^{imp} = \left\\ \begin{array}{l l} \\\\ & \quad \mbox{if
\$\textbf{y}^{mis}\$}\\ \textbf{y}^{obs} \\ \end{array} \right. \\
\\\textbf{X}^{imp} = \left\\ \begin{array}{l l} \\\\ & \quad \mbox{if
\$\textbf{X}^{mis}\$}\\ \textbf{X}^{obs} \\ \end{array} \right. \\
\\\alpha\_{j,l} \sim \mathcal{N}(0, 1000), \quad j=1,\dots,(J-1), \quad
l=1,\dots,(J-1)\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\gamma_m \sim \mathcal{N}(0, 1000), \quad m=1,\dots,M\\
\\\delta_p \sim \mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\sigma_j
\sim \mathcal{HC}(25), \quad j=1,\dots,J\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because missing-data augmentation requires auxiliary parameter vectors
(\\\gamma\\, \\\delta\\) whose dimensions depend on the missingness
pattern of the observed data. The DSL cannot express conditional
imputation logic or dynamic-length parameter vectors.

#### Ground truth and data

We simulate a dataset with \\N = 100\\ observations and \\J = 5\\
columns (intercept plus four predictors). Approximately 10% of covariate
entries and 5% of response values are set to `NA` after generation, so
the model must recover both the regression coefficients and the missing
values.

``` r

set.seed(42152)
N <- 100
J <- 5
true.beta <- c(2.0, -1.2, 0.5, 1.8, -0.6)
true.sigma <- rep(0.5, J)
X.complete <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y.complete <- as.vector(X.complete %*% true.beta + rnorm(N, 0, true.sigma[J]))
### Introduce missingness in covariates (~10% per column, never in intercept)
X <- X.complete
for (j in 2:J) {
    miss.idx <- sample(1:N, round(N * 0.10))
    X[miss.idx, j] <- NA
}
### Introduce missingness in response (~5%)
y <- y.complete
y[sample(1:N, round(N * 0.05))] <- NA
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=matrix(0,J-1,J-1),
    beta=rep(0,J),
    gamma=rep(0,sum(is.na(X))),
    delta=rep(0,sum(is.na(y))),
    sigma=rep(0,J)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm((Data$J-1)*(Data$J-1))
    beta <- rnorm(Data$J)
    gamma <- rnorm(sum(is.na(Data$X)))
    delta <- rnorm(sum(is.na(Data$y)), mean(Data$y, na.rm=TRUE), 1)
    sigma <- runif(Data$J)
    return(c(alpha, beta, gamma, delta, sigma))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.gamma=pos.gamma, pos.delta=pos.delta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- matrix(parm[Data$pos.alpha], Data$J-1, Data$J-1)
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    delta <- parm[Data$pos.delta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    delta.prior <- sum(dnormv(delta, 0, 1000, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- X.imputed <- Data$X
    X.imputed[which(is.na(X.imputed))] <- gamma
    y.imputed <- Data$y
    y.imputed[which(is.na(y.imputed))] <- delta
    for (j in 2:Data$J) {mu[,j] <- tcrossprod(X.imputed[,-j],
        t(alpha[,(j-1)]))}
    nu <- tcrossprod(X.imputed, t(beta))
    LL <- sum(dnorm(X.imputed[,-1], mu[,-1],
        matrix(sigma[1:(Data$J-1)], Data$N, Data$J-1), log=TRUE),
        dnorm(y.imputed, nu, sigma[Data$J], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + gamma.prior + delta.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(nu), nu, sigma[Data$J]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, (J - 1)^2), rep(0, J), rep(0, sum(is.na(X))),
    rep(0, sum(is.na(y))), rep(1, J))
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare regression coefficients to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma[J]: ", true.sigma[J], "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma[J]]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression with Missing Response

This example introduces missing-value handling via data augmentation
with auxiliary variables. The response vector contains both observed
values \\\textbf{y}^{obs}\\ and missing values \\\textbf{y}^{mis}\\. The
auxiliary parameter vector \\\alpha\\ holds imputed values for the
missing responses; by including these in the state, the sampler explores
the full-likelihood surface and produces valid posterior inference
despite the incomplete data. In the model form, \\M\\ denotes the number
of missing values, though in the data list it serves as a binary
indicator vector. Missing response imputation within MCMC was developed
by Tanner and Wong [\[102\]](#ref102) and formalised in the data
augmentation framework of Little and Rubin [\[74\]](#ref74). Missing
response models are common in survey research where some respondents
fail to answer the outcome question.

#### Form

\\\textbf{y}^{imp} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\textbf{y}^{imp}
= \left\\ \begin{array}{l l} \\\\ & \quad \mbox{if
\$\textbf{y}^{mis}\$}\\ \textbf{y}^{obs} \\ \end{array} \right. \\ \\\mu
= \textbf{X}\beta\\ \\\alpha_m \sim \mathcal{N}(0, 1000), \quad
m=1,\dots,M\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\
\\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the auxiliary imputation parameters \\\alpha\\ are
dynamic-length vectors determined by the missingness pattern, and the
conditional replacement logic for \\\textbf{y}^{imp}\\ cannot be
expressed in the declarative DSL.

#### Ground truth and data

We simulate data with known regression coefficients and then introduce
5% missingness in the response. The design matrix has an intercept and
three predictors (\\J = 4\\), with true coefficients \\\beta = (3.0,
-1.5, 0.7, 1.2)\\ and \\\sigma = 0.8\\.

``` r

set.seed(42153)
N <- 150
J <- 4
true.beta <- c(3.0, -1.5, 0.7, 1.2)
true.sigma <- 0.8
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y.complete <- rnorm(N, X %*% true.beta, true.sigma)
y <- y.complete
y[sample(1:N, round(N * 0.05))] <- NA
M <- ifelse(is.na(y), 1, 0)
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,sum(M)), beta=rep(0,J),
    sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(sum(Data$M), mean(Data$y, na.rm=TRUE), 1)
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(alpha, beta, sigma))
    }
Data <- list(J=J, M=M, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    y.imputed <- Data$y
    y.imputed[which(is.na(Data$y))] <- alpha
    LL <- sum(dnorm(y.imputed, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, sum(M)), rep(0, J), 1)
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare estimates to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Linear Regression with Missing Response via ABB

The Approximate Bayesian Bootstrap (ABB), using the `ABB` function,
imputes missing values in the response given a propensity score
estimated from the covariates. Vector \\\alpha\\ parameterizes the
propensity score model, vector \\\beta\\ captures the regression
effects, and vector \\\gamma\\ (monitored, not sampled directly)
contains the imputed missing values generated by ABB within each MCMC
iteration. This approach is useful when the missingness mechanism
depends on observed covariates. The Approximate Bayesian Bootstrap for
missing data was proposed by Rubin and Schenker [\[96\]](#ref96) as a
computationally efficient imputation mechanism. ABB-based imputation is
used in large-scale federal surveys where fully Bayesian approaches are
computationally infeasible.

#### Form

\\\textbf{y}^{imp} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\textbf{y}^{imp}
= \left\\ \begin{array}{l l} \\\\ & \quad \mbox{if
\$\textbf{y}^{mis}\$}\\ \textbf{y}^{obs} \\ \end{array} \right. \\ \\\mu
= \textbf{X}\beta\\ \\\gamma \sim p(\textbf{y}^{obs} \|
\textbf{y}^{obs}, \textbf{y}^{mis}, \eta)\\ \\\eta = \frac{1}{1 +
\exp(-\nu)}\\ \\\nu = \textbf{X} \alpha\\ \\\alpha_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because ABB imputation is a stochastic within-iteration procedure that
modifies the response vector using propensity-score stratification. The
DSL cannot express mid-likelihood stochastic imputation steps or the
[`ABB()`](https://robustecologies.github.io/lucifer/reference/ABB.md)
function call.

#### Ground truth and data

We simulate data with an intercept and three predictors (\\J = 4\\) and
introduce 5% missingness in the response. The propensity-score model
uses the same design matrix as the regression, so \\\alpha\\ and
\\\beta\\ each have \\J\\ elements. The monitored quantities include the
imputed missing values \\\gamma\\.

``` r

set.seed(42154)
N <- 150
J <- 4
true.beta <- c(2.5, -0.8, 1.3, -0.5)
true.sigma <- 0.9
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
y.complete <- rnorm(N, X %*% true.beta, true.sigma)
y <- y.complete
y[sample(1:N, round(N * 0.05))] <- NA
M <- ifelse(is.na(y), 1, 0)
mon.names <- c("LP", paste("gamma[", 1:sum(is.na(y)), "]", sep=""))
parm.names <- as.parm.names(list(alpha=rep(0,J), beta=rep(0,J), sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(alpha, beta, sigma))
    }
Data <- list(J=J, M=M, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    y.imputed <- Data$y
    mu <- tcrossprod(Data$X, t(beta))
    nu <- as.vector(tcrossprod(Data$X, t(alpha)))
    eta <- invlogit(nu)
    breaks <- as.vector(quantile(eta, probs=c(0,0.2,0.4,0.6,0.8,1)))
    B <- matrix(breaks[-length(breaks)], length(Data$y), 5, byrow=TRUE)
    z <- rowSums(eta >= B)
    for (i in 1:5) {
        if(any(is.na(Data$y[which(z == i)]))) {
    imp <- unlist(ABB(Data$y[which(z == i)]))
    y.imputed[which({z == i} & is.na(Data$y))] <- imp}}
    gamma <- y.imputed[which(is.na(Data$y))]
    LL <- sum(dbern(Data$M, eta, log=TRUE),
        dnorm(y.imputed, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,gamma),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), rep(0, J), 1)
```

#### Fitting and recovery

``` r

### Prescribe the best algorithm for this model
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### fit with lucifer using the top-ranked MCMC method
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

### Posterior recovery: compare regression coefficients to ground truth
cat("True beta:     ", true.beta, "\n")
cat("Posterior mean:", round(colMeans(fit$Posterior1[, pos.beta]), 3), "\n")
cat("True sigma:    ", true.sigma, "\n")
cat("Posterior mean:", round(mean(fit$Posterior1[, pos.sigma]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Revision, Normal

This example provides both an analytic solution and numerical
approximation of the revision of a normal distribution. Given a normal
prior distribution (\\\alpha\\) and data distribution (\\\beta\\), the
posterior (\\\gamma\\) is the revised normal distribution. This is an
introductory example of Bayesian inference, and allows the user to
experiment with numerical approximation, such as with MCMC in `lucifer`.
Note that, regardless of the data sample size \\N\\ in this example,
Laplace Approximation is inappropriate due to asymptotics since the data
(\\\beta\\) is perceived by the algorithm as a single datum rather than
a collection of data. MCMC, on the other hand, is biased only by the
effective number of samples taken of the posterior. Bayesian belief
revision updates a normal prior with new data to produce a posterior,
illustrating the core mechanism of sequential learning (Bernardo and
Smith [\[9\]](#ref9)). Sequential updating is used in quality control
for continuously revising the estimated mean of a manufacturing process
as new measurements arrive.

``` r

set.seed(42209)
#### Analytic Solution
prior.mu <- 0
prior.sigma <- 10
N <- 10
data.mu <- 1
data.sigma <- 2
posterior.mu <- (prior.sigma^ -2 * prior.mu + N * data.sigma^ -2 * data.mu) /
    (prior.sigma^ -2 + N * data.sigma^ -2)
posterior.sigma <- sqrt(1/(prior.sigma^ -2 + data.sigma^ -2))
posterior.mu
posterior.sigma
```

#### Form

\\\alpha \sim \mathcal{N}(0,10)\\ \\\beta \sim \mathcal{N}(1,2)\\
\\\gamma = \frac{\alpha^{-2}\_\sigma \alpha + N \beta^{-2}\_\sigma
\beta}{\alpha^{-2}\_\sigma + N \beta^{-2}\_\sigma}\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it computes a derived posterior quantity \\\gamma\\ from two
separate prior-like distributions rather than specifying a standard
likelihood-prior structure, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42209)
N <- 10
mon.names <- c("LP","gamma")
parm.names <- c("alpha","beta")
PGF <- function(Data) {
    alpha <- rnorm(1,0,10)
    beta <- rnorm(1,1,2)
    return(c(alpha, beta))
    }
Data <- list(N=N, PGF=PGF, mon.names=mon.names, parm.names=parm.names)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    alpha.mu <- 0
    alpha.sigma <- 10
    beta.mu <- 1
    beta.sigma <- 2
    ### Parameters
    alpha <- parm[1]
    beta <- parm[2]
    ### Log-Prior
    alpha.prior <- dnorm(alpha, alpha.mu, alpha.sigma, log=TRUE)
    ### Log-Likelihood
    LL <- dnorm(beta, beta.mu, beta.sigma, log=TRUE)
    ### Posterior
    gamma <- (alpha.sigma^ -2 * alpha + N * beta.sigma^ -2 * beta) /
        (alpha.sigma^ -2 + N * beta.sigma^ -2)
    ### Log-Posterior
    LP <- LL + alpha.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,gamma),
        yhat=rnorm(1, beta.mu, beta.sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0,0)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Compare analytic vs. MCMC posterior
cat("Analytic posterior mu:   ", round(posterior.mu, 4), "\n")
cat("MCMC posterior mean(gamma):", round(mean(fit$Monitor[, 2]), 4), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Time series: autoregressive and moving average

### AR(p)

The autoregressive model of order \\p\\ predicts the current value of a
time series as a linear combination of the \\p\\ most recent lagged
values plus a constant; the lag structure is specified through the
vector \\L\\ of lag positions. The residual variance \\\sigma^2\\ is
modelled with a half-Cauchy prior. The autoregressive model was
formalised by Yule [\[117\]](#ref117) and developed within the
Box-Jenkins framework (Box and Jenkins [\[16\]](#ref16)). AR models are
fundamental in macroeconomics for forecasting GDP growth and inflation
from lagged quarterly observations.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p}, \quad
t=1,\dots,T\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi_p \sim
\mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The AR(p) model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean function involves recursive time indexing
(\\\mu_t\\ depends on lagged values \\\textbf{y}\_{t-p}\\), which
requires imperative loop constructs that fall outside the declarative
scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42106)
T <- 300
L <- c(1, 5)          # Autoregressive lags
P <- length(L)         # Autoregressive order

#### True parameter values
true.alpha <- 0.1
true.phi <- c(0.6, -0.2)
true.sigma <- 0.5

#### Simulate AR(p) process from the true model
y <- numeric(T)
y[1:max(L)] <- rnorm(max(L), true.alpha, true.sigma)
for (t in (max(L) + 1):T) {
    y[t] <- true.alpha
    for (p in 1:P) y[t] <- y[t] + true.phi[p] * y[t - L[p]]
    y[t] <- y[t] + rnorm(1, 0, true.sigma)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    sigma <- rhalfcauchy(1,5)
    return(c(alpha, phi, sigma))
    }
Data <- list(L=L, PGF=PGF, P=P, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L[p])] <- mu[-c(1:Data$L[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L[p])]
    LL <- sum(dnorm(Data$y[-c(1:Data$L[Data$P])], mu[-c(1:Data$L[Data$P])],
        sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,P+1), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### ARIMAX **\[NEW\]**

The ARIMAX model extends the classical ARMA framework by incorporating
exogenous regressors alongside autoregressive and moving-average
components, allowing the conditional mean to depend on both the process
history and external covariates. The moving-average term introduces
residual recursion: at each time step the one-step-ahead prediction
error feeds back into the mean, making the likelihood a sequential
computation that cannot be vectorised in closed form. The ARMA framework
was established by Box and Jenkins [\[16\]](#ref16), and the extension
to exogenous inputs was formalised by Pankratz [\[150\]](#ref150) under
the transfer function paradigm. ARIMAX models are widely used in energy
demand forecasting, where electricity consumption depends on its own
momentum, weather covariates, and economic indicators.

#### Form

\\z_t \sim \mathcal{N}(\mu_t, \sigma^2)\\ \\\mu_t = \alpha +
\sum\_{i=1}^{p} \phi_i z\_{t-i} + \sum\_{j=1}^{q} \theta_j
\epsilon\_{t-j} + \textbf{x}\_t^\top \boldsymbol{\beta}\\ \\\epsilon_t =
z_t - \mu_t\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi \sim
\mathcal{N}(0, 1000), \quad \|\phi\| \< 1\\ \\\theta \sim \mathcal{N}(0,
1000)\\ \\\beta_k \sim \mathcal{N}(0, 1000)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the moving-average residual recursion \\\epsilon\_{t-j} =
z\_{t-j} - \mu\_{t-j}\\ requires an imperative loop where each time step
depends on the fitted residual from the previous step. This kind of
sequential state update falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42202)
T <- 200
p <- 1   # AR order
q <- 1   # MA order
k <- 2   # exogenous variables

#### True parameter values
true.alpha <- 0.5
true.phi <- 0.6
true.theta.ma <- 0.3
true.beta <- c(2, 1)
true.sigma <- 1.0

#### Exogenous covariates
x <- matrix(rnorm(T * k), nrow=T, ncol=k)

#### Simulate ARIMAX process
z <- numeric(T)
eps <- numeric(T)
z[1] <- rnorm(1, true.alpha + x[1,] %*% true.beta, true.sigma)
eps[1] <- 0
for (t in 2:T) {
    mu.t <- true.alpha + true.phi * z[t-1] +
        true.theta.ma * eps[t-1] + x[t,] %*% true.beta
    eps[t] <- rnorm(1, 0, true.sigma)
    z[t] <- mu.t + eps[t]
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=0, theta.ma=0,
    beta=rep(0,k), sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.theta.ma <- grep("theta.ma", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(1, -0.9, 0.9)
    theta.ma <- rnorm(1)
    beta <- rnorm(Data$k)
    sigma <- rhalfcauchy(1, 5)
    return(c(alpha, phi, theta.ma, beta, sigma))
    }
Data <- list(T=T, p=p, q=q, k=k, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.theta.ma=pos.theta.ma, pos.beta=pos.beta, pos.sigma=pos.sigma,
    x=x, z=z)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- interval(parm[Data$pos.phi], -0.999, 0.999)
    parm[Data$pos.phi] <- phi
    theta.ma <- parm[Data$pos.theta.ma]
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- dnormv(phi, 0, 1000, log=TRUE)
    theta.ma.prior <- dnormv(theta.ma, 0, 1000, log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- numeric(Data$T)
    eps <- numeric(Data$T)
    mu[1] <- alpha + Data$x[1,] %*% beta
    eps[1] <- Data$z[1] - mu[1]
    for (t in 2:Data$T) {
        mu[t] <- alpha + phi * Data$z[t-1] +
            theta.ma * eps[t-1] + Data$x[t,] %*% beta
        eps[t] <- Data$z[t] - mu[t]
    }
    LL <- sum(dnorm(Data$z[-1], mu[-1], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + theta.ma.prior +
        beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$T, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0.3, 0, rep(0, k), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("phi -- true:", true.phi,
    " post. mean:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")
cat("theta.ma -- true:", true.theta.ma,
    " post. mean:", round(fit$Summary2[pos.theta.ma, "Mean"], 3), "\n")
for (j in 1:k)
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(fit$Summary2[pos.beta[j], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Autoregressive Moving Average, ARMA(p,q)

The autoregressive moving average model extends the AR(p) by adding
\\q\\ lagged residual terms to the mean equation, enabling more flexible
modelling of both the autoregressive and moving-average components of
stationary time series. The moving-average coefficients \\\psi_q\\
capture short-run dynamics that pure autoregressive models require many
more lags to approximate. The ARMA model was formalised by Box and
Jenkins [\[16\]](#ref16), combining autoregressive and moving average
components for flexible short-memory time series modelling. ARMA models
are used in hydrology to forecast river discharge from historical
streamflow measurements.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p} +
\sum^Q\_{q=1} \theta_q \epsilon\_{t-q}\\ \\\epsilon_t = \textbf{y}\_t -
\mu_t\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi_p \sim
\mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\sigma \sim
\mathcal{HC}(25)\\ \\\theta_q \sim \mathcal{N}(0, 1000), \quad
q=1,\dots,Q\\

#### model_spec() notation

The ARMA(p,q) model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean equation involves recursive time indexing
(\\\mu_t\\ depends on \\\textbf{y}\_{t-p}\\ and on lagged residuals
\\\epsilon\_{t-q} = \textbf{y}\_{t-q} - \mu\_{t-q}\\), which requires
imperative loop constructs that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42113)
T <- 260
L.P <- c(1, 5)        # Autoregressive lags
L.Q <- c(1, 2)        # Moving average lags
P <- length(L.P)      # Autoregressive order
Q <- length(L.Q)      # Moving average order

#### True parameter values
true.alpha <- 0.05
true.phi   <- c(0.50, -0.15)
true.theta <- c(0.30, -0.10)
true.sigma <- 0.4

#### Simulate ARMA(p,q) process
y <- numeric(T)
eps <- numeric(T)
y[1:max(L.P)] <- rnorm(max(L.P), true.alpha, true.sigma)
eps[1:max(L.P)] <- rnorm(max(L.P), 0, true.sigma)
for (t in (max(L.P) + 1):T) {
    mu.t <- true.alpha
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L.P[p]]
    for (q in 1:Q)
        if (t > L.Q[q]) mu.t <- mu.t + true.theta[q] * eps[t - L.Q[q]]
    eps[t] <- rnorm(1, 0, true.sigma)
    y[t] <- mu.t + eps[t]
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), sigma=0,
    theta=rep(0,Q)))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    sigma <- rhalfcauchy(1,5)
    theta <- rnorm(Data$Q)
    return(c(alpha, phi, sigma, theta))
    }
Data <- list(L.P=L.P, L.Q=L.Q, PGF=PGF, P=P, Q=Q, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.sigma=pos.sigma, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    theta <- parm[Data$pos.theta]
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    theta.prior <- sum(dnormv(theta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L.P[p])] <- mu[-c(1:Data$L.P[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L.P[p])]
    epsilon <- Data$y - mu
    for (q in 1:Data$Q)
        mu[-c(1:Data$L.Q[q])] <- mu[-c(1:Data$L.Q[q])] +
    theta[q]*epsilon[1:(Data$T-Data$L.Q[q])]
    LL <- sum(dnorm(Data$y[-c(1:Data$L.P[Data$P])],
        mu[-c(1:Data$L.P[Data$P])], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + sigma.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,P), 1, rep(0,Q))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
for (q in 1:Q)
    cat("theta[", q, "] -- true:", true.theta[q],
        " post. mean:", round(fit$Summary2[pos.theta[q], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Beta AR(1) **\[NEW\]**

The Beta AR(1) model captures serial dependence in a time series of
proportions by linking the conditional mean of a Beta distribution to
the logit of the previous observation. The precision parameter \\\phi\\
controls how tightly each observation concentrates around its
conditional mean, while the autoregressive coefficient \\\psi\\ governs
persistence on the logit scale. This formulation was proposed by Rocha
and Cribari-Neto [\[149\]](#ref149) as an extension of Beta regression
to dynamic settings where temporal correlation cannot be ignored. Beta
AR(1) models are applied in hydrology to forecast reservoir storage
ratios, which are bounded proportions exhibiting strong temporal
persistence.

#### Form

\\y_t \sim \text{Beta}(\mu_t \phi, (1 - \mu_t) \phi)\\
\\\text{logit}(\mu_t) = \alpha + \psi \\ \text{logit}(y\_{t-1})\\
\\\alpha \sim \mathcal{N}(0, 1000)\\ \\\psi \sim \mathcal{N}(0, 1000)\\
\\\phi \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the conditional mean depends on \\\text{logit}(y\_{t-1})\\,
the logit-transformed lagged observation. This autoregressive structure
requires recursive indexing into the observed data vector, which falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42203)
T <- 200

#### True parameter values
true.alpha <- 0.2
true.psi <- 0.5
true.phi <- 20

#### Simulate Beta AR(1) process
y <- numeric(T)
y[1] <- 0.5
for (t in 2:T) {
    logit.mu <- true.alpha + true.psi * qlogis(y[t-1])
    mu.t <- plogis(logit.mu)
    a <- mu.t * true.phi
    b <- (1 - mu.t) * true.phi
    y[t] <- rbeta(1, a, b)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, psi=0, phi=0))
pos.alpha <- grep("alpha", parm.names)
pos.psi <- grep("psi", parm.names)
pos.phi <- grep("phi", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1, 0, 0.5)
    psi <- runif(1, -0.9, 0.9)
    phi <- rhalfcauchy(1, 10)
    return(c(alpha, psi, phi))
    }
Data <- list(T=T, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.psi=pos.psi, pos.phi=pos.phi, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    psi <- parm[Data$pos.psi]
    phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
    parm[Data$pos.phi] <- phi
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    psi.prior <- dnormv(psi, 0, 1000, log=TRUE)
    phi.prior <- dhalfcauchy(phi, 25, log=TRUE)
    ### Log-Likelihood
    y.clamped <- pmin(pmax(Data$y, 0.001), 0.999)
    logit.mu <- alpha + psi * qlogis(y.clamped[1:(Data$T-1)])
    mu <- plogis(logit.mu)
    a <- mu * phi
    b <- (1 - mu) * phi
    LL <- sum(dbeta(y.clamped[2:Data$T], a, b, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + psi.prior + phi.prior
    mu.all <- c(plogis(alpha), mu)
    a.all <- mu.all * phi
    b.all <- (1 - mu.all) * phi
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbeta(Data$T, a.all, b.all), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0.3, 10)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("psi -- true:", true.psi,
    " post. mean:", round(fit$Summary2[pos.psi, "Mean"], 3), "\n")
cat("phi -- true:", true.phi,
    " post. mean:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.psi]   <- true.psi[seq_along(pos.psi)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Dirichlet Autoregressive **\[NEW\]**

The Dirichlet autoregressive model captures temporal dependence in
compositional time series by letting the concentration parameters of the
Dirichlet distribution at each time step depend on the log-transformed
composition at the previous step through an autoregressive coefficient
\\\psi\\, while an overall concentration \\\phi\\ controls the
dispersion of each observation around its conditional mean. The
log-linear link ensures positivity of the Dirichlet parameters, and the
softmax transformation maps the linear predictor to the simplex before
scaling by \\\phi\\. Compositional time series analysis using log-ratio
transformations was formalised by Aitchison [\[160\]](#ref160), whose
framework underpins this autoregressive extension. Dirichlet
autoregressive models are applied in ecology to track temporal changes
in species relative abundance within communities.

#### Form

\\\textbf{y}\_t \sim \mathcal{D}(\alpha_t), \quad t=2,\dots,T\\
\\\alpha\_{t,k} = \phi \cdot \eta\_{t,k}, \quad k=1,\dots,R\\
\\\eta\_{t,k} = \frac{\exp(\mu_k + \psi \log y\_{t-1,k})}{\sum\_{r=1}^R
\exp(\mu_r + \psi \log y\_{t-1,r})}\\ \\\mu_k \sim \mathcal{N}(0, 100),
\quad k=1,\dots,R\\ \\\psi \sim \mathcal{N}(0, 100)\\ \\\phi \sim
\mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the autoregressive structure conditions each time step’s
Dirichlet parameters on the observed composition at the previous step,
requiring a sequential likelihood computation over lagged compositional
observations that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42209)
R <- 3    # Compositional parts
T <- 150  # Time steps

#### True parameter values
true.mu <- c(1.0, 0.5, 0.2)
true.psi <- 0.4
true.phi <- 50

#### Simulate compositional time series
y <- matrix(NA, T, R)
y[1,] <- c(0.40, 0.35, 0.25)
for (t in 2:T) {
    log.eta <- true.mu + true.psi * log(y[t-1,])
    eta <- exp(log.eta) / sum(exp(log.eta))
    alpha <- true.phi * eta
    y[t,] <- rdirichlet(1, alpha)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=rep(0,R), psi=0, phi=0))
pos.mu <- grep("mu", parm.names)
pos.psi <- grep("psi", parm.names)
pos.phi <- grep("phi", parm.names)
PGF <- function(Data) {
    mu <- rnorm(Data$R)
    psi <- rnorm(1, 0, 0.5)
    phi <- runif(1, 5, 100)
    return(c(mu, psi, phi))
    }
Data <- list(R=R, T=T, y=y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.psi=pos.psi,
    pos.phi=pos.phi)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    psi <- parm[Data$pos.psi]
    phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
    parm[Data$pos.phi] <- phi
    ### Log-Prior
    mu.prior <- sum(dnormv(mu, 0, 100, log=TRUE))
    psi.prior <- dnormv(psi, 0, 100, log=TRUE)
    phi.prior <- dhalfcauchy(phi, 25, log=TRUE)
    ### Log-Likelihood
    LL <- 0
    alpha.last <- rep(1, Data$R)
    for (t in 2:Data$T) {
        log.eta <- mu + psi * log(Data$y[t-1,])
        eta <- exp(log.eta) / sum(exp(log.eta))
        alpha.t <- phi * eta
        LL <- LL + ddirichlet(Data$y[t,], alpha.t, log=TRUE)
        if (t == Data$T) alpha.last <- alpha.t
    }
    ### Log-Posterior
    LP <- LL + mu.prior + psi.prior + phi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rdirichlet(Data$T - 1, alpha.last)[,1], parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0.5, R), 0.3, 30)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True mu:  ", true.mu, "\n")
cat("Post. mean:", round(fit$Summary2[pos.mu, "Mean"], 3), "\n")
cat("True psi: ", true.psi, "\n")
cat("Post. mean:", round(fit$Summary2[pos.psi, "Mean"], 3), "\n")
cat("True phi: ", true.phi, "\n")
cat("Post. mean:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]  <- true.mu[seq_along(pos.mu)]
ground_truth[pos.phi] <- true.phi[seq_along(pos.phi)]
ground_truth[pos.psi] <- true.psi[seq_along(pos.psi)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Distributed Lag, Koyck

The Koyck distributed lag model applies geometrically decaying weights
\\\beta\lambda^l\\ to the effects of discrete events in a covariate, so
that the impact of each event diminishes exponentially with lag \\l\\.
The autoregressive parameter \\\phi\\ captures persistence in the
dependent variable itself, and the decay rate \\\lambda \in (0,1)\\
governs how quickly event effects fade. The Koyck [\[66\]](#ref66)
distributed lag model imposes geometric decay on lag coefficients,
reducing an infinite lag structure to two parameters. Koyck models are
used in advertising research to estimate how the effect of marketing
expenditure decays over subsequent sales periods.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_t = \alpha + \phi
\textbf{y}\_{t-1} + \sum^K\_{k=1} \textbf{X}\_{t,k} \beta
\lambda^{\textbf{L}\[t,k\]}, \quad k=1,\dots,K, \quad t=2,\dots,T\\
\\\mu_1 = \alpha + \sum^K\_{k=1} \textbf{X}\_{1,k} \beta
\lambda^{\textbf{L}\[1,k\]}, \quad k=1,\dots,K\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\beta \sim \mathcal{N}(0, 1000)\\ \\\lambda
\sim \mathcal{U}(0, 1)\\ \\\phi \sim \mathcal{N}(0, 1000)\\ \\\sigma
\sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Koyck lag structure involves constructing event-specific lag
matrices with geometric decay weights applied inside a summation, which
requires imperative matrix construction logic that falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42207)
T <- 260

#### True parameter values
true.alpha  <- 0.001
true.phi    <- 0.3
true.beta   <- 0.5
true.lambda <- 0.7
true.sigma  <- 0.02

#### Simulate a return-like series with event effects
y <- numeric(T)
x <- rbinom(T, 1, 0.05)  # sparse events
K <- length(which(x != 0))
L <- X <- matrix(0, T, K)
for (i in 1:K) {
    X[which(x != 0)[i]:T, i] <- x[which(x != 0)[i]]
    L[(which(x != 0)[i]):T, i] <- 0:(T - which(x != 0)[i])
}
y[1] <- true.alpha + sum(X[1,] * true.beta * true.lambda^L[1,]) +
    rnorm(1, 0, true.sigma)
for (t in 2:T) {
    y[t] <- true.alpha + true.phi * y[t-1] +
        sum(X[t,] * true.beta * true.lambda^L[t,]) +
        rnorm(1, 0, true.sigma)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, beta=0, lambda=0, phi=0, sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(1)
    lambda <- runif(1)
    phi <- rnorm(1)
    sigma <- runif(1)
    return(c(alpha, beta, lambda, phi, sigma))
    }
Data <- list(L=L, PGF=PGF, T=T, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.lambda=pos.lambda, pos.phi=pos.phi, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    lambda <- interval(parm[Data$pos.lambda], 0, 1)
    parm[Data$pos.lambda] <- lambda
    phi <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- dnormv(beta, 0, 1000, log=TRUE)
    lambda.prior <- dunif(lambda, 0, 1, log=TRUE)
    phi.prior <- dnormv(phi, 0, 1000, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- c(alpha, alpha + phi*Data$y[-Data$T]) +
        rowSums(Data$X * beta * lambda^ Data$L)
    LL <- sum(dnorm(Data$y[-1], mu[-1], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + lambda.prior + phi.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,2), 0.5, 0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True alpha:  ", true.alpha, " post. mean:", round(fit$Summary2[1, "Mean"], 4), "\n")
cat("True beta:   ", true.beta,  " post. mean:", round(fit$Summary2[2, "Mean"], 4), "\n")
cat("True lambda: ", true.lambda," post. mean:", round(fit$Summary2[3, "Mean"], 4), "\n")
cat("True phi:    ", true.phi,   " post. mean:", round(fit$Summary2[4, "Mean"], 4), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha]  <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.phi]    <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Exponential Smoothing

Exponential smoothing produces forecasts as a weighted average of past
observations with exponentially declining weights controlled by the
smoothing parameter \\\alpha \in (0,1)\\; the Bayesian formulation
places a uniform prior on \\\alpha\\ and a half-Cauchy prior on the
observation noise. Exponential smoothing methods were unified under a
state space framework by Hyndman et al. [\[58\]](#ref58) and Gardner
[\[43\]](#ref43). Exponential smoothing is widely used in supply chain
management for short-term demand forecasting of retail products.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_t = \alpha
\textbf{y}\_{t-1} + (1 - \alpha) \mu\_{t-1}, \quad t=2,\dots,T\\
\\\alpha \sim \mathcal{U}(0,1)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The exponential smoothing model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the smoothed mean \\\mu_t\\ depends recursively on both
\\\textbf{y}\_{t-1}\\ and \\\mu\_{t-1}\\, which requires imperative loop
constructs that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42134)
T <- 250

#### True parameter values
true.alpha <- 0.3
true.sigma <- 0.5

#### Simulate exponential smoothing process
y <- numeric(T)
mu <- numeric(T)
y[1] <- rnorm(1, 0, true.sigma)
mu[1] <- y[1]
for (t in 2:T) {
    mu[t] <- true.alpha * y[t - 1] + (1 - true.alpha) * mu[t - 1]
    y[t] <- rnorm(1, mu[t], true.sigma)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- c("alpha","sigma")
PGF <- function(Data) {
    alpha <- runif(1)
    sigma <- runif(1)
    return(c(alpha, sigma))
    }
Data <- list(PGF=PGF, T=T, mon.names=mon.names, parm.names=parm.names,
    y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[1] <- alpha <- interval(parm[1], 0, 1)
    parm[2] <- sigma <- interval(parm[2], 1e-100, Inf)
    ### Log-Prior
    alpha.prior <- dunif(alpha, 0, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- y
    mu[-1] <- alpha*Data$y[-1]
    mu[-1] <- mu[-1] + (1 - alpha) * mu[-Data$T]
    LL <- sum(dnorm(Data$y[-1], mu[-Data$T], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0.5, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[1, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[2, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("alpha", parm.names)] <- true.alpha[seq_along(grep("alpha", parm.names))]
ground_truth[grep("sigma", parm.names)] <- true.sigma[seq_along(grep("sigma", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Logistic Autoregression **\[NEW\]**

Logistic autoregression extends the binary logit model by including the
lagged binary response \\y\_{t-1}\\ as a predictor alongside exogenous
covariates, creating a first-order autoregressive structure on the logit
scale. The framework was developed by Kedem and Fokianos
[\[152\]](#ref152) as part of a broader treatment of regression models
for time series of counts and binary outcomes. The autoregressive
coefficient \\\phi\\ captures the persistence of the binary state: large
positive values indicate that the outcome tends to repeat, while values
near zero reduce the model to a standard logistic regression. One
exogenous covariate enters through a linear coefficient \\\beta\\, and
the intercept \\\alpha\\ sets the baseline log-odds when both the lagged
response and the covariate are zero. Logistic autoregression is applied
in epidemiology to model the recurrence of disease episodes and in
finance for modelling sequences of trading signals.

#### Form

\\y_t \sim \text{Bernoulli}(p_t), \quad t=2,\dots,T\\
\\\text{logit}(p_t) = \alpha + \phi\\ y\_{t-1} + \beta\\ x_t\\ \\\alpha
\sim \mathcal{N}(0, 100)\\ \\\phi \sim \mathcal{N}(0, 100)\\ \\\beta
\sim \mathcal{N}(0, 100)\\

#### model_spec() notation

The autoregressive structure on lagged binary responses requires
explicit indexing into \\y\_{t-1}\\ within the linear predictor, which
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function handles this directly.

#### Ground truth and data

``` r

set.seed(42220)
TT <- 300

#### True parameter values
true.alpha <- -0.5
true.phi <- 1.2
true.beta <- 0.8

#### Generate exogenous covariate and binary AR series
x <- rnorm(TT)
y <- integer(TT)
y[1] <- 0L
for (t in 2:TT) {
    p.t <- plogis(true.alpha + true.phi * y[t-1] + true.beta * x[t])
    y[t] <- rbinom(1, 1, p.t)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=0, beta=0))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- rnorm(1)
    beta <- rnorm(1)
    return(c(alpha, phi, beta))
    }
Data <- list(TT=TT, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.beta=pos.beta, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    beta <- parm[Data$pos.beta]
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 100, log=TRUE)
    phi.prior <- dnormv(phi, 0, 100, log=TRUE)
    beta.prior <- dnormv(beta, 0, 100, log=TRUE)
    ### Log-Likelihood
    p <- plogis(alpha + phi * Data$y[1:(Data$TT-1)] + beta * Data$x[2:Data$TT])
    LL <- sum(dbinom(Data$y[2:Data$TT], 1, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + beta.prior
    ### yhat: pad to length TT (first observation has no prediction)
    yhat <- c(Data$y[1], rbinom(Data$TT - 1, 1, p))
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0, 0)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("phi   -- true:", true.phi,
    " post. mean:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")
cat("beta  -- true:", true.beta,
    " post. mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### LSTAR

The logistic smooth transition autoregressive (LSTAR) model allows the
autoregressive coefficients to switch smoothly between two regimes as a
function of a lagged value of the series, with the transition governed
by a logistic function parameterized by location \\\theta\\ and scale
\\\gamma\\. This is specified with a transition function that includes
\\\gamma\\ as the shape parameter, \\\textbf{y}\\ as the transition
variable, \\\theta\\ as the location parameter, and \\d\\ as the delay
parameter. The logistic smooth transition autoregressive model was
proposed by Teräsvirta [\[105\]](#ref105) for modelling regime-dependent
dynamics with smooth transitions. LSTAR models are used in
macroeconomics to capture asymmetric business cycle dynamics where
expansions and recessions exhibit different persistence.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad t=1,\dots,T\\
\\\mu_t = \pi_t (\alpha_1 + \phi_1 \textbf{y}\_{t-1}) + (1 - \pi_t)
(\alpha_2 + \phi_2 \textbf{y}\_{t-1}), \quad t=2,\dots,T\\ \\\pi_t =
\frac{1}{1 + \exp(-(\gamma (\textbf{y}\_{t-d} - \theta)))}\\ \\\alpha_j
\sim \mathcal{N}(0, 1000) \in \[\textbf{y}\_{min}, \textbf{y}\_{max}\],
\quad j=1,\dots,2\\ \\\frac{\phi_j+1}{2} \sim \mathcal{BETA}(1, 1),
\quad j=1,\dots,2\\ \\\gamma \sim \mathcal{HC}(25)\\ \\\theta \sim
\mathcal{U}(\textbf{y}\_{min}, \textbf{y}\_{max})\\ \\\pi_1 \sim
\mathcal{U}(0.001, 0.999)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The LSTAR model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean equation involves a logistic transition function
\\\pi_t\\ that depends on lagged values of the series
\\\textbf{y}\_{t-d}\\ and switches between two AR(1) regimes, requiring
imperative constructs that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42137)
T <- 300

#### True parameter values
true.alpha <- c(0.5, -0.3)   # Intercepts for regime 1 and 2
true.phi   <- c(0.8, -0.5)   # AR coefficients for regime 1 and 2
true.gamma <- 5.0             # Transition sharpness
true.theta <- 0.0             # Transition location
true.sigma <- 0.3

#### Simulate LSTAR process
y <- numeric(T)
y[1] <- rnorm(1, true.alpha[1], true.sigma)
for (t in 2:T) {
    pi.t <- 1 / (1 + exp(-(true.gamma * (y[t - 1] - true.theta))))
    mu.t <- pi.t * (true.alpha[1] + true.phi[1] * y[t - 1]) +
        (1 - pi.t) * (true.alpha[2] + true.phi[2] * y[t - 1])
    y[t] <- rnorm(1, mu.t, true.sigma)
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,2), phi=rep(0,2), gamma=0,
    theta=0, pi=0, sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.theta <- grep("theta", parm.names)
pos.pi <- grep("pi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- runif(2,min(Data$y),max(Data$y))
    phi <- runif(2, -1, 1)
    gamma <- runif(1)
    theta <- runif(1,min(Data$y),max(Data$y))
    pi <- runif(1, 0.001, 0.999)
    sigma <- runif(1)
    return(c(alpha, phi, gamma, theta, pi, sigma))
    }
Data <- list(PGF=PGF, T=T, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.phi=pos.phi, pos.gamma=pos.gamma,
    pos.theta=pos.theta, pos.pi=pos.pi, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- interval(parm[Data$pos.alpha], min(Data$y), max(Data$y))
    parm[Data$pos.alpha] <- alpha
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], -1, 1)
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.gamma] <- gamma
    theta <- interval(parm[Data$pos.theta], min(Data$y), max(Data$y))
    parm[Data$pos.theta] <- theta
    parm[Data$pos.pi] <- pi <- interval(parm[Data$pos.pi], 0.001, 0.999)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    phi.prior <- sum(dbeta((phi+1)/2, 1, 1, log=TRUE))
    gamma.prior <- dhalfcauchy(gamma, 25, log=TRUE)
    theta.prior <- dunif(theta, min(Data$y), max(Data$y), log=TRUE)
    pi.prior <- dunif(pi, 0.001, 0.999, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    pi <- c(pi, 1 / (1 + exp(-(gamma*(Data$y[-Data$T]-theta)))))
    mu <- pi * c(alpha[1], alpha[1] + phi[1]*Data$y[-Data$T]) +
        (1-pi) * c(alpha[2], alpha[2] + phi[2]*Data$y[-Data$T])
    LL <- sum(dnorm(Data$y[-1], mu[-1], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + gamma.prior + theta.prior +
        pi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(mean(y),2), rep(0.5,2), 1, mean(y), 0.5, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha[1] -- true:", true.alpha[1],
    " post. mean:", round(fit$Summary2[pos.alpha[1], "Mean"], 3), "\n")
cat("alpha[2] -- true:", true.alpha[2],
    " post. mean:", round(fit$Summary2[pos.alpha[2], "Mean"], 3), "\n")
cat("phi[1] -- true:", true.phi[1],
    " post. mean:", round(fit$Summary2[pos.phi[1], "Mean"], 3), "\n")
cat("phi[2] -- true:", true.phi[2],
    " post. mean:", round(fit$Summary2[pos.phi[2], "Mean"], 3), "\n")
cat("gamma -- true:", true.gamma,
    " post. mean:", round(fit$Summary2[pos.gamma, "Mean"], 3), "\n")
cat("theta -- true:", true.theta,
    " post. mean:", round(fit$Summary2[pos.theta, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Panel Mixed-Effects AR(p) **\[NEW\]**

The Bayesian spatial mixed-effects panel autoregressive model extends
the univariate AR(p) to panel data observed across \\m\\ spatial
locations, where each location has its own AR coefficient vector
\\\theta_k\\ drawn from a common normal distribution centred on a global
mean \\\mu\\, and location-specific residual standard deviations
\\\sigma_k\\ are linked through a log-normal hierarchy. A shared
temporal random effect \\\eta_t\\ captures year-specific shocks common
to all locations, inducing cross-sectional correlation without requiring
an explicit spatial weight matrix. The model was developed by Samia et
al. [\[172\]](#ref172) to analyse multiannual population cycles of
grey-sided voles (*Myodes rufocanus*) across monitoring stations in
Hokkaido, Japan, building on the panel time series framework of Hjellvik
and Tjostheim [\[173\]](#ref173). Mixed-effects panel AR models are
applied in ecology to detect density-dependent regulation and cyclic
dynamics across spatially replicated populations. A companion example
using the RJ algorithm for lag-order selection appears in [Panel
Mixed-Effects AR, RJ](#panel.me.ar.rj).

#### Form

\\Y\_{k,t} \sim \mathcal{N}(\mu\_{k,t},\\ \sigma^2_k), \quad
k=1,\dots,m, \quad t=p+1,\dots,T\\ \\\mu\_{k,t} = \theta\_{k,0} +
\sum\_{j=1}^{p} \theta\_{k,j}\\ Y\_{k,t-j} + \eta\_{t-p}\\ \\\theta_k
\sim \mathcal{N}\_{p+1}(\mu,\\ \xi^2 \textbf{I}), \quad k=1,\dots,m\\
\\\sigma_k \sim \text{LogNormal}(\kappa, \psi), \quad k=1,\dots,m\\
\\\eta_t \sim \mathcal{N}(0, \omega^2), \quad t=1,\dots,n\\ \\\mu_j \sim
\mathcal{N}(0, 1000), \quad j=0,\dots,p\\ \\\xi \sim \mathcal{HC}(25)\\
\\\omega \sim \mathcal{HC}(25)\\ \\\psi \sim \mathcal{HC}(25)\\ \\\kappa
\sim \mathcal{N}(0, 10)\\

where \\n = T - p\\ is the number of usable time points after discarding
the first \\p\\ lags.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it requires panel-indexed autoregressive dynamics where each
location \\k\\ has its own AR coefficient vector \\\theta_k\\ applied to
location-specific lagged observations \\Y\_{k,t-j}\\, combined with a
shared temporal random effect vector \\\eta_t\\ that enters the mean of
every location simultaneously. This structure of crossed location and
time indexing with hierarchical priors on location-specific parameters
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42233)
m <- 5       # Number of locations
TT <- 30     # Time points per location
p <- 2       # AR order
n <- TT - p  # Usable time points after discarding p lags

#### True parameter values
true.mu <- c(0.3, 0.5, -0.3)    # global mean: intercept, AR(1), AR(2)
true.xi <- 0.15                   # hierarchical SD for theta
true.omega <- 0.20                # SD of temporal random effects
true.kappa <- log(0.4)            # mean of log(sigma.k)
true.psi <- 0.30                  # SD of log(sigma.k) across locations

#### Location-specific AR coefficients from the hierarchical prior
true.theta <- matrix(NA, m, p + 1)
for (k in 1:m)
    true.theta[k, ] <- rnorm(p + 1, true.mu, true.xi)

#### Location-specific residual SDs from the log-normal hierarchy
true.sigma.k <- exp(rnorm(m, true.kappa, true.psi))

#### Common temporal random effects
true.eta <- rnorm(n, 0, true.omega)

#### Simulate panel AR(p) data: Y is m x TT
Y <- matrix(NA, m, TT)
for (k in 1:m) {
    Y[k, 1:p] <- rnorm(p, true.theta[k, 1] / (1 - sum(true.theta[k, -1])),
        true.sigma.k[k])
}
for (t in (p + 1):TT) {
    for (k in 1:m) {
        mu.kt <- true.theta[k, 1]
        for (j in 1:p)
            mu.kt <- mu.kt + true.theta[k, j + 1] * Y[k, t - j]
        mu.kt <- mu.kt + true.eta[t - p]
        Y[k, t] <- rnorm(1, mu.kt, true.sigma.k[k])
    }
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=rep(0, p + 1),
    theta=matrix(0, m, p + 1), sigma.k=rep(0, m), eta=rep(0, n),
    xi=0, omega=0, psi=0, kappa=0))
pos.mu <- grep("^mu\\[", parm.names)
pos.theta <- grep("theta", parm.names)
pos.sigma.k <- grep("sigma.k", parm.names)
pos.eta <- grep("^eta\\[", parm.names)
pos.xi <- which(parm.names == "xi")
pos.omega <- which(parm.names == "omega")
pos.psi <- which(parm.names == "psi")
pos.kappa <- which(parm.names == "kappa")
PGF <- function(Data) {
    mu <- rnorm(Data$p + 1)
    xi <- runif(1, 0.05, 0.5)
    theta <- rnorm(Data$m * (Data$p + 1), rep(mu, each=Data$m), xi)
    sigma.k <- exp(rnorm(Data$m, 0, 0.5))
    eta <- rnorm(Data$n, 0, 0.2)
    omega <- runif(1, 0.05, 0.5)
    psi <- runif(1, 0.05, 0.5)
    kappa <- rnorm(1, 0, 1)
    return(c(mu, theta, sigma.k, eta, xi, omega, psi, kappa))
    }
Data <- list(m=m, n=n, p=p, PGF=PGF, TT=TT, Yall=Y,
    Y=Y[, (p + 1):TT], mon.names=mon.names, parm.names=parm.names,
    pos.mu=pos.mu, pos.theta=pos.theta, pos.sigma.k=pos.sigma.k,
    pos.eta=pos.eta, pos.xi=pos.xi, pos.omega=pos.omega,
    pos.psi=pos.psi, pos.kappa=pos.kappa)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    xi <- interval(parm[Data$pos.xi], 1e-100, Inf)
    parm[Data$pos.xi] <- xi
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    psi <- interval(parm[Data$pos.psi], 1e-100, Inf)
    parm[Data$pos.psi] <- psi
    kappa <- parm[Data$pos.kappa]
    ### Parameters
    mu <- parm[Data$pos.mu]
    theta <- matrix(parm[Data$pos.theta], Data$m, Data$p + 1)
    sigma.k <- interval(parm[Data$pos.sigma.k], 1e-100, Inf)
    parm[Data$pos.sigma.k] <- sigma.k
    eta <- parm[Data$pos.eta]
    ### Log-Hyperprior
    xi.prior <- dhalfcauchy(xi, 25, log=TRUE)
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    psi.prior <- dhalfcauchy(psi, 25, log=TRUE)
    kappa.prior <- dnormv(kappa, 0, 10, log=TRUE)
    ### Log-Prior
    mu.prior <- sum(dnormv(mu, 0, 1000, log=TRUE))
    theta.prior <- sum(dnorm(theta, rep(mu, each=Data$m), xi, log=TRUE))
    sigma.k.prior <- sum(dlnorm(sigma.k, kappa, psi, log=TRUE))
    eta.prior <- sum(dnorm(eta, 0, omega, log=TRUE))
    ### Log-Likelihood
    mu.mat <- matrix(theta[, 1], Data$m, Data$n)
    for (j in 1:Data$p)
        mu.mat <- mu.mat + theta[, j + 1] *
            Data$Yall[, (Data$p + 1 - j):(Data$TT - j)]
    mu.mat <- mu.mat + matrix(eta, Data$m, Data$n, byrow=TRUE)
    LL <- sum(dnorm(Data$Y, mu.mat, sigma.k, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + theta.prior + sigma.k.prior + eta.prior +
        xi.prior + omega.prior + psi.prior + kappa.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$m * Data$n, mu.mat, sigma.k), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, p + 1), rep(0, m * (p + 1)), rep(0.5, m),
    rep(0, n), 0.2, 0.2, 0.2, 0)
```

#### Fitting and recovery

``` r

#### fit using AIES (ensemble MCMC, robust to multimodality)
fit <- lucifer(Model, Data, Initial.Values,
    Covar = NULL, Iterations = 100000, Status = 1000, Thinning = 100,
    Algorithm = "AIES",
    Specs = list(Nc = 2 * length(Initial.Values), Z = NULL,
        beta = 2, CPUs = 1, Packages = NULL, Dyn.libs = NULL))
print(fit)
Consort(fit)

#### Parameter recovery: global AR means
S <- if (all(is.na(fit$Summary2[, "Mean"]))) fit$Summary1 else fit$Summary2
for (j in 1:(p + 1))
    cat("mu[", j, "] -- true:", round(true.mu[j], 3),
        " post. mean:", round(S[pos.mu[j], "Mean"], 3), "\n")

#### Parameter recovery: hierarchical SDs and hyperparameters
cat("xi -- true:", true.xi,
    " post. mean:", round(S[pos.xi, "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(S[pos.omega, "Mean"], 3), "\n")
cat("psi -- true:", true.psi,
    " post. mean:", round(S[pos.psi, "Mean"], 3), "\n")
cat("kappa -- true:", round(true.kappa, 3),
    " post. mean:", round(S[pos.kappa, "Mean"], 3), "\n")

#### Parameter recovery: location-specific AR coefficients (first location)
cat("\nLocation 1 AR coefficients:\n")
theta.idx <- matrix(pos.theta, m, p + 1)
for (j in 1:(p + 1))
    cat("  theta[1,", j, "] -- true:", round(true.theta[1, j], 3),
        " post. mean:", round(S[theta.idx[1, j], "Mean"], 3), "\n")

#### Parameter recovery: location-specific SDs
cat("\nLocation-specific residual SDs:\n")
for (k in 1:m)
    cat("  sigma.k[", k, "] -- true:", round(true.sigma.k[k], 3),
        " post. mean:", round(S[pos.sigma.k[k], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.eta]   <- true.eta[seq_along(pos.eta)]
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
ground_truth[grep("kappa", parm.names)] <- true.kappa[seq_along(grep("kappa", parm.names))]
ground_truth[grep("omega", parm.names)] <- true.omega[seq_along(grep("omega", parm.names))]
ground_truth[grep("psi", parm.names)]   <- true.psi[seq_along(grep("psi", parm.names))]
ground_truth[grep("xi", parm.names)]    <- true.xi[seq_along(grep("xi", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Panel Mixed-Effects AR, RJ **\[NEW\]**

This example extends the fixed-order [Panel Mixed-Effects
AR(p)](#panel.me.arp) by using the RJ (Reversible-Jump) algorithm for
automatic lag-order selection. Instead of fixing \\p\\, the model is
specified at the maximum order \\p\_{\max} = 5\\, and the RJ mechanism
proposes birth and death moves that include or exclude individual global
lag means \\\mu_j\\ from the active set, effectively determining which
lags contribute to the dynamics. When a global mean \\\mu_j\\ is set to
zero by the RJ mechanism, the hierarchical prior \\\theta\_{k,j} \sim
\mathcal{N}(\mu_j, \xi^2)\\ centres the location-specific coefficients
at zero, providing soft shrinkage rather than hard exclusion. This
approach differs from the paper’s global order selection: while Samia et
al. [\[172\]](#ref172) use dimension-jumping RJMCMC to select a single
\\p\\ shared across all locations in a subgroup, lucifer’s RJ selects on
the global lag means and lets the hierarchical prior propagate the
selection to the location-specific coefficients. Posterior inclusion
frequencies for each \\\mu_j\\ reveal the effective AR order, and lags
excluded by the RJ mechanism have posterior means near zero at both
global and location levels. The recommended initialisation strategy is
to start with all lag means selected (active) and let the RJ death moves
prune those that do not contribute to the likelihood; this works well
because death proposals deterministically set a coefficient to zero,
which the Metropolis-Hastings step can evaluate exactly, whereas birth
proposals perturb from zero by a small random amount and may require
many attempts to reach the true coefficient magnitude.

#### Form

\\Y\_{k,t} \sim \mathcal{N}(\mu\_{k,t},\\ \sigma^2_k), \quad
k=1,\dots,m, \quad t=p\_{\max}+1,\dots,T\\ \\\mu\_{k,t} =
\theta\_{k,0} + \sum\_{j=1}^{p\_{\max}} \theta\_{k,j}\\ Y\_{k,t-j} +
\eta\_{t-p\_{\max}}\\ \\\theta_k \sim \mathcal{N}\_{p\_{\max}+1}(\mu,\\
\xi^2 \textbf{I}), \quad k=1,\dots,m\\ \\\sigma_k \sim
\text{LogNormal}(\kappa, \psi), \quad k=1,\dots,m\\ \\\eta_t \sim
\mathcal{N}(0, \omega^2), \quad t=1,\dots,n\\ \\\mu_j \sim
\mathcal{N}(0, 1000), \quad j=0,\dots,p\_{\max}\\ \\\xi \sim
\mathcal{HC}(25), \quad \omega \sim \mathcal{HC}(25), \quad \psi \sim
\mathcal{HC}(25)\\ \\\kappa \sim \mathcal{N}(0, 10)\\

The RJ algorithm selects which global lag means \\\mu_j\\ for \\j \geq
1\\ are active; inactive means are set to zero, and the hierarchical
prior \\\theta\_{k,j} \sim \mathcal{N}(\mu_j, \xi^2)\\ shrinks the
corresponding location-specific coefficients toward zero. Intercepts
(\\\mu_0\\), all location-specific coefficients \\\theta\_{k,j}\\,
location SDs, temporal effects, and hyperparameters are always active
(not selectable). The chain is initialised with all five lag means
selected, allowing the RJ death moves to remove superfluous lags during
burn-in.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it requires reversible-jump specifications (birth/death move
probabilities, model-size prior, and a selectable indicator vector) that
operate on the discrete model-indicator space, in addition to the
panel-indexed autoregressive dynamics and hierarchical location-specific
parameters. The manual `Model` function below handles this directly.

#### Ground truth and data

The data are generated from an AR(2) process (true lag coefficients for
orders 3, 4, and 5 are exactly zero), so the RJ algorithm should
preferentially include lags 1 and 2 while excluding lags 3 through 5.

``` r

set.seed(42234)
m <- 10       # Number of locations
TT <- 50      # Time points per location
p.max <- 5    # Maximum AR order
p.true <- 2   # True AR order
n <- TT - p.max  # Usable time points

#### True parameter values (AR(2) embedded in p.max=5 space)
true.mu <- c(0.3, 0.5, -0.3, 0, 0, 0)  # intercept + 5 lags; lags 3-5 zero
true.xi <- 0.15
true.omega <- 0.20
true.kappa <- log(0.4)
true.psi <- 0.30

#### Location-specific AR coefficients
true.theta <- matrix(NA, m, p.max + 1)
for (k in 1:m) {
    true.theta[k, 1:(p.true + 1)] <- rnorm(p.true + 1,
        true.mu[1:(p.true + 1)], true.xi)
    true.theta[k, (p.true + 2):(p.max + 1)] <- 0  # inactive lags
}

#### Location-specific residual SDs
true.sigma.k <- exp(rnorm(m, true.kappa, true.psi))

#### Common temporal random effects
true.eta <- rnorm(n, 0, true.omega)

#### Simulate panel data (using all p.max lags for indexing consistency)
Y <- matrix(NA, m, TT)
for (k in 1:m)
    Y[k, 1:p.max] <- rnorm(p.max, true.theta[k, 1], true.sigma.k[k])
for (t in (p.max + 1):TT) {
    for (k in 1:m) {
        mu.kt <- true.theta[k, 1]
        for (j in 1:p.true)
            mu.kt <- mu.kt + true.theta[k, j + 1] * Y[k, t - j]
        mu.kt <- mu.kt + true.eta[t - p.max]
        Y[k, t] <- rnorm(1, mu.kt, true.sigma.k[k])
    }
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=rep(0, p.max + 1),
    theta=matrix(0, m, p.max + 1), sigma.k=rep(0, m), eta=rep(0, n),
    xi=0, omega=0, psi=0, kappa=0))
pos.mu <- grep("^mu\\[", parm.names)
pos.theta <- grep("theta", parm.names)
pos.sigma.k <- grep("sigma.k", parm.names)
pos.eta <- grep("^eta\\[", parm.names)
pos.xi <- which(parm.names == "xi")
pos.omega <- which(parm.names == "omega")
pos.psi <- which(parm.names == "psi")
pos.kappa <- which(parm.names == "kappa")
PGF <- function(Data) {
    mu <- rnorm(Data$p.max + 1)
    xi <- runif(1, 0.05, 0.5)
    theta <- rnorm(Data$m * (Data$p.max + 1), rep(mu, each=Data$m), xi)
    sigma.k <- exp(rnorm(Data$m, 0, 0.5))
    eta <- rnorm(Data$n, 0, 0.2)
    omega <- runif(1, 0.05, 0.5)
    psi <- runif(1, 0.05, 0.5)
    kappa <- rnorm(1, 0, 1)
    return(c(mu, theta, sigma.k, eta, xi, omega, psi, kappa))
    }
Data <- list(m=m, n=n, p.max=p.max, PGF=PGF, TT=TT, Yall=Y,
    Y=Y[, (p.max + 1):TT], mon.names=mon.names, parm.names=parm.names,
    pos.mu=pos.mu, pos.theta=pos.theta, pos.sigma.k=pos.sigma.k,
    pos.eta=pos.eta, pos.xi=pos.xi, pos.omega=pos.omega,
    pos.psi=pos.psi, pos.kappa=pos.kappa)

#### Reversible-Jump specifications
LIV <- length(parm.names)
selectable <- rep(0, LIV)
## Selectable: only global lag means mu[2] through mu[p.max+1]
## The hierarchical prior shrinks theta[k,j] toward mu[j]; when mu[j]=0,
## the location-specific coefficients cluster near zero automatically
selectable[pos.mu[2:(p.max + 1)]] <- 1
## Start with ALL lags selected; RJ death moves prune inactive lags
selected <- rep(0, LIV)
selected[pos.mu[2:(p.max + 1)]] <- 1
bin.n <- sum(selectable)  # max selectable = p.max = 5
bin.p <- 0.5              # symmetric model size prior
parm.p <- rep(0.5, LIV)  # flat inclusion prior (no Bernoulli penalty)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    xi <- interval(parm[Data$pos.xi], 1e-100, Inf)
    parm[Data$pos.xi] <- xi
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    psi <- interval(parm[Data$pos.psi], 1e-100, Inf)
    parm[Data$pos.psi] <- psi
    kappa <- parm[Data$pos.kappa]
    ### Parameters
    mu <- parm[Data$pos.mu]
    theta <- matrix(parm[Data$pos.theta], Data$m, Data$p.max + 1)
    sigma.k <- interval(parm[Data$pos.sigma.k], 1e-100, Inf)
    parm[Data$pos.sigma.k] <- sigma.k
    eta <- parm[Data$pos.eta]
    ### Log-Hyperprior
    xi.prior <- dhalfcauchy(xi, 25, log=TRUE)
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    psi.prior <- dhalfcauchy(psi, 25, log=TRUE)
    kappa.prior <- dnormv(kappa, 0, 10, log=TRUE)
    ### Log-Prior
    mu.prior <- sum(dnormv(mu, 0, 1000, log=TRUE))
    theta.prior <- sum(dnorm(theta, rep(mu, each=Data$m), xi, log=TRUE))
    sigma.k.prior <- sum(dlnorm(sigma.k, kappa, psi, log=TRUE))
    eta.prior <- sum(dnorm(eta, 0, omega, log=TRUE))
    ### Log-Likelihood
    mu.mat <- matrix(theta[, 1], Data$m, Data$n)
    for (j in 1:Data$p.max)
        mu.mat <- mu.mat + theta[, j + 1] *
            Data$Yall[, (Data$p.max + 1 - j):(Data$TT - j)]
    mu.mat <- mu.mat + matrix(eta, Data$m, Data$n, byrow=TRUE)
    LL <- sum(dnorm(Data$Y, mu.mat, sigma.k, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + theta.prior + sigma.k.prior + eta.prior +
        xi.prior + omega.prior + psi.prior + kappa.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$m * Data$n, mu.mat, sigma.k), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0.1, rep(0.1, p.max), rep(0, m * (p.max + 1)),
    rep(0.5, m), rep(0, n), 0.2, 0.2, 0.2, 0)
```

#### Fitting and recovery

``` r

#### fit using the RJ algorithm
fit <- lucifer(Model, Data, Initial.Values,
    Covar = NULL, Iterations = 20000, Status = 2000, Thinning = 20,
    Algorithm = "RJ",
    Specs = list(bin.n = bin.n, bin.p = bin.p,
                 parm.p = parm.p, selectable = selectable,
                 selected = selected))
print(fit)
summary(fit)

#### Lag-order selection: identify active lags from posterior means
## RJ sets excluded parameters to exactly zero; nonzero posterior mean
## indicates the lag was included in some proportion of posterior samples.
## Lags 1-2 (true AR order) should retain substantial posterior means,
## while lags 3-5 (inactive) should have means near zero.
S <- if (all(is.na(fit$Summary2[, "Mean"]))) fit$Summary1 else fit$Summary2
post.mu <- S[pos.mu, "Mean"]
cat("Lag-order selection (global means):\n")
for (j in 1:p.max)
    cat("  Lag", j, "-- mu post. mean:", round(post.mu[j + 1], 4),
        " (true:", round(true.mu[j + 1], 3),
        ifelse(j <= p.true, ", active", ", inactive"), ")\n")

#### Parameter recovery: global AR means
cat("\nGlobal AR means:\n")
for (j in 1:(p.max + 1))
    cat("  mu[", j, "] -- true:", round(true.mu[j], 3),
        " post. mean:", round(post.mu[j], 3), "\n")

#### Parameter recovery: hierarchical SDs
cat("\nxi -- true:", true.xi,
    " post. mean:", round(S[pos.xi, "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(S[pos.omega, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.eta]   <- true.eta[seq_along(pos.eta)]
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
ground_truth[grep("kappa", parm.names)] <- true.kappa[seq_along(grep("kappa", parm.names))]
ground_truth[grep("omega", parm.names)] <- true.omega[seq_along(grep("omega", parm.names))]
ground_truth[grep("psi", parm.names)]   <- true.psi[seq_along(grep("psi", parm.names))]
ground_truth[grep("xi", parm.names)]    <- true.xi[seq_along(grep("xi", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Threshold Autoregression (TAR)

The threshold autoregression model allows abrupt regime switching in the
autoregressive coefficients at an unknown change point \\\theta\\;
observations before and after the threshold follow distinct AR(1)
processes. The threshold autoregressive model was proposed by Tong
[\[107\]](#ref107) and Tong and Lim [\[108\]](#ref108) for modelling
regime-dependent dynamics. TAR models are applied in ecology to detect
threshold responses of lake ecosystems to nutrient loading, where the
system switches between oligotrophic and eutrophic regimes.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\nu_t, \sigma^2), \quad t=1,\dots,T\\
\\\nu_t = \left\\ \begin{array}{l l} \alpha_1 + \phi_1
\textbf{y}\_{t-1}, \quad t=1,\dots,T & \quad \mbox{if \$t \ge
\theta\$}\\ \alpha_2 + \phi_2 \textbf{y}\_{t-1}, \quad t=1,\dots,T &
\quad \mbox{if \$t \< \theta\$} \\ \end{array} \right. \\ \\\alpha_j
\sim \mathcal{N}(0, 1000) \in \[-1,1\], \quad j=1,\dots,2\\ \\\phi_j
\sim \mathcal{N}(0, 1000) \in \[-1,1\], \quad j=1,\dots,2\\ \\\theta
\sim \mathcal{U}(2, T-1)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The TAR model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean equation involves a regime indicator that depends
on a latent change point \\\theta\\ applied to the time index, and the
autoregressive structure requires recursive time indexing. These
imperative constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42157)
T <- 260

#### True parameter values
true.alpha <- c(0.1, -0.1)    # Intercepts for regime 1 (t < theta) and regime 2
true.phi   <- c(0.7, -0.4)    # AR(1) coefficients for each regime
true.theta <- 130              # Change point
true.sigma <- 0.3

#### Simulate TAR process
y <- numeric(T)
y[1] <- rnorm(1, true.alpha[1], true.sigma)
for (t in 2:T) {
    if (t < true.theta) {
        y[t] <- rnorm(1, true.alpha[1] + true.phi[1] * y[t - 1], true.sigma)
    } else {
        y[t] <- rnorm(1, true.alpha[2] + true.phi[2] * y[t - 1], true.sigma)
    }
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,2), phi=rep(0,2), theta=0,
    sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.theta <- grep("theta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rtrunc(2, "norm", a=-1, b=1, mean=0, sd=1)
    phi <- rtrunc(2, "norm", a=-1, b=1, mean=0, sd=1)
    theta <- runif(1,2,Data$T-1)
    sigma <- runif(1)
    return(c(alpha, phi, theta, sigma))
    }
Data <- list(PGF=PGF, T=T, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.phi=pos.phi, pos.theta=pos.theta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    parm[Data$pos.alpha] <- alpha <- interval(parm[Data$pos.alpha], -1, 1)
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], -1, 1)
    theta <- interval(parm[Data$pos.theta], 2, Data$T-1)
    parm[Data$pos.theta] <- theta
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dtrunc(alpha, "norm", a=-1, b=1, mean=0,
        sd=sqrt(1000), log=TRUE))
    phi.prior <- sum(dtrunc(phi, "norm", a=-1, b=1, mean=0,
        sd=sqrt(1000), log=TRUE))
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    theta.prior <- dunif(theta, 2, Data$T-1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- matrix(0, Data$T, 2)
    mu[,1] <- c(alpha[1], alpha[1] + phi[1]*Data$y[-Data$T])
    mu[,2] <- c(alpha[2], alpha[2] + phi[2]*Data$y[-Data$T])
    nu <- mu[,2]; temp <- which(1:Data$T < theta)
    nu[temp] <- mu[temp,1]
    LL <- sum(dnorm(Data$y[-1], nu[-1], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + theta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(nu), nu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,4), T/2, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha[1] -- true:", true.alpha[1],
    " post. mean:", round(fit$Summary2[pos.alpha[1], "Mean"], 3), "\n")
cat("alpha[2] -- true:", true.alpha[2],
    " post. mean:", round(fit$Summary2[pos.alpha[2], "Mean"], 3), "\n")
cat("phi[1] -- true:", true.phi[1],
    " post. mean:", round(fit$Summary2[pos.phi[1], "Mean"], 3), "\n")
cat("phi[2] -- true:", true.phi[2],
    " post. mean:", round(fit$Summary2[pos.phi[2], "Mean"], 3), "\n")
cat("theta -- true:", true.theta,
    " post. mean:", round(fit$Summary2[pos.theta, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Time Varying AR(1) with Chebyshev Series

This example consists of a first-order autoregressive model, AR(1), with
a time-varying parameter (TVP) \\\phi\\, that is a Chebyshev series
constructed from a linear combination of orthonormal Chebyshev time
polynomials (CTPs) and parameter vector \\\beta\\. The user creates
basis matrix **P**, specifying polynomial degree \\D\\ and time \\T\\.
Each column is a CTP of a different degree, and the first column is
restricted to 1, the linear basis. CTPs are very flexible for TVPs, and
estimate quickly because each is orthogonal, unlike simple polynomials
and splines. Time-varying autoregressive models with Chebyshev
polynomial coefficients follow the functional coefficient approach of
Chen and Tsay [\[22\]](#ref22). Chebyshev-based time-varying AR models
are used in seismology to capture the non-stationary spectral properties
of earthquake ground motion records.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \phi\_{t-1} \textbf{y}\_{t-1}\\ \\\phi_t = \textbf{P}
\beta\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\beta_d \sim
\mathcal{N}(0, 1000), \quad d=1,\dots,(D+1)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### Data

``` r

data(demonfx)
y <- as.vector(diff(log(as.matrix(demonfx[1:261,1]))))
D <- 6 #Maximum degree of Chebyshev time polynomials
T <- length(y)
P <- matrix(1, T, D+1)
for (d in 1:D) {P[,d+1] <- sqrt(2)*cos(d*pi*(c(1:T)-0.5)/T)}
mon.names <- c("LP", "ynew", as.parm.names(list(phi=rep(0,T-1))))
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,D+1), sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(Data$D+1)
    sigma <- runif(1)
    return(c(alpha, beta, sigma))
    }
Data <- list(D=D, P=P, PGF=PGF, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    phi <- tcrossprod(Data$P[-Data$T,], t(beta))
    mu <- c(alpha, alpha + phi*Data$y[-Data$T])
    ynew <- rnorm(1, alpha + tcrossprod(Data$P[Data$T,], t(beta))*
        Data$y[Data$T], sigma)
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew,phi),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,D+2), 1)
```

#### model_spec() notation

The TVAR-Chebyshev model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the time-varying autoregressive coefficient \\\phi_t =
\textbf{P}\beta\\ is constructed from a user-supplied Chebyshev basis
matrix **P** and applied element-wise to lagged observations, requiring
imperative matrix algebra that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Inspect the time-varying AR(1) coefficient
phi.post <- fit$Monitor[, grep("phi", mon.names)]
phi.mean <- colMeans(phi.post)
plot(1:(T-1), phi.mean, type = "l", xlab = "Time", ylab = expression(phi[t]),
    main = "Posterior mean of time-varying AR(1) coefficient")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Time series: volatility models

### AR(p)-ARCH(q)

The AR(p)-ARCH(q) model combines an autoregressive mean equation with an
autoregressive conditional heteroscedasticity (ARCH) specification for
the variance, allowing the conditional variance to depend on past
squared residuals. It captures both serial dependence in the mean and
volatility clustering in financial and economic time series. The ARCH
model was introduced by Engle [\[35\]](#ref35), who received the Nobel
Prize in Economics for this contribution to modelling time-varying
volatility. AR-ARCH models are widely used in finance to capture the
clustering of volatile periods in daily stock returns.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p}, \quad
t=1,\dots,T\\ \\\epsilon_t = \textbf{y}\_t - \mu_t\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\phi_p \sim \mathcal{N}(0, 1000), \quad
p=1,\dots,P\\ \\\sigma^2_t = \omega + \sum^Q\_{q=1} \theta_q
\epsilon^2\_{t-q}, \quad t=2,\dots,T\\ \\\omega \sim \mathcal{HC}(25)\\
\\\theta_q \sim \mathcal{U}(0, 1), \quad q=1,\dots,Q\\

#### model_spec() notation

The AR(p)-ARCH(q) model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because both the mean and variance equations involve recursive time
indexing (\\\mu_t\\ depends on \\\textbf{y}\_{t-p}\\ and \\\sigma^2_t\\
depends on \\\epsilon^2\_{t-q}\\), which requires imperative loop
constructs that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42107)
T <- 260
L.P <- c(1, 5)        # Autoregressive lags
L.Q <- c(1, 2)        # Volatility lags
P <- length(L.P)      # Autoregressive order
Q <- length(L.Q)      # Volatility order

#### True parameter values
true.alpha <- 0.001
true.phi   <- c(0.15, -0.08)
true.omega <- 0.0001
true.theta <- c(0.25, 0.15)

#### Simulate AR(p)-ARCH(q) process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
y[1:max(L.P)] <- rnorm(max(L.P), true.alpha, sqrt(true.omega))
eps[1:max(L.P)] <- y[1:max(L.P)] - true.alpha
sigma2[1:max(L.P)] <- true.omega
for (t in (max(L.P) + 1):T) {
    sigma2[t] <- true.omega
    for (q in 1:Q)
        if (t > L.Q[q]) sigma2[t] <- sigma2[t] + true.theta[q] * eps[t - L.Q[q]]^2
    mu.t <- true.alpha
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L.P[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), omega=0,
    theta=rep(0,Q)))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(Data$Q, 1e-10, 1-1e-5)
    return(c(alpha, phi, omega, theta))
    }
Data <- list(L.P=L.P, L.Q=L.Q, PGF=PGF, P=P, Q=Q, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- interval(parm[Data$pos.theta], 1e-10, 1-1e-5)
    parm[Data$pos.theta] <- theta
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 1e-10, 1-1e-5, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L.P[p])] <- mu[-c(1:Data$L.P[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L.P[p])]
    epsilon <- Data$y - mu
    sigma2 <- rep(omega, Data$T)
    for (q in 1:Data$Q)
        sigma2[-c(1:Data$L.Q[q])] <- sigma2[-c(1:Data$L.Q[q])] +
    theta[q]*epsilon[1:(Data$T-Data$L.Q[q])]^ 2
    LL <- sum(dnormv(Data$y[-c(1:Data$L.P[Data$P])],
        mu[-c(1:Data$L.P[Data$P])], sigma2[-c(1:Data$L.P[Data$P])],
        log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + omega.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,P+1), 1, rep(0.5,Q))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
for (q in 1:Q)
    cat("theta[", q, "] -- true:", true.theta[q],
        " post. mean:", round(fit$Summary2[pos.theta[q], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.omega] <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### AR(p)-ARCH(q)-M

The ARCH-in-mean extension of AR(p)-ARCH(q) includes the conditional
variance as an additional regressor in the mean equation through the
parameter \\\delta\\, capturing the risk-return tradeoff commonly
hypothesized in financial time series. The ARCH-in-mean specification
was proposed by Engle et al. [\[37\]](#ref37) to model risk premia in
asset returns. In financial economics, ARCH-M models quantify how the
conditional volatility of an asset directly affects its expected return,
testing the risk-return tradeoff.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p} + \delta
\sigma^2\_{t-1}, \quad t=1,\dots,T\\ \\\epsilon_t = \textbf{y}\_t -
\mu_t\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi_p \sim
\mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\delta \sim \mathcal{N}(0,
1000)\\ \\\sigma^2_t = \omega + \sum^Q\_{q=1} \theta_q
\epsilon^2\_{t-q}, \quad t=2,\dots,T\\ \\\omega \sim \mathcal{HC}(25)\\
\\\theta_q \sim \mathcal{U}(0, 1), \quad q=1,\dots,Q\\

#### model_spec() notation

The AR(p)-ARCH(q)-M model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because both the mean and variance equations involve recursive time
indexing (\\\mu_t\\ depends on \\\textbf{y}\_{t-p}\\ and \\\sigma^2_t\\
depends on \\\epsilon^2\_{t-q}\\), and the in-mean term feeds the
conditional variance back into the mean equation. These imperative loop
constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42108)
T <- 260
L.P <- c(1, 5)        # Autoregressive lags
L.Q <- c(1, 2)        # Volatility lags
P <- length(L.P)      # Autoregressive order
Q <- length(L.Q)      # Volatility order

#### True parameter values
true.alpha <- 0.002
true.phi   <- c(0.12, -0.06)
true.delta <- 0.5
true.omega <- 0.0001
true.theta <- c(0.20, 0.15)

#### Simulate AR(p)-ARCH(q)-M process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
y[1:max(L.P)] <- rnorm(max(L.P), true.alpha, sqrt(true.omega))
eps[1:max(L.P)] <- y[1:max(L.P)] - true.alpha
sigma2[1:max(L.P)] <- true.omega
for (t in (max(L.P) + 1):T) {
    sigma2[t] <- true.omega
    for (q in 1:Q)
        if (t > L.Q[q]) sigma2[t] <- sigma2[t] + true.theta[q] * eps[t - L.Q[q]]^2
    mu.t <- true.alpha + true.delta * sigma2[t]
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L.P[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), delta=0, omega=0,
    theta=rep(0,Q)))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.delta <- grep("delta", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    delta <- rnorm(1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(Data$Q, 1e-10, 1-1e-5)
    return(c(alpha, phi, delta, omega, theta))
    }
Data <- list(L.P=L.P, L.Q=L.Q, PGF=PGF, P=P, Q=Q, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.delta=pos.delta, pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    delta <- parm[Data$pos.delta]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- interval(parm[Data$pos.theta], 1e-10, 1-1e-5)
    parm[Data$pos.theta] <- theta
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    delta.prior <- dnormv(delta, 0, 1000, log=TRUE)
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 1e-10, 1-1e-5, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L.P[p])] <- mu[-c(1:Data$L.P[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L.P[p])]
    epsilon <- Data$y - mu
    sigma2 <- rep(omega, Data$T)
    for (q in 1:Data$Q)
        sigma2[-c(1:Data$L.Q[q])] <- sigma2[-c(1:Data$L.Q[q])] +
    theta[q]*epsilon[1:(Data$T-Data$L.Q[q])]^ 2
    mu <- mu + delta*sigma2
    LL <- sum(dnormv(Data$y[-c(1:Data$L.P[Data$P])],
        mu[-c(1:Data$L.P[Data$P])], sigma2[-c(1:Data$L.P[Data$P])],
        log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + delta.prior + omega.prior +
        theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
        return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,P+2), 1, rep(0.5,Q))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("delta -- true:", true.delta,
    " post. mean:", round(fit$Summary2[pos.delta, "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
for (q in 1:Q)
    cat("theta[", q, "] -- true:", true.theta[q],
        " post. mean:", round(fit$Summary2[pos.theta[q], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.omega] <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### AR(p)-GARCH(1,1)

The AR(p)-GARCH(1,1) model extends ARCH by adding a lagged conditional
variance term \\\sigma^2\_{t-1}\\ to the variance equation, providing a
more parsimonious representation of persistent volatility. The
GARCH(1,1) specification often fits financial returns better than
high-order ARCH models. The GARCH(1,1) model was introduced by
Bollerslev [\[15\]](#ref15) as a parsimonious extension of ARCH that
captures persistent volatility dynamics. GARCH models dominate
volatility forecasting in foreign exchange markets and are used by
central banks for risk management.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p}, \quad
t=1,\dots,T\\ \\\epsilon_t = \textbf{y}\_t - \mu_t\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\phi_p \sim \mathcal{N}(0, 1000), \quad
p=1,\dots,P\\ \\\sigma^2_t = \theta_1 + \theta_2 \epsilon^2\_{t-1} +
\theta_3 \sigma^2\_{t-1}\\ \\\omega \sim \mathcal{HC}(25)\\ \\\theta_k =
\frac{1}{1 + \exp(-\theta_k)}, \quad k=1,\dots,3\\ \\\theta_k \sim
\mathcal{N}(0, 1000) \in \[-10,10\], \quad k=1,\dots,3\\

#### model_spec() notation

The AR(p)-GARCH(1,1) model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because both the mean and variance equations involve recursive time
indexing (\\\mu_t\\ depends on \\\textbf{y}\_{t-p}\\, and \\\sigma^2_t\\
depends on both \\\epsilon^2\_{t-1}\\ and \\\sigma^2\_{t-1}\\), which
requires imperative loop constructs that fall outside the declarative
scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42109)
T <- 260
L <- c(1, 5)           # Autoregressive lags
P <- length(L)          # Autoregressive order

#### True parameter values
true.alpha  <- 0.001
true.phi    <- c(0.10, -0.05)
true.omega  <- 0.00005
true.theta1 <- 0.15     # ARCH coefficient
true.theta2 <- 0.75     # GARCH coefficient (theta1 + theta2 < 1)

#### Simulate AR(p)-GARCH(1,1) process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
y[1:max(L)] <- rnorm(max(L), true.alpha, sqrt(true.omega / (1 - true.theta1 - true.theta2)))
eps[1:max(L)] <- y[1:max(L)] - true.alpha
sigma2[1:max(L)] <- true.omega / (1 - true.theta1 - true.theta2)
for (t in (max(L) + 1):T) {
    sigma2[t] <- true.omega + true.theta1 * eps[t - 1]^2 + true.theta2 * sigma2[t - 1]
    mu.t <- true.alpha
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), omega=0,
    theta=rep(0,2)))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(2, 1e-10, 1-1e-5)
    return(c(alpha, phi, omega, theta))
    }
Data <- list(L=L, P=P, PGF=PGF, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- interval(parm[Data$pos.theta], 1e-10, 1-1e-5)
    if(sum(theta) >= 1) theta[2] <- 1 - 1e-5 - theta[1]
    parm[Data$pos.theta] <- theta
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 0, 1, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L[p])] <- mu[-c(1:Data$L[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L[p])]
    epsilon <- Data$y - mu
    sigma2 <- c(omega, omega + theta[1]*epsilon[-Data$T]^ 2)
    sigma2[-1] <- sigma2[-1] + theta[2]*sigma2[-Data$T]
    LL <- sum(dnormv(Data$y[-c(1:Data$L[Data$P])],
        mu[-c(1:Data$L[Data$P])], sigma2[-c(1:Data$L[Data$P])],
    log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + omega.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,P), rep(0.4,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
cat("theta[1] (ARCH) -- true:", true.theta1,
    " post. mean:", round(fit$Summary2[pos.theta[1], "Mean"], 3), "\n")
cat("theta[2] (GARCH) -- true:", true.theta2,
    " post. mean:", round(fit$Summary2[pos.theta[2], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha]  <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.omega]  <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]    <- true.phi[seq_along(pos.phi)]
ground_truth[grep("theta1", parm.names)] <- true.theta1[seq_along(grep("theta1", parm.names))]
ground_truth[grep("theta2", parm.names)] <- true.theta2[seq_along(grep("theta2", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### AR(p)-GARCH(1,1)-M

The GARCH-in-mean variant includes the conditional variance in the mean
equation, which is useful for modelling risk premia in financial
returns. The GARCH-in-mean model combines Bollerslev’s (1986) volatility
specification with the risk premium mechanism of Engle et
al. [\[37\]](#ref37). GARCH-M models are applied in portfolio allocation
to estimate how expected bond returns vary with interest rate
volatility.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=1,\dots,T\\
\\\mu_t = \alpha + \sum^P\_{p=1} \phi_p \textbf{y}\_{t-p} + \delta
\sigma^2\_{t-1}, \quad t=1,\dots,(T+1)\\ \\\epsilon_t = \textbf{y}\_t -
\mu_t\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi_p \sim
\mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\sigma^2_t = \omega +
\theta_1 \epsilon^2\_{t-1} + \theta_2 \sigma^2\_{t-1}\\ \\\omega \sim
\mathcal{HC}(25)\\ \\\theta_k \sim \mathcal{U}(0, 1), \quad
k=1,\dots,2\\

#### model_spec() notation

The AR(p)-GARCH(1,1)-M model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because both the mean and variance equations involve recursive time
indexing (\\\mu_t\\ depends on \\\textbf{y}\_{t-p}\\ and \\\sigma^2_t\\
depends on both \\\epsilon^2\_{t-1}\\ and \\\sigma^2\_{t-1}\\), and the
in-mean term feeds the conditional variance back into the mean equation.
These imperative loop constructs fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42110)
T <- 260
L <- c(1, 5)           # Autoregressive lags
P <- length(L)          # Autoregressive order

#### True parameter values
true.alpha  <- 0.001
true.delta  <- 0.4
true.phi    <- c(0.08, -0.05)
true.omega  <- 0.00005
true.theta1 <- 0.12     # ARCH coefficient
true.theta2 <- 0.78     # GARCH coefficient

#### Simulate AR(p)-GARCH(1,1)-M process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
uncond.var <- true.omega / (1 - true.theta1 - true.theta2)
y[1:max(L)] <- rnorm(max(L), true.alpha, sqrt(uncond.var))
eps[1:max(L)] <- y[1:max(L)] - true.alpha
sigma2[1:max(L)] <- uncond.var
for (t in (max(L) + 1):T) {
    sigma2[t] <- true.omega + true.theta1 * eps[t - 1]^2 + true.theta2 * sigma2[t - 1]
    mu.t <- true.alpha + true.delta * sigma2[t]
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, delta=0, phi=rep(0,P), omega=0,
    theta=rep(0,2)))
pos.alpha <- grep("alpha", parm.names)
pos.delta <- grep("delta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    delta <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(2, 1e-10, 1-1e-5)
    return(c(alpha, delta, phi, omega, theta))
    }
Data <- list(L=L, P=P, PGF=PGF, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.delta=pos.delta,
    pos.phi=pos.phi, pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    delta <- parm[Data$pos.delta]
    phi <- parm[Data$pos.phi]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- interval(parm[Data$pos.theta], 1e-10, 1-1e-5)
    if(sum(theta) >= 1) theta[2] <- 1 - 1e-5 - theta[1]
    parm[Data$pos.theta] <- theta
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    delta.prior <- dnormv(delta, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 0, 1, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L[p])] <- mu[-c(1:Data$L[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L[p])]
    epsilon <- Data$y - mu
    sigma2 <- c(omega, omega + theta[1]*epsilon[-Data$T]^ 2)
    sigma2[-1] <- sigma2[-1] + theta[2]*sigma2[-Data$T]
    mu <- mu + delta*sigma2
    LL <- sum(dnormv(Data$y[-c(1:Data$L[Data$P])],
        mu[-c(1:Data$L[Data$P])], sigma2[-c(1:Data$L[Data$P])],
        log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + delta.prior + phi.prior + omega.prior +
        theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,2), rep(0,P), rep(0.4,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("delta -- true:", true.delta,
    " post. mean:", round(fit$Summary2[pos.delta, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
cat("theta[1] (ARCH) -- true:", true.theta1,
    " post. mean:", round(fit$Summary2[pos.theta[1], "Mean"], 3), "\n")
cat("theta[2] (GARCH) -- true:", true.theta2,
    " post. mean:", round(fit$Summary2[pos.theta[2], "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha]  <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.delta]  <- true.delta[seq_along(pos.delta)]
ground_truth[pos.omega]  <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]    <- true.phi[seq_along(pos.phi)]
ground_truth[grep("theta1", parm.names)] <- true.theta1[seq_along(grep("theta1", parm.names))]
ground_truth[grep("theta2", parm.names)] <- true.theta2[seq_along(grep("theta2", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### AR(p)-TARCH(q)

The threshold ARCH (TARCH) model allows positive and negative shocks to
have asymmetric effects on conditional variance through separate
coefficients \\\theta\_{q,1}\\ for positive and \\\theta\_{q,2}\\ for
negative innovations. This captures the leverage effect observed in
equity markets where negative returns tend to increase volatility more
than positive returns. The Threshold ARCH model was proposed by Zakoian
[\[118\]](#ref118), with a related asymmetric formulation by Glosten et
al. [\[51\]](#ref51). TARCH models capture the leverage effect in equity
markets, where negative returns increase future volatility more than
positive returns of equal magnitude.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=2,\dots,T\\
\\\mu_t = \alpha + \phi^P\_{p=1} \textbf{y}\_{t-p}, \quad
t=(p+1),\dots,T\\ \\\epsilon = \textbf{y} - \mu\\ \\\delta = (\epsilon
\> 0) \times 1\\ \\\sigma^2_t = \omega + \sum^Q\_{q=1} \theta\_{q,1}
\delta\_{t-1} \epsilon^2\_{t-1} + \theta\_{q,2} (1-\delta\_{t-1})
\epsilon^2\_{t-1}\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\phi_p \sim
\mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\omega \sim
\mathcal{HC}(25)\\ \\\theta\_{q,j} \sim \mathcal{U}(0, 1), \quad
q=1\dots,Q, \quad j=1,\dots,2\\

#### model_spec() notation

The AR(p)-TARCH(q) model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the variance equation involves recursive time indexing with
an asymmetric indicator function (\\\delta_t = \mathbb{1}(\epsilon_t \>
0)\\ applied to lagged squared residuals), which requires imperative
loop constructs that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42111)
T <- 260
L.P <- c(1, 5)        # Autoregressive lags
L.Q <- c(1, 2)        # Volatility lags
P <- length(L.P)      # Autoregressive order
Q <- length(L.Q)      # Volatility order

#### True parameter values
true.alpha <- 0.001
true.phi   <- c(0.10, -0.06)
true.omega <- 0.0001
true.theta <- matrix(c(0.10, 0.08,   # theta[,1]: positive shock coefficients
                        0.30, 0.20),  # theta[,2]: negative shock coefficients
                     Q, 2)

#### Simulate AR(p)-TARCH(q) process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
y[1:max(L.P)] <- rnorm(max(L.P), true.alpha, sqrt(true.omega))
eps[1:max(L.P)] <- y[1:max(L.P)] - true.alpha
sigma2[1:max(L.P)] <- true.omega
for (t in (max(L.P) + 1):T) {
    sigma2[t] <- true.omega
    for (q in 1:Q) {
        if (t > L.Q[q]) {
            d <- as.numeric(eps[t - L.Q[q]] > 0)
            sigma2[t] <- sigma2[t] +
                d * true.theta[q, 1] * eps[t - L.Q[q]]^2 +
                (1 - d) * true.theta[q, 2] * eps[t - L.Q[q]]^2
        }
    }
    mu.t <- true.alpha
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L.P[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, phi=rep(0,P), omega=0,
    theta=matrix(0,Q,2)))
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    phi <- runif(Data$P,-1,1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(Data$Q*2, 1e-10, 1-1e-5)
    return(c(alpha, phi, omega, theta))
    }
Data <- list(L.P=L.P, L.Q=L.Q, PGF=PGF, P=P, Q=Q, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    phi <- parm[Data$pos.phi]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- matrix(interval(parm[Data$pos.theta], 1e-10, 1-1e-5), Data$Q,
        2)
    parm[Data$pos.theta] <- as.vector(theta)
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 1e-10, 1-1e-5, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L.P[p])] <- mu[-c(1:Data$L.P[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L.P[p])]
    epsilon <- Data$y - mu
    delta <- (epsilon > 0) * 1
    sigma2 <- rep(omega, Data$T)
    for (q in 1:Data$Q)
        sigma2[-c(1:Data$L.Q[q])] <- sigma2[-c(1:Data$L.Q[q])] +
    delta[1:(Data$T-Data$L.Q[q])] * theta[q,1] *
    epsilon[1:(Data$T-Data$L.Q[q])]^ 2 +
    (1 - delta[1:(Data$T-Data$L.Q[q])]) * theta[q,2] *
    epsilon[1:(Data$T-Data$L.Q[q])]^ 2
    LL <- sum(dnormv(Data$y[-c(1:Data$L.P[Data$P])],
        mu[-c(1:Data$L.P[Data$P])], sigma2[-c(1:Data$L.P[Data$P])],
        log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + phi.prior + omega.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0,P), 1, rep(0.5,Q*2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
cat("theta (true):\n"); print(true.theta)
post.theta <- matrix(fit$Summary2[pos.theta, "Mean"], Q, 2)
cat("theta (post. mean):\n"); print(round(post.theta, 4))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.omega] <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### AR(p)-TARCH(q)-M

The TARCH-in-mean variant adds the conditional variance to the mean
equation, combining asymmetric volatility response with a risk premium.
The mean equation includes \\\gamma_1 \sigma^2_t\\ when the previous
innovation was positive and \\\gamma_2 \sigma^2_t\\ when it was
negative, allowing the risk premium itself to be sign-dependent. The
TARCH-in-mean specification extends the asymmetric volatility model of
Zakoian [\[118\]](#ref118) with a direct mean effect. In energy markets,
TARCH-M models capture how asymmetric oil price volatility feeds back
into expected returns on energy stocks.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_t), \quad t=2,\dots,T\\
\\\mu_t = \alpha + \phi^P\_{p=1} \textbf{y}\_{t-p} + \delta\_{t-1}
\gamma_1 \sigma^2\_{t-1} + (1 - \delta\_{t-1}) \gamma_2 \sigma^2\_{t-1},
\quad t=(p+1),\dots,T\\ \\\epsilon = \textbf{y} - \mu\\ \\\delta =
(\epsilon \> 0) \times 1\\ \\\sigma^2_t = \omega + \sum^Q\_{q=1}
\theta\_{q,1} \delta\_{t-1} \epsilon^2\_{t-1} + \theta\_{q,2}
(1-\delta\_{t-1}) \epsilon^2\_{t-1}\\ \\\alpha \sim \mathcal{N}(0,
1000)\\ \\\gamma_k \sim \mathcal{N}(0, 1000), \quad k=1,\dots,2\\
\\\phi_p \sim \mathcal{N}(0, 1000), \quad p=1,\dots,P\\ \\\omega \sim
\mathcal{HC}(25)\\ \\\theta\_{q,j} \sim \mathcal{U}(0, 1), \quad
q=1\dots,Q, \quad j=1,\dots,2\\

#### model_spec() notation

The AR(p)-TARCH(q)-M model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the variance equation involves recursive time indexing with
an asymmetric indicator function and the in-mean term feeds the
conditional variance back into the mean equation with sign-dependent
coefficients. These imperative loop constructs fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42112)
T <- 260
L.P <- c(1, 5)        # Autoregressive lags
L.Q <- c(1, 2)        # Volatility lags
P <- length(L.P)      # Autoregressive order
Q <- length(L.Q)      # Volatility order

#### True parameter values
true.alpha <- 0.001
true.gamma <- c(0.3, 0.8)   # Risk premium: positive vs negative shocks
true.phi   <- c(0.10, -0.05)
true.omega <- 0.0001
true.theta <- matrix(c(0.10, 0.08,   # theta[,1]: positive shock coefficients
                        0.30, 0.20),  # theta[,2]: negative shock coefficients
                     Q, 2)

#### Simulate AR(p)-TARCH(q)-M process
y <- numeric(T)
eps <- numeric(T)
sigma2 <- numeric(T)
y[1:max(L.P)] <- rnorm(max(L.P), true.alpha, sqrt(true.omega))
eps[1:max(L.P)] <- y[1:max(L.P)] - true.alpha
sigma2[1:max(L.P)] <- true.omega
for (t in (max(L.P) + 1):T) {
    sigma2[t] <- true.omega
    for (q in 1:Q) {
        if (t > L.Q[q]) {
            d <- as.numeric(eps[t - L.Q[q]] > 0)
            sigma2[t] <- sigma2[t] +
                d * true.theta[q, 1] * eps[t - L.Q[q]]^2 +
                (1 - d) * true.theta[q, 2] * eps[t - L.Q[q]]^2
        }
    }
    d.prev <- as.numeric(eps[t - 1] > 0)
    mu.t <- true.alpha +
        d.prev * true.gamma[1] * sigma2[t] +
        (1 - d.prev) * true.gamma[2] * sigma2[t]
    for (p in 1:P) mu.t <- mu.t + true.phi[p] * y[t - L.P[p]]
    y[t] <- rnorm(1, mu.t, sqrt(sigma2[t]))
    eps[t] <- y[t] - mu.t
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, gamma=rep(0,2), phi=rep(0,P),
    omega=0, theta=matrix(0,Q,2)))
pos.alpha <- grep("alpha", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.phi <- grep("phi", parm.names)
pos.omega <- grep("omega", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    gamma <- rnorm(2)
    phi <- runif(Data$P,-1,1)
    omega <- rhalfcauchy(1,5)
    theta <- runif(Data$Q*2, 1e-10, 1-1e-5)
    return(c(alpha, gamma, phi, omega, theta))
    }
Data <- list(L.P=L.P, L.Q=L.Q, PGF=PGF, P=P, Q=Q, T=T, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.gamma=pos.gamma,
    pos.phi=pos.phi, pos.omega=pos.omega, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    gamma <- parm[Data$pos.gamma]
    phi <- parm[Data$pos.phi]
    omega <- interval(parm[Data$pos.omega], 1e-100, Inf)
    parm[Data$pos.omega] <- omega
    theta <- matrix(interval(parm[Data$pos.theta], 1e-10, 1-1e-5), Data$Q,
        2)
    parm[Data$pos.theta] <- as.vector(theta)
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    phi.prior <- sum(dnormv(phi, 0, 1000, log=TRUE))
    omega.prior <- dhalfcauchy(omega, 25, log=TRUE)
    theta.prior <- sum(dunif(theta, 1e-10, 1-1e-5, log=TRUE))
    ### Log-Likelihood
    mu <- rep(alpha, Data$T)
    for (p in 1:Data$P)
        mu[-c(1:Data$L.P[p])] <- mu[-c(1:Data$L.P[p])] +
    phi[p]*Data$y[1:(Data$T-Data$L.P[p])]
    epsilon <- Data$y - mu
    delta <- (epsilon > 0) * 1
    sigma2 <- rep(omega, Data$T)
    for (q in 1:Data$Q)
        sigma2[-c(1:Data$L.Q[q])] <- sigma2[-c(1:Data$L.Q[q])] +
    delta[1:(Data$T-Data$L.Q[q])] * theta[q,1] *
    epsilon[1:(Data$T-Data$L.Q[q])]^ 2 +
    (1 - delta[1:(Data$T-Data$L.Q[q])]) * theta[q,2] *
    epsilon[1:(Data$T-Data$L.Q[q])]^ 2
    mu <- mu + delta*gamma[1]*sigma2 + (1 - delta)*gamma[2]*sigma2
    LL <- sum(dnormv(Data$y[-c(1:Data$L.P[Data$P])],
        mu[-c(1:Data$L.P[Data$P])], sigma2[-c(1:Data$L.P[Data$P])],
        log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + gamma.prior + phi.prior + omega.prior +
        theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(mu), mu, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,3), rep(0,P), 1, rep(0.5,Q*2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("gamma[1] -- true:", true.gamma[1],
    " post. mean:", round(fit$Summary2[pos.gamma[1], "Mean"], 3), "\n")
cat("gamma[2] -- true:", true.gamma[2],
    " post. mean:", round(fit$Summary2[pos.gamma[2], "Mean"], 3), "\n")
for (p in 1:P)
    cat("phi[", p, "] -- true:", true.phi[p],
        " post. mean:", round(fit$Summary2[pos.phi[p], "Mean"], 3), "\n")
cat("omega -- true:", true.omega,
    " post. mean:", round(fit$Summary2[pos.omega, "Mean"], 3), "\n")
cat("theta (true):\n"); print(true.theta)
post.theta <- matrix(fit$Summary2[pos.theta, "Mean"], Q, 2)
cat("theta (post. mean):\n"); print(round(post.theta, 4))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.omega] <- true.omega[seq_along(pos.omega)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Time series: vector autoregression

### PVAR(p)

The Poisson vector autoregression of order \\p\\ models multivariate
count time series by combining a log-linear Poisson likelihood with VAR
dynamics in the latent intensity. Unit-specific random effects
\\\gamma\_{t,j}\\ account for overdispersion, and a Wishart prior on the
precision matrix \\\Omega\\ induces correlation across the \\J\\ series.
Panel vector autoregression was developed by Canova and Ciccarelli
[\[18\]](#ref18) for analysing dynamic interdependencies across
cross-sectional units. Panel VAR models are used in international
economics to study the transmission of monetary policy shocks across
countries.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{P}(\lambda\_{t,j}), \quad t=1,\dots,T
\quad j=1,\dots,J\\ \\\lambda\_{t,j} = \sum^P\_{p=1} \Phi\_{1:J,j,p}
\textbf{Y}\_{t-p,j} + \exp(\alpha_j + \gamma\_{t,j})\\ \\\alpha_j \sim
\mathcal{N}(0, 1000)\\ \\\Phi\_{i,k,p} \sim
\mathcal{N}(\Phi^\mu\_{i,k,p}, \Sigma\_{i,k,p}), \quad i=1,\dots,J,
\quad k=1,\dots,J, \quad p=1,\dots,P\\ \\\gamma\_{t,1:J} \sim
\mathcal{N}\_J(0, \Omega^{-1}), \quad t=1,\dots,T\\ \\\Omega \sim
\mathcal{W}\_{J+1}(\textbf{S}), \quad \textbf{S} = \textbf{I}\_J\\ where
\\\Phi^\mu\\ and \\\Sigma\\ are set according to the Minnesota prior.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it combines Poisson likelihoods with VAR lag structure,
Minnesota priors on the coefficient arrays, and a Wishart-distributed
precision matrix, all of which require imperative array construction
that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42208)
T <- 50
J <- 3
L <- c(1, 2)
P <- length(L)

#### Simulate count panel from Poisson VAR
alpha.true <- c(1.5, 1.0, 0.8)
Y <- matrix(0, T, J)
Y[1:max(L), ] <- matrix(rpois(max(L) * J, exp(alpha.true)), max(L), J, byrow = TRUE)
for (t in (max(L) + 1):T) {
    lambda.t <- exp(alpha.true + rnorm(J, 0, 0.3))
    for (p in 1:P) lambda.t <- lambda.t + 0.05 * Y[t - L[p], ]
    Y[t, ] <- rpois(J, pmax(lambda.t, 0.1))
}

#### Assemble Data list
Phi.mu <- array(0, dim=c(J,J,P))
Phi.mu[, , 1] <- diag(J)
S <- diag(J)
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J),
    Phi=array(0, dim=c(J,J,P)), gamma=matrix(0,T-L[P],J), U=S),
    uppertri=c(0,0,0,1))
pos.alpha <- grep("alpha", parm.names)
pos.Phi <- grep("Phi", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    Phi <- runif(Data$J*Data$J*Data$P, -1e-10, 1e-10)
    gamma <- rnorm((Data$T-Data$L[Data$P])*Data$J)
    U <- rwishartc(Data$J+1, diag(Data$J))
    return(c(alpha, Phi, gamma, U[upper.tri(U, diag=TRUE)]))
    }
Data <- list(J=J, L=L, P=P, PGF=PGF, Phi.mu=Phi.mu, S=S, T=T, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.Phi=pos.Phi, pos.gamma=pos.gamma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    Phi <- array(parm[Data$pos.Phi], dim=c(Data$J, Data$J, Data$P))
    gamma <- matrix(parm[Data$pos.gamma], Data$T-Data$L[Data$P], Data$J)
    U <- as.parm.matrix(U, Data$J, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    Omega <- t(U) %*% U
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    Sigma <- MinnesotaPrior(Data$J, lags=Data$L, lambda=1, theta=0.5,
        diag(as.inverse(Omega)))
    Phi.prior <- sum(dnormv(Phi, Data$Phi.mu, Sigma, log=TRUE))
    gamma.prior <- sum(dmvnp(gamma, 0, Omega, log=TRUE))
    U.prior <- dwishart(Omega, Data$J+1, Data$S, log=TRUE)
    ### Log-Likelihood
    lambda <- exp(matrix(alpha, Data$T, Data$J, byrow=TRUE) +
        rbind(matrix(0, Data$L[Data$P], Data$J), gamma))
    for (p in 1:Data$P)
        lambda[(1+Data$L[p]):Data$T,] <- lambda[(1+Data$L[p]):Data$T,] +
    Data$Y[1:(Data$T-Data$L[p]),] %*% Phi[, , p]
    LL <- sum(dpois(Data$Y[(1+Data$L[Data$P]):Data$T,],
        lambda[(1+Data$L[Data$P]):Data$T,], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + Phi.prior + gamma.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(prod(dim(lambda)), lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(0,J*J*P), rep(0,(T-L[P])*J),
    rep(0,J*(J+1)/2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### VAR(p) - Minnesota Prior

This is a vector autoregression of order \\p\\ with the Minnesota prior
applied to the coefficient arrays \\\Phi\\. The Minnesota prior shrinks
off-diagonal and higher-lag coefficients more strongly toward zero,
encoding the belief that each series is best predicted by its own recent
past, with cross-series effects and distant lags having diminishing
influence. The Minnesota prior for vector autoregressions was introduced
by Litterman [\[73\]](#ref73) and Doan et al. [\[31\]](#ref31) to
regularise the large number of VAR parameters. Minnesota-prior VARs are
the standard tool at central banks for producing macroeconomic
forecasts.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{N}(\mu\_{t,j}, \sigma^2_j), \quad
t=1,\dots,T, \quad j=1,\dots,J\\ \\\mu\_{t,j} = \alpha_j + \sum^P\_{p=1}
\Phi\_{1:J,j,p} \textbf{Y}\_{t-p,j}\\ \\\textbf{y}^{new}\_j = \alpha_j +
\Phi\_{1:J,j} \textbf{Y}\_{T,j}\\ \\\alpha_j \sim \mathcal{N}(0, 1000)\\
\\\Phi\_{i,k,p} \sim \mathcal{N}(\Phi^\mu\_{i,k,p}, \Sigma\_{i,k,p}),
\quad i=1,\dots,J, \quad k=1,\dots,J, \quad p=1,\dots,P\\ \\\sigma_j
\sim \mathcal{HC}(25)\\ where \\\Phi^\mu\\ and \\\Sigma\\ are set
according to the Minnesota prior.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the Minnesota prior requires constructing a structured prior
variance array \\\Sigma\\ that depends on lag positions, own-
vs. cross-series status, and residual variances, which falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42212)
T <- 100
J <- 3
L <- c(1, 5)
P <- length(L)

#### Simulate a VAR process
Y <- matrix(rnorm(T * J, 0, 0.5), T, J)
for (t in (max(L) + 1):T) {
    Y[t, ] <- Y[t, ] + 0.3 * Y[t - 1, ] - 0.1 * Y[t - 5, ]
}
Y.scales <- sqrt(apply(Y, 2, var))
Y <- Y / matrix(Y.scales, T, J, byrow = TRUE)

#### Assemble Data list
Phi.mu <- array(0, dim=c(J,J,P))
Phi.mu[, , 1] <- diag(J)
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J),
    Phi=array(0, dim=c(J,J,P)), sigma=rep(0,J)))
pos.alpha <- grep("alpha", parm.names)
pos.Phi <- grep("Phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    Phi <- runif(Data$J*Data$J*Data$P, -1, 1)
    sigma <- runif(Data$J)
    return(c(alpha, Phi, sigma))
    }
Data <- list(J=J, L=L, P=P, PGF=PGF, Phi.mu=Phi.mu, T=T, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.Phi=pos.Phi, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    Phi <- array(parm[Data$pos.Phi], dim=c(Data$J, Data$J, Data$P))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    Sigma <- MinnesotaPrior(Data$J, lags=Data$L, lambda=1, theta=0.5,
        sigma)
    Phi.prior <- sum(dnormv(Phi, Data$Phi.mu, Sigma, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(alpha, Data$T, Data$J, byrow=TRUE)
    for (p in 1:Data$P) {
        mu[(1+Data$L[p]):Data$T,] <- mu[(1+Data$L[p]):Data$T,] +
        Data$Y[1:(Data$T-Data$L[p]),] %*% Phi[ , , p]}
    Sigma <- matrix(sigma, Data$T, Data$J, byrow=TRUE)
    LL <- sum(dnorm(Data$Y[(1+Data$L[Data$P]):Data$T,],
        mu[(1+Data$L[Data$P]):Data$T,],
        Sigma[(1+Data$L[Data$P]):Data$T,], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + Phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(as.vector(colMeans(Y)), rep(0,J*J*P), rep(1,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### VAR(p) - SSVS

Stochastic search variable selection (SSVS) is applied to the VAR
autoregressive coefficients using Bernoulli indicator variables
\\\Gamma\\ that switch each coefficient between a tight spike
(\\\sigma=0.1\\) and a diffuse slab (\\\sigma=10\\). This allows
automatic identification of which cross-series and lagged dependencies
are supported by the data. VAR models with SSVS priors combine
Lutkepohl’s (2005) VAR framework with the stochastic search approach of
George et al. [\[50\]](#ref50) for parsimonious lag selection. VAR-SSVS
models are applied in monetary economics to identify which lagged
variables carry predictive power for interest rate dynamics.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{N}(\mu\_{t,j}, \sigma^2_j), \quad
t=1,\dots,T, \quad j=1,\dots,J\\ \\\mu\_{t,j} = \alpha_j + \sum^P\_{p=1}
\Gamma\_{1:J,j,p}\Phi\_{1:J,j,p}\textbf{Y}\_{t-p,j}\\ \\\alpha_j \sim
\mathcal{N}(0, 1000)\\ \\\Gamma\_{i,k,p} \sim \mathcal{BERN}(0.5), \quad
i=1,\dots,J, \quad k=1,\dots,J, \quad p=1,\dots,P\\ \\(\Phi\_{i,k,p} \|
\Gamma\_{i,k,p}) \sim (1 - \Gamma\_{i,k,p})\mathcal{N}(0, 0.01) +
\Gamma\_{i,k,p}\mathcal{N}(0, 10), \quad i=1,\dots,J, \quad k=1,\dots,J,
\quad p=1,\dots,P\\ \\\sigma_j \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because SSVS requires Bernoulli indicator variables that switch between
spike and slab prior variances for each VAR coefficient, combined with
multivariate lag structure, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42213)
T <- 100
J <- 3
L <- c(1, 5)
P <- length(L)

#### Simulate a sparse VAR process
Y <- matrix(rnorm(T * J, 0, 0.5), T, J)
for (t in (max(L) + 1):T) {
    Y[t, ] <- Y[t, ] + 0.4 * Y[t - 1, ] - 0.15 * Y[t - 5, ]
}
Y.scales <- sqrt(apply(Y, 2, var))
Y <- Y / matrix(Y.scales, T, J, byrow = TRUE)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J),
     Gamma=array(0, dim=c(J,J,P)), Phi=array(0, dim=c(J,J,P)),
     sigma=rep(0,J)))
pos.alpha <- grep("alpha", parm.names)
pos.Gamma <- grep("Gamma", parm.names)
pos.Phi <- grep("Phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    Gamma <- rep(1, Data$J*Data$J*Data$P)
    Phi <- runif(Data$J*Data$J*Data$P, -1, 1)
    sigma <- runif(Data$J)
    return(c(alpha, Gamma, Phi, sigma))
    }
Data <- list(J=J, L=L, P=P, PGF=PGF, T=T, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.Gamma=pos.Gamma,
    pos.Phi=pos.Phi, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    Gamma <- array(parm[Data$pos.Gamma], dim=c(Data$J, Data$J, Data$P))
    Phi.Sigma <- Gamma * 10
    Phi.Sigma[Gamma == 0] <- 0.1
    Phi <- array(parm[Data$pos.Phi], dim=c(Data$J, Data$J, Data$P))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    Gamma.prior <- sum(dbern(Gamma, 0.5, log=TRUE))
    Phi.prior <- sum(dnorm(Phi, 0, Phi.Sigma, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(alpha, Data$T, Data$J, byrow=TRUE)
    for (p in 1:Data$P)
        mu[(1+Data$L[p]):Data$T,] <- mu[(1+Data$L[p]):Data$T,] +
        Data$Y[1:(Data$T-Data$L[p]),] %*% (Gamma[, , p]*Phi[, , p])
    Sigma <- matrix(sigma, Data$T, Data$J, byrow=TRUE)
    LL <- sum(dnorm(Data$Y[(1+Data$L[Data$P]):Data$T,],
        mu[(1+Data$L[Data$P]):Data$T,],
        Sigma[(1+Data$L[Data$P]):Data$T,], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + Gamma.prior + Phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(colMeans(Y), rep(1,J*J*P), runif(J*J*P,-1,1), rep(1,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### VAR(p)-GARCH(1,1)-M

The VAR(p)-GARCH(1,1)-M model combines a multivariate vector
autoregression with asymmetric BEKK multivariate GARCH and a
GARCH-in-mean term. The Minnesota prior governs the VAR coefficients,
and the conditional covariance matrix \\H_t\\ enters the mean equation
through coefficients \\\delta\\, allowing the risk premium to depend on
time-varying volatility. The VAR-GARCH-in-mean model extends the
multivariate GARCH framework of Engle and Kroner [\[36\]](#ref36) to
vector autoregressions with volatility spillovers. VAR-GARCH-M models
are applied in international finance to study volatility transmission
between stock markets.

#### Form

\\\textbf{Y}\_{t,1:J} \sim \mathcal{N}\_J(\mu\_{t,1:J},
H\_{1:J,1:J,t})\\ \\\mu\_{t,j} = \alpha_j + \sum^P\_{p=1}
\Phi\_{1:J,j,p} \textbf{Y}\_{t-p,j} + \sum \textbf{H}\_{1:J,j,t-1}
\delta\_{1:J,j}\\ \\\textbf{H}\_{,,t} = \Omega + \textbf{A}^T
\epsilon\_{t-1,}\epsilon^T\_{t-1} \textbf{A} + \textbf{B}^T
\textbf{H}\_{,,t-1}\textbf{B} +
\textbf{D}^T\zeta\_{t-1,}\zeta^T\_{t-1,}\textbf{D}, \quad t=2,\dots,T\\
\\\Omega = \textbf{C}\textbf{C}^T\\ \\\alpha_j \sim \mathcal{N}(0,
1000)\\ \\\delta\_{i,k} \sim \mathcal{N}(0, 1000), \quad i=1,\dots,J,
\quad k=1,\dots,J\\ \\\Phi\_{i,k,p} \sim \mathcal{N}(\Phi^\mu\_{i,k,p},
\Sigma\_{i,k,p}), \quad i=1,\dots,J, \quad k=1,\dots,J, \quad
p=1,\dots,P\\ \\\textbf{C}\_{i,j} \sim \mathcal{N}(0, 100)\\
\\\textbf{A}\_{i,j} \sim \mathcal{N}(0, 100)\\ \\\textbf{B}\_{i,j} \sim
\mathcal{N}(0, 100)\\ \\\textbf{D}\_{i,j} \sim \mathcal{N}(0, 100)\\
where \\\Phi\\ has a Minnesota prior, **C** is lower-triangular with
positive-only diagonal elements, and \\\textbf{A}\_{1,1}\\,
\\\textbf{B}\_{1,1}\\, and \\\textbf{D}\_{1,1}\\ must be positive.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the BEKK GARCH recursion for \\H_t\\, the asymmetric leverage
term via \\\zeta\\, and the GARCH-in-mean feedback into \\\mu_t\\ all
require imperative matrix-level computations that fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42211)
T <- 100
J <- 3
L <- c(1, 5)
P <- length(L)

#### Simulate a VAR process with constant variance (simplified)
Y <- matrix(rnorm(T * J, 0, 0.5), T, J)
for (t in (max(L) + 1):T) {
    Y[t, ] <- Y[t, ] + 0.2 * Y[t - 1, ] - 0.05 * Y[t - 5, ]
}
Y.scales <- sqrt(apply(Y, 2, var))
Y <- Y / matrix(Y.scales, T, J, byrow = TRUE)

#### Assemble Data list
Phi.mu <- array(0, dim=c(J,J,P))
Phi.mu[, , 1] <- diag(J)
C <- matrix(NA, J, J)
C[lower.tri(C, diag=TRUE)] <- 0
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J), delta=matrix(0,J,J),
    Phi=array(0, dim=c(J,J,P)), C=C, A=matrix(0,J,J), B=matrix(0,J,J),
    D=matrix(0,J,J)))
pos.alpha <- grep("alpha", parm.names)
pos.delta <- grep("delta", parm.names)
pos.Phi <- grep("Phi", parm.names)
pos.C <- grep("C", parm.names)
pos.A <- grep("A", parm.names)
pos.B <- grep("B", parm.names)
pos.D <- grep("D", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    delta <- rnorm(Data$J*Data$J)
    Phi <- runif(Data$J*Data$J*Data$P, -1, 1)
    C <- runif(Data$J*(Data$J+1)/2)
    A <- as.vector(diag(Data$J)) + runif(Data$J*Data$J, -0.1, 0.1)
    B <- as.vector(diag(Data$J)) + runif(Data$J*Data$J, -0.1, 0.1)
    D <- as.vector(diag(Data$J)) + runif(Data$J*Data$J, -0.1, 0.1)
    return(c(alpha, delta, Phi, C, A, B, D))
    }
Data <- list(J=J, L=L, P=P, PGF=PGF, Phi.mu=Phi.mu, T=T, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.delta=pos.delta, pos.Phi=pos.Phi, pos.C=pos.C, pos.A=pos.A,
    pos.B=pos.B, pos.D=pos.D)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    delta <- matrix(parm[Data$pos.delta], Data$J, Data$J)
    Phi <- array(parm[Data$pos.Phi], dim=c(Data$J, Data$J, Data$P))
    C <- matrix(0, Data$J, Data$J)
    C[lower.tri(C, diag=TRUE)] <- parm[Data$pos.C]
    diag(C) <- abs(diag(C))
    parm[Data$pos.C] <- C[lower.tri(C, diag=TRUE)]
    Omega <- C %*% t(C)
    A <- matrix(parm[Data$pos.A], Data$J, Data$J)
    A[1,1] <- abs(A[1,1])
    parm[Data$pos.A] <- as.vector(A)
    B <- matrix(parm[Data$pos.B], Data$J, Data$J)
    B[1,1] <- abs(B[1,1])
    parm[Data$pos.B] <- as.vector(B)
    D <- matrix(parm[Data$pos.D], Data$J, Data$J)
    D[1,1] <- abs(D[1,1])
    parm[Data$pos.D] <- as.vector(D)
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    delta.prior <- sum(dnormv(delta, 0, 1000, log=TRUE))
    Sigma <- MinnesotaPrior(Data$J, lags=Data$L, lambda=1,
        theta=0.5, sqrt(diag(Omega)))
    Phi.prior <- sum(dnormv(Phi, Data$Phi.mu, Sigma, log=TRUE))
    C.prior <- sum(dnormv(C[lower.tri(C, diag=TRUE)], 0, 100, log=TRUE))
    A.prior <- sum(dnormv(A, 0, 100, log=TRUE))
    B.prior <- sum(dnormv(B, 0, 100, log=TRUE))
    D.prior <- sum(dnormv(D, 0, 100, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(alpha, Data$T, Data$J, byrow=TRUE)
    for (p in 1:Data$P)
        mu[(1+Data$L[p]):Data$T,] <- mu[(1+Data$L[p]):Data$T,] +
    Data$Y[1:(Data$T-Data$L[p]),] %*% Phi[, , p]
    LL <- 0
    Yhat <- Data$Y
    H <- array(Omega, dim=c(Data$J, Data$J, Data$T))
    for (t in 2:Data$T) {
        eps <- Data$Y - mu
        zeta <- matrix(interval(eps, -Inf, 0, reflect=FALSE), Data$T,
    Data$J)
        part1 <- t(A) %*% eps[t-1,] %*% t(eps[t-1,]) %*% A
        part2 <- t(B) %*% H[, , t-1] %*% B
        part3 <- t(D) %*% zeta[t-1,] %*% t(zeta[t-1,]) %*% D
        H0 <- Omega + part1 + part2 + part3
        H0[upper.tri(H0, diag=TRUE)] <- t(H0)[upper.tri(H0, diag=TRUE)]
        H[, , t] <- H0
        mu[t-1,] <- mu[t-1,] + colMeans(H[, , t-1]*delta)
        Sigma <- MinnesotaPrior(Data$J, lags=Data$L, lambda=1,
    theta=0.5, sqrt(diag(H[, , t])))
        Phi.prior <- Phi.prior + sum(dnormv(Phi, Data$Phi.mu, Sigma,
    log=TRUE))
        LL <- LL + dmvn(Y[t,], mu[t,], H[, , t], log=TRUE)
        Yhat[t,] <- rmvn(1, mu[t,], H[, , t])
        }
    Phi.prior <- Phi.prior / Data$T
    ### Log-Posterior
    LP <- LL + alpha.prior + delta.prior + Phi.prior + C.prior +
        A.prior + B.prior + D.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=Yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(colMeans(Y), rnorm(J*J), runif(J*J*P,-1,1),
    runif(J*(J+1)/2), as.vector(diag(J)), as.vector(diag(J)),
    as.vector(diag(J)))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### VARMA(p,q) - SSVS

Stochastic search variable selection (SSVS) is applied to VARMA
parameters, using spike-and-slab priors on both the autoregressive
\\\Phi\\ and moving-average \\\Theta\\ coefficient arrays. When
\\\Gamma=1\\ the coefficient receives a diffuse slab prior; when
\\\Gamma=0\\ it is shrunk toward zero by a tight spike. Note that the
mixture variance constants are typically multiplied by posterior
standard deviations from an unrestricted VARMA estimated previously;
this example omits that calibration step for simplicity. VARMA models
with stochastic search variable selection combine the multivariate time
series framework of Lütkepohl [\[77\]](#ref77) with the spike-and-slab
approach of George et al. [\[50\]](#ref50). VARMA-SSVS models are used
in neuroscience to identify directed connectivity patterns among brain
regions from fMRI time series.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{N}(\mu\_{t,j}, \sigma^2_j), \quad
t=1,\dots,T, \quad j=1,\dots,J\\ \\\mu\_{t,j} = \alpha_j + \sum^P\_{p=1}
\Gamma^\Phi\_{1:J,j,p}\Phi\_{1:J,j,p}\textbf{Y}\_{t-p,j} + \sum^Q\_{q=1}
\Gamma^\Theta\_{1:J,j,q}\Theta\_{1:J,j,q} \epsilon\_{t-q,j}\\ \\\alpha_j
\sim \mathcal{N}(0, 1000)\\ \\\Gamma^\Phi\_{i,k,p} \sim
\mathcal{BERN}(0.5), \quad i=1,\dots,J, \quad k=1,\dots,J, \quad
p=1,\dots,P\\ \\(\Phi\_{i,k,p} \| \Gamma^\Phi\_{i,k,p}) \sim (1 -
\Gamma^\Phi\_{i,k,p})\mathcal{N}(0, 0.01) +
\Gamma^\Phi\_{i,k,p}\mathcal{N}(0, 10), \quad i=1,\dots,J, \quad
k=1,\dots,J, \quad p=1,\dots,P\\ \\\Gamma^\Theta\_{i,k,q} \sim
\mathcal{BERN}(0.5), \quad i=1,\dots,J, \quad k=1,\dots,J, \quad
q=1,\dots,Q\\ \\(\Theta\_{i,k,q} \| \Gamma^\Theta\_{i,k,q}) \sim (1 -
\Gamma^\Theta\_{i,k,q})\mathcal{N}(0, 0.01) +
\Gamma^\Theta\_{i,k,q}\mathcal{N}(0, 10), \quad i=1,\dots,J, \quad
k=1,\dots,J, \quad q=1,\dots,Q\\ \\\sigma_j \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because SSVS requires Bernoulli indicator variables that switch between
spike and slab prior variances, combined with VARMA lag structure and
moving-average residual feedback, all of which fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42210)
T <- 100
J <- 3
L.P <- c(1, 5)
L.Q <- c(1, 2)
P <- length(L.P)
Q <- length(L.Q)

#### Simulate a VAR process (simplified; MA component is zero in truth)
Y <- matrix(rnorm(T * J, 0, 0.5), T, J)
for (t in (max(L.P) + 1):T) {
    Y[t, ] <- Y[t, ] + 0.3 * Y[t - 1, ] - 0.1 * Y[t - 5, ]
}
Y.scales <- sqrt(apply(Y, 2, var))
Y <- Y / matrix(Y.scales, T, J, byrow = TRUE)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,J),
    Gamma.phi=array(0, dim=c(J,J,P)), Phi=array(0, dim=c(J,J,P)),
    Gamma.theta=array(0, dim=c(J,J,Q)), Theta=array(0, dim=c(J,J,Q)),
    sigma=rep(0,J)))
pos.alpha <- grep("alpha", parm.names)
pos.Gamma.phi <- grep("Gamma.phi", parm.names)
pos.Phi <- grep("Phi", parm.names)
pos.Gamma.theta <- grep("Gamma.theta", parm.names)
pos.Theta <- grep("Theta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$J)
    Gamma.phi <- rep(1, Data$J*Data$J*Data$P)
    Phi <- runif(Data$J*Data$J*Data$P, -1, 1)
    Gamma.theta <- rep(1, Data$J*Data$J*Data$Q)
    Theta <- rnorm(Data$J*Data$J*Data$Q)
    sigma <- runif(Data$J)
    return(c(alpha, Gamma.phi, Phi, Gamma.theta, Theta, sigma))
    }
Data <- list(J=J, L.P=L.P, L.Q=L.Q, P=P, Q=Q, PGF=PGF, T=T, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.alpha=pos.alpha,
    pos.Gamma.phi=pos.Gamma.phi, pos.Phi=pos.Phi,
    pos.Gamma.theta=pos.Gamma.theta, pos.Theta=pos.Theta,
    pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    Gamma.phi <- array(parm[Data$pos.Gamma.phi],
        dim=c(Data$J, Data$J, Data$P))
    Phi.Sigma <- Gamma.phi * 10
    Phi.Sigma[Gamma.phi == 0] <- 0.1
    Phi <- array(parm[Data$pos.Phi], dim=c(Data$J, Data$J, Data$P))
    Gamma.theta <- array(parm[Data$pos.Gamma.theta],
        dim=c(Data$J, Data$J, Data$Q))
    Theta.Sigma <- Gamma.theta * 10
    Theta.Sigma[Gamma.theta == 0] <- 0.1
    Theta <- array(parm[Data$pos.Theta], dim=c(Data$J, Data$J, Data$Q))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 1000, log=TRUE))
    Gamma.phi.prior <- sum(dbern(Gamma.phi, 0.5, log=TRUE))
    Phi.prior <- sum(dnorm(Phi, 0, Phi.Sigma, log=TRUE))
    Gamma.theta.prior <- sum(dbern(Gamma.theta, 0.5, log=TRUE))
    Theta.prior <- sum(dnorm(Theta, 0, Theta.Sigma, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- matrix(alpha, Data$T, Data$J, byrow=TRUE)
    for (p in 1:Data$P)
        mu[(1+Data$L.P[p]):Data$T,] <- mu[(1+Data$L.P[p]):Data$T,] +
    Data$Y[1:(Data$T-Data$L.P[p]),] %*%
    (Gamma.phi[, , p] * Phi[, , p])
    epsilon <- Data$Y - mu
    for (q in 1:Data$Q)
        mu[(1+Data$L.Q[q]):Data$T,] <- mu[(1+Data$L.Q[q]):Data$T,] +
    epsilon[1:(Data$T-Data$L.Q[q]),] %*%
    (Gamma.theta[, , q] * Theta[, , q])
    Sigma <- matrix(sigma, Data$T, Data$J, byrow=TRUE)
    LL <- sum(dnorm(Data$Y[(1+Data$L.P[Data$P]):Data$T,],
        mu[(1+Data$L.P[Data$P]):Data$T,],
        Sigma[(1+Data$L.P[Data$P]):Data$T,], log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + Gamma.phi.prior + Phi.prior +
    Gamma.theta.prior + Theta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(colMeans(Y), rep(1,J*J*P), runif(J*J*P,-1,1),
    rep(1,J*J*Q), rep(0,J*J*Q), rep(1,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## State space, latent process, and change point models

### Brownian Motion with Drift **\[NEW\]**

Brownian motion with drift is the simplest continuous-time stochastic
process, describing a particle whose position increments are independent
Gaussian random variables with a constant drift rate \\\alpha\\ and
diffusion coefficient \\\sigma\\. In its discretised form the transition
density is normal with mean \\y\_{t-1} + \alpha \Delta t\\ and variance
\\\sigma^2 \Delta t\\, yielding a Gaussian random walk with systematic
trend. The physical theory was developed by Einstein [\[124\]](#ref124),
who derived the diffusion equation from first principles and showed that
the mean squared displacement grows linearly with time. Brownian motion
with drift is fundamental in quantitative finance for modelling
geometric asset price evolution under the Black-Scholes framework.

#### Form

\\y_t \mid y\_{t-1} \sim \mathcal{N}(y\_{t-1} + \alpha \\ \Delta t, \\
\sigma^2 \\ \Delta t)\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\sigma
\sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the likelihood conditions each observation on the previous
observation \\y\_{t-1}\\, requiring recursive dependence on the data
vector indexed by \\t\\. The declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
does not support this kind of sequential transition density. The manual
`Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42204)
T <- 300
dt <- 0.1

#### True parameter values
true.alpha <- 0.5
true.sigma <- 1.0

#### Simulate Brownian motion with drift
y <- numeric(T)
y[1] <- 0
for (t in 2:T)
    y[t] <- y[t-1] + true.alpha * dt + true.sigma * sqrt(dt) * rnorm(1)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, sigma=0))
pos.alpha <- grep("alpha", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    sigma <- rhalfcauchy(1, 5)
    return(c(alpha, sigma))
    }
Data <- list(T=T, dt=dt, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- Data$y[-Data$T] + alpha * Data$dt
    sd.inc <- sigma * sqrt(Data$dt)
    LL <- sum(dnorm(Data$y[-1], mu, sd.inc, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + sigma.prior
    yhat <- c(rnorm(Data$T - 1, mu, sd.inc),
        rnorm(1, Data$y[Data$T] + alpha * Data$dt, sd.inc))
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Change Point in Variance **\[NEW\]**

This model detects an unknown point \\\tau\\ in a time series at which
the variance shifts from \\\sigma_1^2\\ to \\\sigma_2^2\\ while the mean
remains constant. The change point is treated as a continuous parameter
constrained to \\\[1, T\]\\, and the floor of its posterior distribution
indicates the most probable transition time. Bayesian change point
methods were developed by Chib [\[24\]](#ref24) using Gibbs sampling for
multiple change points, building on the product-partition framework of
Barry and Hartigan [\[136\]](#ref136) for detecting distributional
shifts. Change point variance models are used in seismology to detect
shifts in background microseismic activity rates that may precede major
tectonic events.

#### Form

\\y_t \sim \mathcal{N}(\mu, \sigma_t^2)\\ \\\sigma_t = \begin{cases}
\sigma_1 & \text{if } t \< \tau \\ \sigma_2 & \text{if } t \geq \tau
\end{cases}\\ \\\mu \sim \mathcal{N}(0, 1000)\\ \\\sigma_1 \sim
\mathcal{HC}(25)\\ \\\sigma_2 \sim \mathcal{HC}(25)\\ \\\tau \sim
\mathcal{U}(1, T)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the observation variance is selected by a piecewise rule
that depends on the estimated change point parameter \\\tau\\. This kind
of parameter-dependent conditional switching between variance components
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42205)
T <- 200

#### True parameter values
true.mu <- 5
true.sigma1 <- 0.5
true.sigma2 <- 2.0
true.tau <- 100

#### Simulate change point in variance
sigma.true <- ifelse(1:T < true.tau, true.sigma1, true.sigma2)
y <- rnorm(T, true.mu, sigma.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, sigma1=0, sigma2=0, tau=0))
pos.mu <- grep("^mu$", parm.names)
pos.sigma1 <- grep("sigma1", parm.names)
pos.sigma2 <- grep("sigma2", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    mu <- rnorm(1, mean(Data$y), sd(Data$y))
    sigma1 <- rhalfcauchy(1, 5)
    sigma2 <- rhalfcauchy(1, 5)
    tau <- runif(1, 1, Data$T)
    return(c(mu, sigma1, sigma2, tau))
    }
Data <- list(T=T, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.mu=pos.mu, pos.sigma1=pos.sigma1, pos.sigma2=pos.sigma2,
    pos.tau=pos.tau, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    sigma1 <- interval(parm[Data$pos.sigma1], 1e-100, Inf)
    parm[Data$pos.sigma1] <- sigma1
    sigma2 <- interval(parm[Data$pos.sigma2], 1e-100, Inf)
    parm[Data$pos.sigma2] <- sigma2
    tau <- interval(parm[Data$pos.tau], 1, Data$T)
    parm[Data$pos.tau] <- tau
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 1000, log=TRUE)
    sigma1.prior <- dhalfcauchy(sigma1, 25, log=TRUE)
    sigma2.prior <- dhalfcauchy(sigma2, 25, log=TRUE)
    tau.prior <- dunif(tau, 1, Data$T, log=TRUE)
    ### Log-Likelihood
    sigma.vec <- ifelse(1:Data$T < floor(tau), sigma1, sigma2)
    LL <- sum(dnorm(Data$y, mu, sigma.vec, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + sigma1.prior + sigma2.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$T, mu, sigma.vec), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(mean(y), 1, 1, T/2)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery
cat("mu -- true:", true.mu,
    " post. mean:", round(fit$Summary2[pos.mu, "Mean"], 3), "\n")
cat("sigma1 -- true:", true.sigma1,
    " post. mean:", round(fit$Summary2[pos.sigma1, "Mean"], 3), "\n")
cat("sigma2 -- true:", true.sigma2,
    " post. mean:", round(fit$Summary2[pos.sigma2, "Mean"], 3), "\n")
cat("tau -- true:", true.tau,
    " post. mean:", round(fit$Summary2[pos.tau, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]     <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma1] <- true.sigma1[seq_along(pos.sigma1)]
ground_truth[pos.sigma2] <- true.sigma2[seq_along(pos.sigma2)]
ground_truth[pos.tau]    <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Change Point Regression

Change point regression fits a piecewise linear model with two slopes
joined at an unknown breakpoint \\\theta\\. Below \\\theta\\ the
response follows one linear trend; above it the slope changes by
\\\beta_2\\. Both the regression coefficients and the change point
location are estimated jointly, making this a useful model for detecting
structural breaks in continuous data. Bayesian change point analysis was
developed by Carlin et al. [\[19\]](#ref19) and extended by Chib
[\[24\]](#ref24) using Markov chain methods. Change point models detect
regime shifts in climate records, such as abrupt transitions in Arctic
sea ice extent.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \alpha + \beta_1
\textbf{x} + \beta_2 (\textbf{x} - \theta)\[(\textbf{x} - \theta) \>
0\]\\ \\\alpha \sim \mathcal{N}(0, 1000)\\ \\\beta \sim \mathcal{N}(0,
1000)\\ \\\sigma \sim \mathcal{HC}(25)\\ \\\theta \sim \mathcal{U}(-1.3,
1.1)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the mean function involves a conditional indicator
\\(\textbf{x} - \theta) \cdot \mathbb{1}\[(\textbf{x} - \theta) \> 0\]\\
that depends on the estimated change point parameter \\\theta\\. This
kind of parameter-dependent piecewise logic falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42118)
N <- 60

#### True parameter values
true.alpha <- 1.0
true.beta <- c(-0.5, -1.2)
true.sigma <- 0.15
true.theta <- 0.3

#### Generate covariate and response
x <- sort(runif(N, -1.3, 1.1))
mu.true <- true.alpha + true.beta[1] * x +
    true.beta[2] * (x - true.theta) * ((x - true.theta) > 0)
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, beta=rep(0,2), sigma=0, theta=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    beta <- rnorm(2)
    sigma <- runif(1)
    theta <- runif(1)
    return(c(alpha, beta, sigma, theta))
    }
Data <- list(N=N, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.beta=pos.beta, pos.sigma=pos.sigma,
    pos.theta=pos.theta, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    theta <- interval(parm[Data$pos.theta], -1.3, 1.1)
    parm[Data$pos.theta] <- theta
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    theta.prior <- dunif(theta, -1.3, 1.1, log=TRUE)
    ### Log-Likelihood
    mu <- alpha + beta[1]*Data$x + beta[2]*(Data$x - theta)*{(Data$x - theta) > 0}
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + sigma.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0.2, -0.45, 0, 0.2, 0)
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery
cat("alpha -- true:", true.alpha,
    " post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")
for (j in 1:2)
    cat("beta[", j, "] -- true:", true.beta[j],
        " post. mean:", round(fit$Summary2[pos.beta[j], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
cat("theta -- true:", true.theta,
    " post. mean:", round(fit$Summary2[pos.theta, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Dynamic Sparse Factor Model (DSFM)

The dynamic sparse factor model extends standard dynamic factor analysis
by allowing both the factor loadings and the observation-error variances
to evolve over time through AR(1) processes, with a sparsity threshold
\\\lambda^d\\ that sets small loadings to zero. This produces a
parsimonious time-varying decomposition suitable for high-dimensional
financial or macroeconomic panel data where the number of active factors
and their influence on individual series may shift over the sample
period. Dynamic factor models extend the state space framework of West
and Harrison [\[112\]](#ref112), with sparse factor loadings following
the shrinkage approach of Bhattacharya and Dunson [\[12\]](#ref12).
Dynamic sparse factor models are applied in macroeconomic nowcasting to
extract a small number of latent factors from hundreds of economic
indicators.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{N}(\alpha{t,j} + \textbf{F}\_{t,1:P}
\Lambda\_{1:P,1:j,t}, \Sigma^2\_{t,j}), \quad t=1,\dots,T, \quad
j=1,\dots,J\\ \\\alpha\_{t,j} \sim \mathcal{N}(\alpha^\mu_j +
\alpha^\phi_j(\alpha\_{t-1,j} - \alpha^mu_j), \alpha^\sigma2_j)\\
\\\textbf{F}\_{t,1:P} \sim \mathcal{N}\_P(\textbf{F}^\phi
\textbf{F}\_{t-1,1:P}, \textbf{f}^\Sigma\_{t,1:P})\\
\\\textbf{f}^\Sigma\_{t,1:P} =
t(\textbf{f}^\textbf{U}\_{1:P,1:P,t})\textbf{f}^\textbf{U}\_{1:P,1:P,t}\\
\\\textbf{f}^\textbf{U}\_{p,q,t} \sim
\mathcal{N}(\textbf{f}^{\textbf{u}\_\mu}\_{p,q} +
\textbf{f}^{\textbf{u}\_\phi}\_{p,q}(\textbf{f}^\textbf{U}\_{p,q,t-1} -
\textbf{f}^{\textbf{u}\_\mu}\_{p,q}),
\textbf{f}^{\textbf{u}\_\sigma^2}\_{p,q})\\ \\\Lambda\_{p,j,t} \sim
\mathcal{N}(\lambda^\mu\_{p,j} +
\lambda^\phi\_{p,j}(\Lambda\_{p,j,t-1} - \lambda^mu\_{p,j}),
\lambda^\sigma2\_{p,j})\\ \\\Sigma\_{t,j} = \exp(\log(\Sigma\_{t,j}))\\
\\log(\Sigma\_{t,j}) \sim \mathcal{N}(\sigma^\mu_j +
\sigma^\phi_j(log(\Sigma\_{t-1,j}) - \sigma^mu_j), \sigma^\sigma2_j)\\
\\\alpha^0_j \sim \mathcal{N}(0, 1), \quad j=1,\dots,J\\ \\\alpha^\mu_j
\sim \mathcal{N}(0, 1), \quad j=1,\dots,J\\ \\\frac{\alpha^\phi_j+1}{2}
\sim \mathcal{BETA}(20, 1.5), \quad j=1,\dots,J\\ \\\alpha^\sigma_j \sim
\mathcal{HC}(5), \quad j=1,\dots,J\\ \\\textbf{f}^0_j \sim
\mathcal{N}(0, 1), \quad j=1,\dots,J\\ \\\frac{\textbf{f}^\phi_j+1}{2}
\sim \mathcal{BETA}(1, 1), \quad j=1,\dots,J\\
\\\textbf{f}^{\textbf{u}0}\_j \sim \mathcal{N}(0, 1), \quad
j=1,\dots,J\\ \\\textbf{f}^{\textbf{u}\mu}\_j \sim \mathcal{N}(0, 1),
\quad j=1,\dots,J\\ \\\frac{\textbf{f}^{\textbf{u}\phi}\_j+1}{2} \sim
\mathcal{BETA}(20, 1.5), \quad j=1,\dots,J\\
\\\textbf{f}^{\textbf{u}\sigma}\_j \sim \mathcal{HC}(1), \quad
j=1,\dots,J\\ \\\lambda^0_j \sim \mathcal{N}(0, 1), \quad j=1,\dots,J\\
\\\lambda^d_j \sim \mathcal{U}(0, \|\lambda^\mu_j\| +
3\sqrt{\frac{\lambda^\sigma_j}{1 - \lambda^\phi_j\lambda^\phi_j}}),
\quad j=1,\dots,J\\ \\\lambda^\mu_j \sim \mathcal{N}(0, 1), \quad
j=1,\dots,J\\ \\\frac{\lambda^\phi_j+1}{2} \sim \mathcal{BETA}(20, 1.5),
\quad j=1,\dots,J\\ \\\lambda^\sigma_j \sim \mathcal{HC}(1), \quad
j=1,\dots,J\\ \\\log(\sigma^0_j) \sim \mathcal{N}(0, 1), \quad
j=1,\dots,J\\ \\\log(\sigma^\mu_j) \sim \mathcal{N}(0, 1), \quad
j=1,\dots,J\\ \\\frac{\log(\sigma^\phi_j)+1}{2} \sim \mathcal{BETA}(20,
1.5), \quad j=1,\dots,J\\ \\\log(\sigma^\sigma_j) \sim \mathcal{HC}(1),
\quad j=1,\dots,J\\

#### model_spec() notation

The DSFM model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the time-varying factor loadings, dynamic observation-error
variances, stochastic volatility in the factor covariance, and the
sparsity threshold mechanism all require imperative loops and dynamic
array construction that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42185)

#### Simulation dimensions
T <- 20  # Number of time-periods
J <- 3   # Number of time-series
P <- 2   # Number of dynamic factors

#### Generate simple time-varying factor data
Y <- matrix(rnorm(T * J), T, J)
F.true <- matrix(rnorm(T * P), T, P)
Lambda.true <- matrix(c(1, 0.5, 0.5, 1, 0.3, 0.7), P, J)
for (t in 1:T) {
    Y[t, ] <- F.true[t, ] %*% Lambda.true + rnorm(J, 0, 0.5)
}

#### Center and scale
Y.means <- colMeans(Y)
Y <- Y - matrix(Y.means, T, J, byrow=TRUE)
Y.scales <- sqrt(apply(Y, 2, var))
Y <- Y / matrix(Y.scales, T, J, byrow=TRUE)
mon.names <- "LP"
U1 <- matrix(NA,P,P); U2 <- matrix(NA,P,J)
U1[upper.tri(U1, diag=TRUE)] <- 0; U2[upper.tri(U2)] <- 0
Lambda <- array(NA, dim=c(P,J,T))
U <- array(NA, dim=c(P,P,T))
for (t in 1:T) {
    U[ , , t] <- U1
    Lambda[ , , t] <- U2}
parm.names <- as.parm.names(list(alpha0=rep(0,J), Alpha=matrix(0,T,J),
    alpha.mu=rep(0,J), alpha.phi=rep(0,J), alpha.sigma=rep(0,J),
    f0=rep(0,P), F=matrix(0,T,P), f.phi=rep(0,P), f.u0=U1, f.U=U,
    f.u.mu=U1, f.u.phi=U1, f.u.sigma=U1, lambda0=U2, Lambda=Lambda,
    lambda.d=U2, lambda.mu=U2, lambda.phi=U2, lambda.sigma=U2,
    lsigma0=rep(0,J), lSigma=matrix(0,T,J),
    lsigma.mu=rep(0,J), lsigma.phi=rep(0,J), lsigma.sigma=rep(0,J)))
pos.alpha0 <- grep("alpha0", parm.names)
pos.Alpha <- grep("Alpha", parm.names)
pos.alpha.mu <- grep("alpha.mu", parm.names)
pos.alpha.phi <- grep("alpha.phi", parm.names)
pos.alpha.sigma <- grep("alpha.sigma", parm.names)
pos.f0 <- grep("f0", parm.names)
pos.F <- grep("F", parm.names)
pos.f.phi <- grep("f.phi", parm.names)
pos.f.u0 <- grep("f.u0", parm.names)
pos.f.U <- grep("f.U", parm.names)
pos.f.u.mu <- grep("f.u.mu", parm.names)
pos.f.u.phi <- grep("f.u.phi", parm.names)
pos.f.u.sigma <- grep("f.u.sigma", parm.names)
pos.lambda0 <- grep("lambda0", parm.names)
pos.Lambda <- grep("Lambda", parm.names)
pos.lambda.d <- grep("lambda.d", parm.names)
pos.lambda.mu <- grep("lambda.mu", parm.names)
pos.lambda.phi <- grep("lambda.phi", parm.names)
pos.lambda.sigma <- grep("lambda.sigma", parm.names)
pos.lsigma0 <- grep("lsigma0", parm.names)
pos.lSigma <- grep("lSigma", parm.names)
pos.lsigma.mu <- grep("lsigma.mu", parm.names)
pos.lsigma.phi <- grep("lsigma.phi", parm.names)
pos.lsigma.sigma <- grep("lsigma.sigma", parm.names)
PGF <- function(Data) {
    alpha0 <- rnorm(Data$J)
    Alpha <- rnorm(Data$T*Data$J)
    alpha.mu <- rnorm(Data$J)
    alpha.phi <- rbeta(Data$J, 20, 1.5) * 2 - 1
    alpha.sigma <- runif(Data$J)
    f0 <- rnorm(Data$P)
    F <- rnorm(Data$T*Data$P)
    f.phi <- rbeta(Data$P, 1, 1) * 2 - 1
    f.u0 <- rnorm(length(Data$pos.f.u0))
    f.U <- rnorm(length(Data$pos.f.U))
    f.u.mu <- rnorm(length(Data$pos.f.u.mu))
    f.u.phi <- runif(length(Data$pos.f.u.phi))
    f.u.sigma <- runif(length(Data$pos.f.u.sigma))
    lambda0 <- rnorm(length(Data$pos.lambda0))
    Lambda <- rnorm(length(Data$pos.Lambda))
    lambda.mu <- rnorm(length(Data$pos.lambda.mu))
    lambda.phi <- rbeta(length(Data$pos.lambda.phi), 20, 1.5)
    lambda.sigma <- runif(length(Data$pos.lambda.sigma))
    lambda.d <- runif(length(Data$pos.lambda.d), 0, abs(lambda.mu) +
        3*sqrt(lambda.sigma/(1-lambda.phi^ 2)))
    lsigma0 <- rnorm(Data$J)
    lSigma <- rnorm(Data$T*Data$J)
    lsigma.mu <- rnorm(Data$J)
    lsigma.phi <- rbeta(Data$J, 20, 1.5) * 2 - 1
    lsigma.sigma <- runif(Data$J)
    return(c(alpha0, Alpha, alpha.mu, alpha.phi, alpha.sigma, f0, F,
        f.phi, f.u0, f.U, f.u.mu, f.u.phi, f.u.sigma, lambda0, Lambda,
        lambda.d, lambda.mu, lambda.phi, lambda.sigma, lsigma0, lSigma,
        lsigma.mu, lsigma.phi, lsigma.sigma))
    }
Data <- list(J=J, P=P, PGF=PGF, T=T, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.alpha0=pos.alpha0, pos.Alpha=pos.Alpha,
    pos.alpha.mu=pos.alpha.mu, pos.alpha.phi=pos.alpha.phi,
    pos.alpha.sigma=pos.alpha.sigma, pos.f0=pos.f0, pos.F=pos.F,
    pos.f.phi=pos.f.phi, pos.f.u0=pos.f.u0, pos.f.U=pos.f.U,
    pos.f.u.mu=pos.f.u.mu, pos.f.u.phi=pos.f.u.phi,
    pos.f.u.sigma=pos.f.u.sigma, pos.lambda0=pos.lambda0,
    pos.Lambda=pos.Lambda, pos.lambda.d=pos.lambda.d,
    pos.lambda.mu=pos.lambda.mu, pos.lambda.phi=pos.lambda.phi,
    pos.lambda.sigma=pos.lambda.sigma, pos.lsigma0=pos.lsigma0,
    pos.lSigma=pos.lSigma, pos.lsigma.mu=pos.lsigma.mu,
    pos.lsigma.phi=pos.lsigma.phi, pos.lsigma.sigma=pos.lsigma.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha0 <- parm[Data$pos.alpha0]
    Alpha <- matrix(parm[Data$pos.Alpha], Data$T, Data$J)
    alpha.mu <- parm[Data$pos.alpha.mu]
    alpha.phi <- interval(parm[Data$pos.alpha.phi], -1, 1)
    parm[Data$pos.alpha.phi] <- alpha.phi
    alpha.sigma <- interval(parm[Data$pos.alpha.sigma], 1e-100, Inf)
    parm[Data$pos.alpha.sigma] <- alpha.sigma
    f0 <- parm[Data$pos.f0]
    F <- matrix(parm[Data$pos.F], Data$T, Data$P)
    f.phi <- interval(parm[Data$pos.f.phi], -1, 1)
    parm[Data$pos.f.phi] <- f.phi
    f.u0 <- parm[Data$pos.f.u0]
    f.U <- parm[Data$pos.f.U]
    f.u.mu <- parm[Data$pos.f.u.mu]
    f.u.phi <- interval(parm[Data$pos.f.u.phi], -1, 1)
    parm[Data$pos.f.u.phi] <- f.u.phi
    f.u.sigma <- interval(parm[Data$pos.f.u.sigma], 1e-100, Inf)
    parm[Data$pos.f.u.sigma] <- f.u.sigma
    lambda0 <- parm[Data$pos.lambda0]
    Lambda <- parm[Data$pos.Lambda]
    lambda.mu <- parm[Data$pos.lambda.mu]
    lambda.phi <- interval(parm[Data$pos.lambda.phi], -1, 1)
    parm[Data$pos.lambda.phi] <- lambda.phi
    lambda.sigma <- interval(parm[Data$pos.lambda.sigma], 1e-100, Inf)
    parm[Data$pos.lambda.sigma] <- lambda.sigma
    lambda.d <- parm[Data$pos.lambda.d]
    for (i in 1:length(lambda.d))
        lambda.d[i] <- interval(lambda.d[i], 0, abs(lambda.mu[i]) +
        3*sqrt(lambda.sigma[i]/(1-lambda.phi[i]^ 2)))
    parm[Data$pos.lambda.d] <- lambda.d
    lsigma0 <- parm[Data$pos.lsigma0]
    lSigma <- matrix(parm[Data$pos.lSigma], Data$T, Data$J)
    lsigma.mu <- parm[Data$pos.lsigma.mu]
    lsigma.phi <- interval(parm[Data$pos.lsigma.phi], -1, 1)
    parm[Data$pos.lsigma.phi] <- lsigma.phi
    lsigma.sigma <- interval(parm[Data$pos.lsigma.sigma], 1e-100, Inf)
    parm[Data$pos.lsigma.sigma] <- lsigma.sigma
    ### Log-Prior
    alpha0.prior <- sum(dnorm(alpha0, 0, 1, log=TRUE))
    Alpha.prior <- sum(dnorm(Alpha,
        matrix(alpha.mu, Data$T, Data$J, byrow=TRUE) +
        matrix(alpha.phi, Data$T, Data$J, byrow=TRUE) *
        (rbind(alpha0, Alpha[-Data$T,]) -
        matrix(alpha.mu, Data$T, Data$J, byrow=TRUE)),
        matrix(alpha.sigma, Data$T, Data$J, byrow=TRUE), log=TRUE))
    alpha.mu.prior <- sum(dnorm(alpha.mu, 0, 1, log=TRUE))
    alpha.phi.prior <- sum(dbeta((alpha.phi + 1) / 2, 20, 1.5, log=TRUE))
    alpha.sigma.prior <- sum(dhalfcauchy(alpha.sigma, 5, log=TRUE))
    f0.prior <- sum(dnorm(f0, 0, 1, log=TRUE))
    f.phi.prior <- sum(dbeta((f.phi + 1) / 2, 1, 1, log=TRUE))
    f.u0.prior <- sum(dnorm(f.u0, 0, 1, log=TRUE))
    f.U.prior <- sum(dnorm(matrix(f.U, nrow=Data$T, byrow=TRUE),
        matrix(f.u.mu, Data$T, Data$P*(Data$P-1)/2+Data$P, byrow=TRUE) +
        matrix(f.u.phi, Data$T, Data$P*(Data$P-1)/2+Data$P, byrow=TRUE) *
        (rbind(f.u0, matrix(f.U, nrow=Data$T, byrow=TRUE)[-Data$T,]) -
        matrix(f.u.mu, Data$T, Data$P*(Data$P-1)/2+Data$P, byrow=TRUE)),
        matrix(f.u.sigma, Data$T, Data$P*(Data$P-1)/2+Data$P, byrow=TRUE),
        log=TRUE))
    f.u.mu.prior <- sum(dnorm(f.u.mu, 0, 1, log=TRUE))
    f.u.phi.prior <- sum(dbeta((f.u.phi + 1) / 2, 20, 1.5, log=TRUE))
    f.u.sigma.prior <- sum(dhalfcauchy(f.u.sigma, 1, log=TRUE))
    lambda0.prior <- sum(dnorm(lambda0, 0, 1, log=TRUE))
    Lambda.prior <- sum(dnorm(matrix(Lambda, nrow=Data$T, byrow=TRUE),
        matrix(lambda.mu, Data$T, length(lambda.mu), byrow=TRUE) +
        (rbind(lambda0, matrix(Lambda, nrow=Data$T, byrow=TRUE))[-(Data$T+1),] -
        matrix(lambda.mu, Data$T, length(lambda.mu), byrow=TRUE)),
        matrix(lambda.sigma, Data$T, length(lambda.sigma), byrow=TRUE),
        log=TRUE))
    lambda.d.prior <- sum(dunif(lambda.d, 0, abs(lambda.mu) +
        3*sqrt(lambda.sigma/(1-lambda.phi^ 2)), log=TRUE))
    lambda.mu.prior <- sum(dnorm(lambda.mu, 0, 1, log=TRUE))
    lambda.phi.prior <- sum(dbeta((lambda.phi + 1) / 2, 20, 1.5, log=TRUE))
    lambda.sigma.prior <- sum(dhalfcauchy(lambda.sigma, 1, log=TRUE))
    lsigma0.prior <- sum(dnorm(lsigma0, 0, 1, log=TRUE))
    lSigma.prior <- sum(dnorm(lSigma,
        matrix(lsigma.mu, Data$T, Data$J, byrow=TRUE) +
        matrix(lsigma.phi, Data$T, Data$J, byrow=TRUE) *
        (rbind(lsigma0, lSigma[-Data$T,]) -
        matrix(lsigma.mu, Data$T, Data$J, byrow=TRUE)),
        matrix(lsigma.sigma, Data$T, Data$J, byrow=TRUE), log=TRUE))
    lsigma.mu.prior <- sum(dnorm(lsigma.mu, 0, 1, log=TRUE))
    lsigma.phi.prior <- sum(dbeta((lsigma.phi + 1) / 2, 20, 1.5, log=TRUE))
    lsigma.sigma.prior <- sum(dhalfcauchy(lsigma.sigma, 1, log=TRUE))
    ### Log-Likelihood
    LL <- 0; Yhat <- Data$Y; F.prior <- 0
    for (t in 1:Data$T) {
        f.U.temp <- matrix(0, Data$P, Data$P)
        f.U.temp[upper.tri(f.U.temp, diag=TRUE)] <- matrix(f.U, nrow=Data$T,
    byrow=TRUE)[t,]
        diag(f.U.temp) <- exp(diag(f.U.temp))
        f.Sigma <- as.symmetric.matrix(t(f.U.temp) %*% f.U.temp)
        F.prior <- F.prior + dmvn(F[t,], rbind(f0, F)[t,] %*% diag(f.phi),
    f.Sigma, log=TRUE)
        Lambda.temp <- matrix(1, Data$P, Data$J)
        Lambda.temp[lower.tri(Lambda.temp)] <- 0
        Lambda.temp[upper.tri(Lambda.temp)] <- matrix(Lambda,
    nrow=Data$T, byrow=TRUE)[t,]*(abs(matrix(Lambda,
    nrow=Data$T, byrow=TRUE)[t,]) > lambda.d)
        mu <- Alpha[t,] + F[t,] %*% Lambda.temp
        LL <- LL + sum(dnorm(Data$Y[t,], mu, exp(lSigma[t,]), log=TRUE))
        Yhat[t,] <- rnorm(Data$J, mu, exp(lSigma[t,])) #Fitted
        }
    ### Log-Posterior
    LP <- LL + alpha0.prior + Alpha.prior + alpha.mu.prior +
        alpha.phi.prior + alpha.sigma.prior + f0.prior + F.prior +
        f.phi.prior + f.u0.prior + f.U.prior + f.u.mu.prior +
        f.u.phi.prior + f.u.sigma.prior + lambda0.prior +
        Lambda.prior + lambda.d.prior + lambda.mu.prior +
        lambda.phi.prior + lambda.sigma.prior + lsigma0.prior +
        lSigma.prior + lsigma.mu.prior + lsigma.phi.prior +
        lsigma.sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=Yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rnorm(J), rnorm(T*J), rnorm(J), runif(J), runif(J),
    rnorm(P), rnorm(T*P), rbeta(P,1,1)*2-1, rnorm(P*(P-1)/2+P),
    rnorm((P*(P-1)/2+P)*T), rnorm(P*(P-1)/2+P),
    rbeta(P*(P-1)/2+P,1,1)*2-1, runif(P*(P-1)/2+P),
    rnorm(P*J-P-P*(P-1)/2), rnorm((P*J-P-P*(P-1)/2)*T),
    runif(P*J-P-P*(P-1)/2,0,1e-3), rnorm(P*J-P-P*(P-1)/2),
    rbeta(P*J-P-P*(P-1)/2,20,1.5)*2-1, runif(P*J-P-P*(P-1)/2),
    rnorm(J), rnorm(T*J), rnorm(J), rbeta(J,20,1.5)*2-1, runif(J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Hidden Markov Model

Hidden Markov models were developed by Baum et al. [\[5\]](#ref5), with
Bayesian inference following Scott [\[98\]](#ref98) and
Frühwirth-Schnatter [\[42\]](#ref42). HMMs are fundamental in speech
recognition and are applied in ecology to classify animal behavioural
states (foraging, resting, travelling) from GPS tracking data. \### Form

This introductory hidden Markov model (HMM) includes \\N\\ discrete
states. \\\textbf{y}\_t \sim \mathcal{N}(\mu\_\theta, \sigma^2\_\theta),
\quad t=1,\dots,T\\ \\\mu \sim \mathcal{N}(\mu_0, \sigma^2)\\ \\\sigma^2
\sim \mathcal{HC}(25)\\ \\\theta_t \sim
\mathcal{CAT}(\phi\_{\theta\_{t-1},1:N}), \quad t=1,\dots,T\\
\\\phi\_{i,1:N} \sim \mathcal{D}(\alpha\_{1:N}), \quad i=1,\dots,N\\
\\\mu_0 \sim \mathcal{N}(0, 1000)\\ \\\sigma^2_0 \sim \mathcal{HC}(25)\\

#### Data

``` r

data(demonfx)
y <- as.vector(log(as.matrix(demonfx[1:50,1])))
T <- length(y) #Number of time-periods
N <- 2 #Number of discrete (hidden) states
alpha <- matrix(1,N,N) #Concentration hyperparameter
mon.names <- "LP"
parm.names <- as.parm.names(list(mu0=rep(0,N), mu1=rep(0,N),
    phi=matrix(0,N,N), sigma2=rep(0,N), theta=rep(0,T)))
pos.mu0 <- grep("mu0", parm.names)
pos.mu1 <- grep("mu1", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma2 <- grep("sigma2", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    mu0 <- sort(runif(Data$N, min(Data$y), max(Data$y)))
    mu1 <- sort(runif(Data$N, min(Data$y), max(Data$y)))
    phi <- matrix(runif(Data$N*Data$N), Data$N, Data$N)
    phi <- as.vector(phi / rowSums(phi))
    sigma2 <- runif(Data$N)
    theta <- rcat(Data$T, rep(1/Data$N,Data$N))
    return(c(mu0, mu1, phi, sigma2, theta))
    }
Data <- list(N=N, PGF=PGF, T=T, alpha=alpha, mon.names=mon.names,
    parm.names=parm.names, pos.mu0=pos.mu0, pos.mu1=pos.mu1,
    pos.phi=pos.phi, pos.sigma2=pos.sigma2, pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu0 <- interval(parm[Data$pos.mu0], min(Data$y), max(Data$y))
    parm[Data$pos.mu0] <- mu0
    mu <- interval(parm[Data$pos.mu1], min(Data$y), max(Data$y))
    parm[Data$pos.mu1] <- mu <- sort(mu)
    phi <- matrix(abs(parm[Data$pos.phi]), Data$N, Data$N)
    parm[Data$pos.phi] <- phi <- phi / rowSums(phi)
    sigma2 <- interval(parm[Data$pos.sigma2], 1e-100, Inf)
    parm[Data$pos.sigma2] <- sigma2
    theta <- parm[Data$pos.theta]
    ### Log-Hyperprior
    mu0.prior <- sum(dnormv(mu0, 0, 1000, log=TRUE))
    ### Log-Prior
    mu.prior <- sum(dnormv(mu, mu0, sigma2, log=TRUE))
    phi.prior <- 0
    for (i in 1:Data$N)
        phi.prior <- phi.prior + sum(ddirichlet(phi[i,], Data$alpha[i,],
    log=TRUE))
    sigma2.prior <- sum(dhalfcauchy(sigma2, 25, log=TRUE))
    theta.prior <- sum(dcat(theta, rbind(rep(1/Data$N,Data$N),
        phi[theta[-Data$T],]), log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnormv(Data$y, mu[theta], sigma2[theta], log=TRUE))
    ### Log-Posterior
    LP <- LL + mu0.prior + mu.prior + phi.prior + sigma2.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormv(length(theta), mu[theta], sigma2[theta]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(sort(runif(N, min(y), max(y))),
    sort(runif(N, min(y), max(y))), runif(N*N), runif(N),
    rcat(T, rep(1/N,N)))
```

#### model_spec() notation

The HMM cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the hidden state sequence \\\theta_t\\ is a discrete latent
variable drawn from a transition matrix \\\phi\\ that depends on the
previous state, and the emission parameters are indexed by these latent
states. This sequential discrete-state switching structure with
transition-dependent priors falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Inspect estimated state means and transition matrix
mu.post <- fit$Summary2[pos.mu1, "Mean"]
cat("Estimated state means (mu):", round(mu.post, 3), "\n")
phi.post <- matrix(fit$Summary2[pos.phi, "Mean"], N, N)
cat("Estimated transition matrix (phi):\n")
print(round(phi.post, 3))

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Interrupted Time Series **\[NEW\]**

Interrupted time series (ITS) analysis is a quasi-experimental design
for evaluating the causal effect of an intervention that occurs at a
known time point \\T^\*\\. The model decomposes the conditional mean
into a pre-intervention linear trend, an immediate level shift at the
intervention, and a change in slope after the intervention, following
the segmented regression framework described by Wagner et
al. [\[134\]](#ref134) and reviewed by Bernal et al. [\[133\]](#ref133).
Four regression coefficients capture the baseline level \\\beta_0\\, the
pre-intervention trend \\\beta_1\\, the level change \\\beta_2\\, and
the slope change \\\beta_3\\, while a half-Cauchy prior constrains the
residual standard deviation to be positive. ITS designs are the standard
tool in public health for evaluating the impact of policy changes on
medication use, hospitalisation rates, and mortality.

#### Form

\\y_t \sim \mathcal{N}(\mu_t,\\ \sigma^2), \quad t=1,\dots,T\\ \\\mu_t =
\beta_0 + \beta_1\\ t + \beta_2\\ D_t + \beta_3\\ (t - T^\*)\\ D_t\\
\\D_t = \mathbb{1}(t \ge T^\*)\\ \\\beta_j \sim \mathcal{N}(0, 1000),
\quad j=0,\dots,3\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

With the design matrix precomputed, this is a standard linear regression
that
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
handles directly.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = beta0 + beta1 * t_scaled + beta2 * D + beta3 * tD
  beta0 ~ NormalV(0, 1000)
  beta1 ~ NormalV(0, 1000)
  beta2 ~ NormalV(0, 1000)
  beta3 ~ NormalV(0, 1000)
  sigma ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42218)
TT <- 200
T.star <- 100

#### True parameter values
true.beta <- c(10, 0.5, -3, 0.8)
true.sigma <- 1

#### Generate time series with intervention
t.seq <- 1:TT
D <- as.numeric(t.seq >= T.star)
tD <- (t.seq - T.star) * D
X <- cbind(1, t.seq, D, tD)
mu.true <- as.vector(X %*% true.beta)
y <- rnorm(TT, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0, 4), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(4)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(TT=TT, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(tcrossprod(Data$X, t(beta)))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$TT, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, 4), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("beta0 -- true:", true.beta[1],
    " post. mean:", round(fit$Summary2[pos.beta[1], "Mean"], 3), "\n")
cat("beta1 -- true:", true.beta[2],
    " post. mean:", round(fit$Summary2[pos.beta[2], "Mean"], 3), "\n")
cat("beta2 -- true:", true.beta[3],
    " post. mean:", round(fit$Summary2[pos.beta[3], "Mean"], 3), "\n")
cat("beta3 -- true:", true.beta[4],
    " post. mean:", round(fit$Summary2[pos.beta[4], "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Ornstein-Uhlenbeck Process **\[NEW\]**

The Ornstein-Uhlenbeck process is a continuous-time mean-reverting
stochastic process whose transition density is Gaussian with conditional
mean and variance that depend on the elapsed time between observations.
The parameter \\\theta\\ governs the speed of reversion toward the
long-run level \\\mu\\, while \\\sigma\\ sets the diffusion intensity.
At regular time intervals \\\Delta t\\, the exact conditional
distribution is available in closed form, enabling straightforward
likelihood-based inference without discretisation error. The process was
introduced by Uhlenbeck and Ornstein [\[147\]](#ref147) in the context
of Brownian motion with friction, and Vasicek [\[148\]](#ref148) applied
it to short-rate interest rate dynamics. Ornstein-Uhlenbeck models are
used in quantitative finance for mean-reverting interest rate processes
and in ecology for modelling animal home-range behaviour.

#### Form

\\y_t \mid y\_{t-1} \sim \mathcal{N}(m_t, v)\\ \\m_t = \mu + (y\_{t-1} -
\mu)\exp(-\theta \Delta t)\\ \\v = \frac{\sigma^2}{2\theta}\left(1 -
\exp(-2\theta \Delta t)\right)\\ \\\mu \sim \mathcal{N}(0, 1000)\\
\\\theta \sim \mathcal{HC}(25)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the conditional distribution depends nonlinearly on the previous
observation through the exponential mean-reversion kernel \\\exp(-\theta
\Delta t)\\, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42224)
TT <- 300
dt <- 0.1
true.mu <- 5
true.theta <- 2
true.sigma <- 1

#### Simulate Ornstein-Uhlenbeck trajectory
y <- numeric(TT)
y[1] <- true.mu
for (t in 2:TT) {
    cond.mean <- true.mu + (y[t-1] - true.mu) * exp(-true.theta * dt)
    cond.var <- true.sigma^2 / (2 * true.theta) *
        (1 - exp(-2 * true.theta * dt))
    y[t] <- rnorm(1, cond.mean, sqrt(cond.var))
    }

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, theta=0, sigma=0))
pos.mu <- grep("mu", parm.names)
pos.theta <- grep("theta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    return(c(rnorm(1), runif(1, 0.1, 5), runif(1)))
    }
Data <- list(PGF=PGF, TT=TT, dt=dt, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.sigma=pos.sigma,
    pos.theta=pos.theta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    theta <- interval(parm[Data$pos.theta], 1e-100, Inf)
    parm[Data$pos.theta] <- theta
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 1000, log=TRUE)
    theta.prior <- dhalfcauchy(theta, 25, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    decay <- exp(-theta * Data$dt)
    cond.mean <- mu + (Data$y[1:(Data$TT - 1)] - mu) * decay
    cond.var <- sigma^2 / (2 * theta) * (1 - exp(-2 * theta * Data$dt))
    LL <- sum(dnorm(Data$y[2:Data$TT], cond.mean, sqrt(cond.var),
        log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + theta.prior + sigma.prior
    yhat <- c(Data$y[1], rnorm(Data$TT - 1, cond.mean, sqrt(cond.var)))
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(mean(y), 1, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Random Walk (1st/2nd order) **\[NEW\]**

A random walk of order one defines the first differences \\\Delta_t =
y_t - y\_{t-1}\\ as independent normal draws with drift \\\mu\\ and
innovation variance \\\sigma^2\\, so that after differencing the model
reduces to an iid Gaussian likelihood on the increments. A second-order
random walk applies the same logic to the second differences
\\\Delta^2_t = y_t - 2 y\_{t-1} + y\_{t-2}\\, penalising changes in
slope rather than level. The classical theory of random walks originates
with Pearson [\[125\]](#ref125), who coined the term in the context of a
drunkard’s trajectory; the first-order case is implemented here, with
the second-order variant described for reference. Random walk models
underpin the analysis of financial asset returns and are used in ecology
to characterise population fluctuations driven by stochastic
environmental forcing.

#### Form

\\\Delta_t = y_t - y\_{t-1} \sim \mathcal{N}(\mu, \sigma^2), \quad
t=2,\dots,T\\ \\\mu \sim \mathcal{N}(0, 1000)\\ \\\sigma \sim
\mathcal{HC}(25)\\

For RW(2), the second differences \\\Delta^2_t = y_t - 2 y\_{t-1} +
y\_{t-2}\\ replace \\\Delta_t\\ and the likelihood operates on \\t = 3,
\dots, T\\.

#### model_spec() notation

The first-differenced random walk reduces to an iid normal model on the
increments, which is directly expressible in
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

``` r

spec <- model_spec("
  delta ~ Normal(mu, sigma)
  mu ~ NormalV(0, 1000)
  sigma ~ HalfCauchy(25)
")
```

#### Ground truth and data

``` r

set.seed(42228)
T.obs <- 300

#### True parameter values
true.mu <- 0.1
true.sigma <- 1.0

#### Simulate random walk
increments <- rnorm(T.obs - 1, true.mu, true.sigma)
y.orig <- cumsum(c(0, increments))
delta <- diff(y.orig)
T_diff <- length(delta)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, sigma=0))
pos.mu <- grep("mu", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    mu <- rnorm(1)
    sigma <- runif(1)
    return(c(mu, sigma))
    }
Data <- list(T_diff=T_diff, PGF=PGF, delta=delta, y.orig=y.orig,
    mon.names=mon.names, parm.names=parm.names,
    pos.mu=pos.mu, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 1000, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    LL <- sum(dnorm(Data$delta, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$T_diff, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True mu:         ", true.mu, "\n")
cat("Posterior mean:   ", round(fit$Summary2[pos.mu, "Mean"], 3), "\n")
cat("True sigma:      ", true.sigma, "\n")
cat("Posterior mean:   ", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Random Walk with Measurement Error **\[NEW\]**

The local level model treats the observed series \\y_t\\ as a noisy
measurement of a latent random walk \\x_t\\, separating observation
error \\\sigma\_{\text{obs}}\\ from state innovation
\\\sigma\_{\text{state}}\\. The latent states \\x\_{1:T}\\ enter the
model as explicit parameters, making this a joint inference problem over
the hidden trajectory and both variance components. Because each state
\\x_t\\ depends on \\x\_{t-1}\\, the prior on the latent sequence
embodies the random walk structure, while the likelihood connects
observed values to the corresponding states through a simple Gaussian
measurement equation. This formulation follows the dynamic linear model
framework of West and Harrison [\[112\]](#ref112). Signal extraction via
the local level model is applied in hydrology to separate measurement
noise from the underlying water level process.

#### Form

\\y_t \sim \mathcal{N}(x_t, \sigma^2\_{\text{obs}}), \quad t=1,\dots,T\\
\\x_1 \sim \mathcal{N}(0, 100)\\ \\x_t \sim \mathcal{N}(x\_{t-1},
\sigma^2\_{\text{state}}), \quad t=2,\dots,T\\ \\\sigma\_{\text{obs}}
\sim \mathcal{HC}(25)\\ \\\sigma\_{\text{state}} \sim \mathcal{HC}(25)\\

#### model_spec() notation

The latent state-space structure with \\T\\ unknown states as parameters
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42229)
T.obs <- 100

#### True parameter values
true.sigma.obs <- 0.5
true.sigma.state <- 0.2

#### Simulate latent random walk and noisy observations
x.true <- numeric(T.obs)
x.true[1] <- 0
for (t in 2:T.obs) {
    x.true[t] <- x.true[t - 1] + rnorm(1, 0, true.sigma.state)
}
y <- x.true + rnorm(T.obs, 0, true.sigma.obs)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(x=rep(0, T.obs), sigma.obs=0,
    sigma.state=0))
pos.x <- grep("^x\\[", parm.names)
pos.sigma.obs <- grep("sigma.obs", parm.names)
pos.sigma.state <- grep("sigma.state", parm.names)
PGF <- function(Data) {
    sigma.obs <- runif(1)
    sigma.state <- runif(1)
    x <- cumsum(c(0, rnorm(Data$T.obs - 1, 0, sigma.state)))
    return(c(x, sigma.obs, sigma.state))
    }
Data <- list(T.obs=T.obs, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.x=pos.x,
    pos.sigma.obs=pos.sigma.obs, pos.sigma.state=pos.sigma.state,
    y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    x <- parm[Data$pos.x]
    sigma.obs <- interval(parm[Data$pos.sigma.obs], 1e-100, Inf)
    parm[Data$pos.sigma.obs] <- sigma.obs
    sigma.state <- interval(parm[Data$pos.sigma.state], 1e-100, Inf)
    parm[Data$pos.sigma.state] <- sigma.state
    ### Log-Prior (random walk on latent states IS the prior)
    x1.prior <- dnorm(x[1], 0, 10, log=TRUE)
    x.prior <- sum(dnorm(x[2:Data$T.obs], x[1:(Data$T.obs - 1)],
        sigma.state, log=TRUE))
    sigma.obs.prior <- dhalfcauchy(sigma.obs, 25, log=TRUE)
    sigma.state.prior <- dhalfcauchy(sigma.state, 25, log=TRUE)
    ### Log-Likelihood
    LL <- sum(dnorm(Data$y, x, sigma.obs, log=TRUE))
    ### Log-Posterior
    LP <- LL + x1.prior + x.prior + sigma.obs.prior + sigma.state.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$T.obs, x, sigma.obs), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, T.obs), 1, 0.5)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True sigma.obs:    ", true.sigma.obs, "\n")
cat("Posterior mean:     ", round(fit$Summary2[pos.sigma.obs, "Mean"], 3), "\n")
cat("True sigma.state:  ", true.sigma.state, "\n")
cat("Posterior mean:     ", round(fit$Summary2[pos.sigma.state, "Mean"], 3), "\n")

#### Compare estimated latent states to truth
x.post <- fit$Summary2[pos.x, "Mean"]
plot(1:T.obs, x.true, type = "l", col = "black", lwd = 2,
    xlab = "Time", ylab = "Latent state",
    main = "Latent state recovery")
lines(1:T.obs, x.post, col = "blue", lwd = 2)
points(1:T.obs, y, col = "grey70", pch = 16, cex = 0.5)
legend("topleft", legend = c("True x", "Posterior mean x", "Observed y"),
    col = c("black", "blue", "grey70"), lwd = c(2, 2, NA),
    pch = c(NA, NA, 16), bty = "n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### State Space Model (SSM), Linear Regression

The state space linear regression allows regression coefficients to
evolve over time as first-order autoregressive processes centred on
long-run means \\\mu_j\\ with persistence \\\phi_j\\. This captures
structural change in the relationship between predictors and the
response, making it suitable for financial or macroeconomic series where
coefficients drift. State space regression with time-varying
coefficients follows the dynamic linear model framework of West and
Harrison [\[112\]](#ref112) and Durbin and Koopman [\[33\]](#ref33).
Time-varying coefficient models are applied in finance to track the
evolving relationship between market risk factors and asset returns.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2\_{J+1}), \quad
t=1,\dots,T\\ \\\mu = \textbf{X}\beta\\ \\\beta\_{t,j} \sim
\mathcal{N}(\mu_j + \phi_j(\beta\_{t-1,j} - \mu_j), \sigma^2_j), \quad
t=2,\dots,T, \quad j=1,\dots,J\\ \\\beta\_{1,j} \sim \mathcal{N}(\mu_j +
\phi_j(b^0_j - \mu_j), \sigma^2_j), \quad j=1,\dots,J\\ \\b^0_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\mu_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,J\\ \\\phi_j \sim \mathcal{BETA}(20, 1.5) \quad
j=1,\dots,J\\ \\\sigma_j \sim \mathcal{HC}(25), \quad j=1,\dots,(J+1)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the regression coefficients follow latent AR(1) processes
indexed by time, requiring imperative state-propagation logic that falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42201)
T <- 50
J <- 3   # intercept + 2 predictors

#### True parameter values
true.mu    <- c(2.0, 0.5, -0.3)
true.phi   <- c(0.95, 0.90, 0.85)
true.sigma <- c(0.15, 0.10, 0.08, 0.5)  # J state + 1 obs

#### Simulate time-varying coefficients
beta.true <- matrix(0, T, J)
beta.true[1, ] <- true.mu
for (t in 2:T) {
    beta.true[t, ] <- true.mu + true.phi * (beta.true[t - 1, ] - true.mu) +
        rnorm(J, 0, true.sigma[1:J])
}

#### Design matrix and response
X <- cbind(1, matrix(rnorm(T * (J - 1)), T, J - 1))
mu.true <- rowSums(beta.true * X)
y <- rnorm(T, mu.true, true.sigma[J + 1])

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(b0=rep(0,J), beta=matrix(0,T,J),
    mu=rep(0,J), phi=rep(0,J), sigma=rep(0,J+1)))
pos.b0 <- grep("b0", parm.names)
pos.beta <- grep("beta", parm.names)
pos.mu <- grep("mu", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    b0 <- rnorm(Data$J)
    beta <- c(rnorm(Data$T,mean(Data$y),1), rnorm(Data$T*(Data$J-1)))
    mu <- rnorm(Data$J)
    phi <- runif(Data$J, -1, 1)
    sigma <- runif(Data$J+1)
    return(c(beta, mu, phi, sigma))
    }
Data <- list(J=J, PGF=PGF, T=T, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.b0=pos.b0, pos.beta=pos.beta,
    pos.mu=pos.mu, pos.phi=pos.phi, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    b0 <- parm[Data$pos.b0]
    beta <- matrix(parm[Data$pos.beta], Data$T, Data$J)
    mu <- parm[Data$pos.mu]
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], -1, 1)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    b0.prior <- sum(dnormv(b0, 0, 1000, log=TRUE))
    beta.prior <- sum(dnorm(beta, matrix(mu, Data$T, Data$J, byrow=TRUE) +
        matrix(phi, Data$T, Data$J, byrow=TRUE) *
        (rbind(b0, beta[-Data$T,]) -
        matrix(mu, Data$T, Data$J, byrow=TRUE)),
        matrix(sigma[1:Data$J], Data$T, Data$J, byrow=TRUE), log=TRUE))
    mu.prior <- sum(dnormv(mu, 0, 1000, log=TRUE))
    phi.prior <- sum(dbeta((phi+1)/2, 20, 1.5, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- rowSums(beta*Data$X)
    LL <- sum(dnorm(Data$y, mu, sigma[Data$J+1], log=TRUE))
    yhat <- rnorm(length(mu), mu, sigma[Data$J+1]) #Fitted
    #yhat <- rnorm(length(mu), rowSums(matrix(rnorm(Data$T*Data$J,
        # rbind(b0, beta[-Data$T,]), matrix(sigma[-Data$J], Data$T, Data$J,
        # byrow=TRUE)), Data$T, Data$J) * Data$X), sigma[Data$J+1]) #One-step ahead
    ### Log-Posterior
    LP <- LL + b0.prior + beta.prior + mu.prior + phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(mean(y),T), rep(0,T*(J-1)), rep(0,J),
    rep(0,J), rep(1,J+1))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### State Space Model (SSM), Local Level

The local level model is the simplest non-trivial state space model,
consisting of a latent random walk observed with Gaussian noise. The two
variance parameters govern the smoothness of the latent level and the
observation noise respectively. This version treats both variances as
static and includes out-of-sample forecasting for the final
observations. The local level model is the simplest state space model,
introduced by Muth [\[86\]](#ref86) and systematised by Harvey
[\[55\]](#ref55). Local level models are used in signal processing to
extract a slowly varying trend from noisy sensor measurements.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_1), \quad t=1,\dots,T\\
\\\mu_t \sim \mathcal{N}(\mu\_{t-1}, \sigma^2_2), \quad t=2,\dots,T\\
\\\mu_1 \sim \mathcal{N}(0, 1000)\\ \\\sigma_j \sim \mathcal{HC}(25),
\quad j=1,\dots,2\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the latent states \\\mu_t\\ follow a dynamic random walk with
recursive time dependence that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42202)
T <- 20
T.m <- 14  # last observed time; (T.m+2):T are held out

#### True parameter values
true.sigma <- c(0.1, 1.0)  # obs noise, state noise

#### Simulate local level process
mu.orig <- rep(0, T)
for (t in 2:T) mu.orig[t] <- mu.orig[t - 1] + rnorm(1, 0, true.sigma[2])
y <- mu.orig + rnorm(T, 0, true.sigma[1])
y[(T.m + 2):T] <- NA  # hold out for forecasting

#### Assemble Data list
mon.names <- rep(NA, (T - T.m))
for (i in 1:(T - T.m)) mon.names[i] <- paste("yhat[", (T.m + i), "]", sep = "")
parm.names <- as.parm.names(list(mu=rep(0,T), sigma=rep(0,2)))
pos.mu <- grep("mu", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    mu <- rnorm(Data$T)
    sigma <- runif(2)
    return(c(mu, sigma))
    }
Data <- list(PGF=PGF, T=T, T.m=T.m, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.sigma=pos.sigma, y=y)
Dyn <- matrix(paste("mu[",1:T,"]",sep=""), T, 1)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    mu.prior <- sum(dnormv(mu[1], 0, 1000, log=TRUE),
        dnorm(mu[-1], mu[-Data$T], sigma[2], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnorm(Data$y[1:Data$T.m], mu[1:Data$T.m], sigma[1],
        log=TRUE))
    yhat <- rnorm(length(mu), c(mu[1], rnorm(Data$T-1, mu[-Data$T],
        sigma[2])), sigma[1]) #One-step ahead
    ### Log-Posterior
    LP <- LL + mu.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[(Data$T.m+1):Data$T],
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,T), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### State Space Model (SSM), Local Linear Trend

The local linear trend model extends the local level model (section
[link](#ssm.ll)) by adding a dynamic slope parameter \\\delta_t\\ that
itself follows a random walk. This two-component latent structure can
capture both level shifts and changing growth rates, making it
well-suited for trending economic and demographic time series. This
example has static variance parameters. The local linear trend model
extends the local level with a stochastic slope, following Harvey
[\[55\]](#ref55) and Durbin and Koopman [\[33\]](#ref33). Local linear
trend models are applied in demography to forecast population growth
trajectories that exhibit changing growth rates.

#### Form

\\\textbf{y}\_t \sim \mathcal{N}(\mu_t, \sigma^2_1), \quad t=1,\dots,T\\
\\\mu_t \sim \mathcal{N}(\mu\_{t-1} + \delta\_{t-1}, \sigma^2_2), \quad
t=2,\dots,T\\ \\\mu_1 \sim \mathcal{N}(0, 1000)\\ \\\delta_t \sim
\mathcal{N}(\delta\_{t-1}, \sigma^2_3), \quad t=2,\dots,T\\ \\\delta_1
\sim \mathcal{N}(0, 1000)\\ \\\sigma_j \sim \mathcal{HC}(25), \quad
j=1,\dots,3\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because both the level \\\mu_t\\ and slope \\\delta_t\\ are latent
dynamic states with recursive time dependence that falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42203)
T <- 20
T.m <- 14  # last observed time; (T.m+2):T held out

#### True parameter values
true.sigma <- c(0.1, 1.0, 0.1)  # obs, level, slope

#### Simulate local linear trend process
mu.orig <- delta.orig <- rep(0, T)
for (t in 2:T) {
    delta.orig[t] <- delta.orig[t - 1] + rnorm(1, 0, true.sigma[3])
    mu.orig[t] <- mu.orig[t - 1] + delta.orig[t - 1] + rnorm(1, 0, true.sigma[2])
}
y <- mu.orig + rnorm(T, 0, true.sigma[1])
y[(T.m + 2):T] <- NA  # hold out for forecasting

#### Assemble Data list
mon.names <- rep(NA, (T - T.m))
for (i in 1:(T - T.m)) mon.names[i] <- paste("yhat[", (T.m + i), "]", sep = "")
parm.names <- as.parm.names(list(mu=rep(0,T), delta=rep(0,T),
    sigma=rep(0,3)))
pos.mu <- grep("mu", parm.names)
pos.delta <- grep("delta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    mu <- rnorm(Data$T)
    delta <- rnorm(Data$T)
    sigma <- runif(3)
    return(c(mu, delta, sigma))
    }
Data <- list(PGF=PGF, T=T, T.m=T.m, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.delta=pos.delta,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    delta <- parm[Data$pos.delta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    mu.prior <- sum(dnormv(mu[1], 0, 1000, log=TRUE),
        dnorm(mu[-1], mu[-Data$T]+delta[-Data$T], sigma[2],
        log=TRUE))
    delta.prior <- sum(dnormv(delta[1], 0, 1000, log=TRUE),
        dnorm(delta[-1], delta[-Data$T], sigma[3], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnorm(Data$y[1:Data$T.m], mu[1:Data$T.m], sigma[1],
        log=TRUE))
    yhat <- rnorm(length(mu), c(mu[1], rnorm(Data$T-1, mu[-Data$T],
        sigma[2])), sigma[1]) #One-step ahead
    ### Log-Posterior
    LP <- LL + mu.prior + delta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[(Data$T.m+1):Data$T],
        yhat=yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,T), rep(0,T), rep(1,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### State Space Model (SSM), Stochastic Volatility (SV)

The stochastic volatility model treats the log-variance of a financial
return series as a latent AR(1) process, allowing the volatility to
change smoothly over time. Unlike GARCH models that specify variance as
a deterministic function of past returns, the SV model introduces a
genuine state equation for volatility, making it more flexible but
requiring MCMC or particle filtering for estimation. The stochastic
volatility model was introduced by Taylor [\[104\]](#ref104), with
efficient Bayesian estimation developed by Kim et al. [\[63\]](#ref63).
Stochastic volatility models are used in options pricing where the
log-volatility follows a latent autoregressive process.

#### Form

\\\textbf{y} \sim \mathcal{N}(0, \sigma^2)\\ \\\sigma^2 =
\frac{1}{\exp(\theta)}\\ \\\beta = \exp(\mu / 2)\\ \\\theta_1 \sim
\mathcal{N}(\mu + \phi (\alpha - \mu), \tau)\\ \\\theta_t \sim
\mathcal{N}(\mu + \phi (\theta\_{t-1} - \mu), \tau), \quad t=2,\dots,T\\
\\\alpha \sim \mathcal{N}(\mu, \tau)\\ \\\phi \sim \mathcal{U}(-1, 1)\\
\\\mu \sim \mathcal{N}(0, 10)\\ \\\tau \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the latent log-volatility \\\theta_t\\ follows a dynamic AR(1)
process with recursive time dependence that falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42204)
T <- 50

#### True parameter values
true.mu  <- -1.0
true.phi <- 0.95
true.tau <- 0.2

#### Simulate stochastic volatility process
alpha <- rnorm(1, true.mu, true.tau)
theta <- numeric(T)
theta[1] <- rnorm(1, true.mu + true.phi * (alpha - true.mu), true.tau)
for (t in 2:T) {
    theta[t] <- rnorm(1, true.mu + true.phi * (theta[t - 1] - true.mu), true.tau)
}
sigma2 <- 1 / exp(theta)
y <- rnorm(T, 0, sqrt(sigma2))

#### Assemble Data list
mon.names <- c("LP", paste("sigma2[", 1:T, "]", sep = ""))
parm.names <- as.parm.names(list(theta=rep(0,T), alpha=0, phi=0, mu=0,
    tau=0))
pos.theta <- grep("theta", parm.names)
pos.alpha <- grep("alpha", parm.names)
pos.phi <- grep("phi", parm.names)
pos.mu <- grep("mu", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    phi <- runif(1,-1,1)
    mu <- rnorm(1)
    tau <- runif(1)
    alpha <- rnorm(1, mu, tau)
    theta <- rnorm(Data$T, mu + phi*(alpha - mu), tau)
    return(c(theta, alpha, phi, mu, tau))
    }
Data <- list(PGF=PGF, T=T, mon.names=mon.names, parm.names=parm.names,
    pos.theta=pos.theta, pos.alpha=pos.alpha, pos.phi=pos.phi,
    pos.mu=pos.mu, pos.tau=pos.tau, y=y)
Dyn <- matrix(paste("theta[",1:T,"]",sep=""), T, 1)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    theta <- parm[Data$pos.theta]
    alpha <- parm[Data$pos.alpha]
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], -1, 1)
    mu <- parm[Data$pos.mu]
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    ### Log-Prior
    alpha.prior <- dnormv(alpha, mu, tau, log=TRUE)
    theta.prior <- sum(dnormv(theta[1], mu + phi*(alpha-mu), tau,
        log=TRUE), dnormv(theta[-1], mu + phi*(theta[-Data$T]-mu), tau,
        log=TRUE))
    phi.prior <- dunif(phi, -1, 1, log=TRUE)
    mu.prior <- dnormv(mu, 0, 10, log=TRUE)
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    ### Log-Likelihood
    beta <- exp(mu / 2)
    sigma2 <- 1 / exp(theta)
    LL <- sum(dnormv(Data$y, 0, sigma2, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + theta.prior + phi.prior + mu.prior +
        tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, sigma2),
        yhat=rnormv(length(Data$y), 0, sigma2), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,T), 0, 0, 0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]  <- true.mu[seq_along(pos.mu)]
ground_truth[pos.phi] <- true.phi[seq_along(pos.phi)]
ground_truth[pos.tau] <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Nonparametric and semiparametric regression

### Adaptive Logistic Basis (ALB) Regression

Adaptive Logistic Basis (ALB) regression is an essentially automatic
non-parametric approach to estimating the nonlinear relationship between
each of multiple independent variables (IVs) and the dependent variable
(DV). It is automatic because when using the suggested \\K = 2J + 1\\
components (see below) given \\J\\ IVs, the data determines the
nonlinear relationships, whereas in other methods, such as with splines,
the user must specify the number of knots and possibly consider
placement of the knots. Knots do not exist in ALB. Both the DV and IVs
should be centered and scaled. The softmax basis function approach
follows the framework of Denison et al. [\[29\]](#ref29) for flexible
nonparametric regression in a Bayesian setting. In environmental
science, ALB models have been used to capture nonlinear dose-response
relationships between pollutant concentrations and ecological endpoints.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu =
\textbf{S}\delta\\ \\\textbf{S}\_{i,m} =
\frac{\phi\_{i,m}}{\sum^M\_{m=1} \phi\_{i,m}}\\ \\\phi\_{i,m} =
\exp(\alpha_m + \textbf{X}\_{i,1:J}\beta\_{1:J,m}), \quad i=1,\dots,N,
\quad m=1,\dots,M\\ \\\alpha_m \sim \mathcal{N}(0, 10), \quad
m=1,\dots,(M-1)\\ \\\alpha_M = 0\\ \\\beta\_{j,m} \sim \mathcal{N}(0,
100), \quad j=1,\dots,J, \quad m=1,\dots,(M-1)\\ \\\beta\_{j,M} = 0\\
\\\delta_m \sim \mathcal{N}(\zeta, \tau^2), \quad m=1,\dots,M\\ \\\sigma
\sim \mathcal{HC}(25)\\ \\\zeta \sim \mathcal{N}(0, 10)\\ \\\tau \sim
\mathcal{HC}(25)\\

#### Data

``` r

data(demonsnacks)
N <- nrow(demonsnacks)
y <- log(demonsnacks$Calories)
X <- as.matrix(log(demonsnacks[,c(1,4,10)]+1))
J <- ncol(X)
y <- CenterScale(y)
for (j in 1:J) X[,j] <- CenterScale(X[,j])
K <- 2*J+1
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=rep(0,K-1), beta=matrix(0,J,K-1),
    delta=rep(0,K), zeta=0, sigma=0, tau=0))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.delta <- grep("delta", parm.names)
pos.zeta <- grep("zeta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(Data$K-1)
    beta <- rnorm(Data$J*(Data$K-1))
    delta <- rnorm(Data$K)
    zeta <- rnorm(1)
    sigma <- rhalfcauchy(1,5)
    tau <- rhalfcauchy(1,5)
    return(c(alpha, beta, delta, zeta, sigma, tau))
    }
Data <- list(J=J, K=K, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, y=y, pos.alpha=pos.alpha, pos.beta=pos.beta,
    pos.delta=pos.delta, pos.zeta=pos.zeta, pos.sigma=pos.sigma,
    pos.tau=pos.tau)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    zeta <- parm[Data$pos.zeta]
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    beta <- matrix(parm[Data$pos.beta], Data$J, Data$K-1)
    delta <- parm[Data$pos.delta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperprior
    zeta.prior <- dnormv(zeta, 0, 10, log=TRUE)
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    ### Log-Prior
    alpha.prior <- sum(dnormv(alpha, 0, 10, log=TRUE))
    beta.prior <- sum(dnormv(beta, 0, 100, log=TRUE))
    delta.prior <- sum(dnorm(delta, zeta, tau, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    phi <- cbind(exp(matrix(alpha, Data$N, Data$K-1, byrow=TRUE) +
        tcrossprod(Data$X, t(beta))),1)
    mu <- tcrossprod(phi / rowSums(phi), t(delta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + delta.prior + zeta.prior +
        sigma.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K-1), rep(0,J*(K-1)), rep(0,K), 0, 1, 1)
```

#### model_spec() notation

The ALB model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the softmax basis function \\S\_{i,m} = \phi\_{i,m} / \sum_m
\phi\_{i,m}\\ involves a row-wise normalisation of exponentiated linear
predictors that creates a nonlinear, parameter-dependent design matrix.
This compositional transformation with coupled parameters falls outside
the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Inspect softmax basis weights
delta.post <- fit$Summary2[pos.delta, "Mean"]
cat("Posterior mean delta:", round(delta.post, 3), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Fourier Regression **\[NEW\]**

Fourier regression decomposes a periodic signal into a sum of cosine and
sine harmonics at specified frequencies, estimating the amplitude and
phase of each component together with an overall mean and residual
variance, making it a natural model for data exhibiting cyclical
patterns. The design matrix is precomputed from the trigonometric basis
functions, so the model reduces to standard linear regression on the
basis columns, and posterior inference on the harmonic coefficients
\\a_k\\ and \\b_k\\ quantifies the contribution of each frequency to the
observed signal. Fourier analysis of periodic phenomena traces back to
Schuster [\[129\]](#ref129), and the Bayesian spectral framework was
advocated by Jaynes [\[130\]](#ref130). Fourier regression is applied in
climate science to extract annual and semi-annual cycles from
temperature records and quantify their amplitudes with uncertainty.

#### Form

\\y_t \sim \mathcal{N}(\mu_t, \sigma^2), \quad t=1,\dots,T\\ \\\mu_t =
\alpha + \sum\_{k=1}^{K} \left\[ a_k \cos\\\left(\frac{2\pi f_k
t}{T}\right) + b_k \sin\\\left(\frac{2\pi f_k t}{T}\right) \right\]\\
\\\alpha \sim \mathcal{N}(0, 100)\\ \\a_k \sim \mathcal{N}(0, 100),
\quad k=1,\dots,K\\ \\b_k \sim \mathcal{N}(0, 100), \quad k=1,\dots,K\\
\\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

With the trigonometric basis precomputed and stored in the data as a
design matrix `B`, this model reduces to standard normal regression on
the basis columns. The
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL handles it directly.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma)
  mu = B %*% coef
  coef[j] ~ NormalV(0, 100), j = 1,...,P
  sigma ~ HalfCauchy(25)
")
code(spec)
```

#### Ground truth and data

``` r

set.seed(42213)
T <- 300
K <- 2       # Number of harmonics
f <- c(1, 3) # Frequencies

#### True parameter values
true.alpha <- 2.0
true.a <- c(1.5, 0.8)
true.b <- c(-0.5, 0.3)
true.sigma <- 0.5

#### Precompute trigonometric basis
t.seq <- 1:T
B <- cbind(1,
    cos(2*pi*f[1]*t.seq/T), sin(2*pi*f[1]*t.seq/T),
    cos(2*pi*f[2]*t.seq/T), sin(2*pi*f[2]*t.seq/T))
P <- ncol(B) # 2K + 1 = 5
true.coef <- c(true.alpha, true.a, true.b)

#### Generate response
mu.true <- B %*% true.coef
y <- rnorm(T, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(coef=rep(0,P), sigma=0))
pos.coef <- grep("coef", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    coef <- rnorm(Data$P, 0, 0.5)
    sigma <- runif(1)
    return(c(coef, sigma))
    }
Data <- list(T=T, P=P, B=B, y=y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.coef=pos.coef, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    coef <- parm[Data$pos.coef]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    coef.prior <- sum(dnormv(coef, 0, 100, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- Data$B %*% coef
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + coef.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$T, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, P), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True coef: ", true.coef, "\n")
cat("Post. mean:", round(fit$Summary2[pos.coef, "Mean"], 3), "\n")
cat("True sigma:", true.sigma, "\n")
cat("Post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.coef]  <- true.coef[seq_along(pos.coef)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[grep("a", parm.names)]     <- true.a[seq_along(grep("a", parm.names))]
ground_truth[grep("alpha", parm.names)] <- true.alpha[seq_along(grep("alpha", parm.names))]
ground_truth[grep("b", parm.names)]     <- true.b[seq_along(grep("b", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Gompertz Growth Model **\[NEW\]**

The Gompertz growth curve, introduced by Gompertz [\[131\]](#ref131) as
a mortality law and later adapted to biological growth by Laird
[\[132\]](#ref132), is an asymmetric sigmoid defined by \\\mu_i = K
\exp(-b \exp(-c\\ t_i))\\. Unlike the logistic curve, the Gompertz
inflection point lies at roughly one-third of the asymptote rather than
one-half, making it more appropriate for growth processes that
decelerate early. The three positive parameters govern the upper bound
\\K\\, the initial displacement \\b\\, and the intrinsic growth rate
\\c\\, while observation noise enters through a normal likelihood with
unknown standard deviation \\\sigma\\. Gompertz models are widely
applied in oncology to describe tumour volume dynamics and in population
ecology for density-dependent growth.

#### Form

\\y_i \sim \mathcal{N}(\mu_i,\\ \sigma^2), \quad i=1,\dots,N\\ \\\mu_i =
K \exp\\\bigl(-b \exp(-c\\ t_i)\bigr)\\ \\K \sim \mathcal{HC}(25)\\ \\b
\sim \mathcal{HC}(25)\\ \\c \sim \mathcal{HC}(25)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The Gompertz growth curve involves a double-exponential nonlinear mean
function \\K \exp(-b \exp(-c\\ t))\\ which cannot be expressed as a
linear predictor in
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function handles this directly.

#### Ground truth and data

``` r

set.seed(42215)
N <- 100

#### True parameter values
true.K <- 10
true.b <- 5
true.c <- 0.3
true.sigma <- 0.3

#### Generate time grid and response
t.seq <- seq(0, 20, length.out = N)
mu.true <- true.K * exp(-true.b * exp(-true.c * t.seq))
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(K=0, b=0, c=0, sigma=0))
pos.K <- grep("K", parm.names)
pos.b <- grep("b", parm.names)
pos.c <- grep("c", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    K <- runif(1, 5, 15)
    b <- runif(1, 1, 10)
    c <- runif(1, 0.1, 1)
    sigma <- runif(1, 0.1, 1)
    return(c(K, b, c, sigma))
    }
Data <- list(N=N, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.K=pos.K, pos.b=pos.b, pos.c=pos.c, pos.sigma=pos.sigma,
    t.seq=t.seq, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    K <- interval(parm[Data$pos.K], 1e-100, Inf)
    parm[Data$pos.K] <- K
    b <- interval(parm[Data$pos.b], 1e-100, Inf)
    parm[Data$pos.b] <- b
    cc <- interval(parm[Data$pos.c], 1e-100, Inf)
    parm[Data$pos.c] <- cc
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    K.prior <- dhalfcauchy(K, 25, log=TRUE)
    b.prior <- dhalfcauchy(b, 25, log=TRUE)
    cc.prior <- dhalfcauchy(cc, 25, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- K * exp(-b * exp(-cc * Data$t.seq))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + K.prior + b.prior + cc.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(8, 3, 0.5, 0.5)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("K     -- true:", true.K,
    " post. mean:", round(fit$Summary2[pos.K, "Mean"], 3), "\n")
cat("b     -- true:", true.b,
    " post. mean:", round(fit$Summary2[pos.b, "Mean"], 3), "\n")
cat("c     -- true:", true.c,
    " post. mean:", round(fit$Summary2[pos.c, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.K]     <- true.K[seq_along(pos.K)]
ground_truth[pos.b]     <- true.b[seq_along(pos.b)]
ground_truth[pos.c]     <- true.c[seq_along(pos.c)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Penalized Spline Regression

Penalized spline regression places a set of radial basis functions (or
truncated polynomial bases) at fixed knot locations and controls
overfitting through a hierarchical prior on the basis coefficients. The
smoothing penalty is automatically determined by the prior variance
\\\sigma^2_2\\ of the spline coefficients, making the effective degrees
of freedom a data-driven quantity. Penalized spline regression was
introduced by Eilers and Marx [\[34\]](#ref34), with the Bayesian mixed
model interpretation following Ruppert et al. [\[97\]](#ref97).
Penalized splines are used in growth curve analysis to model the
nonlinear relationship between age and height in paediatric cohorts.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2_1)\\ \\\mu = \textbf{X}
\beta + \textbf{S}\\ \\\textbf{S} = \textbf{Z} \gamma\\
\\\textbf{Z}\_{i,k} = \left\\ \begin{array}{l l} (\textbf{x}\_i - k)^D &
\quad \mbox{if \$\textbf{Z}\_{i,k} \> 0\$}\\ 0 \\ \end{array} \right. \\
\\\beta_d \sim \mathcal{N}(0, 1000), \quad d=1,\dots,(D+1)\\ \\\gamma_k
\sim \mathcal{N}(0, \sigma^2_2), \quad k=1,\dots,K\\ \\\sigma_j \sim
\mathcal{HC}(25), \quad j=1,\dots,2\\

#### model_spec

The penalized spline model requires a custom truncated-polynomial basis
matrix \\\textbf{Z}\\ constructed from user-chosen knots and polynomial
degree, which falls outside the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The manual Model function is therefore the recommended approach.

#### Ground truth and data

``` r

set.seed(42160)
N <- 200
K <- 10 #Number of knots
D <- 2  #Degree of polynomial

#### True function: smooth sinusoidal with known sigma
true.sigma <- 0.3
x.raw <- seq(0, 1, length.out = N)
f.true <- sin(2 * pi * x.raw) + 0.5 * cos(4 * pi * x.raw)
y <- f.true + rnorm(N, 0, true.sigma)

#### Build design and basis matrices
x <- CenterScale(x.raw)
k <- as.vector(quantile(x, probs = (1:K / (K + 1))))
X <- cbind(1, matrix(x, N, D))
for (d in 1:D) {X[, d + 1] <- X[, d + 1]^d}
Z <- matrix(x, N, K) - matrix(k, N, K, byrow = TRUE)
Z <- ifelse(Z > 0, Z, 0); Z <- Z^D
mon.names <- c("LP", paste("S[", 1:nrow(X), "]", sep = ""))
parm.names <- as.parm.names(list(beta = rep(0, 1 + D), gamma = rep(0, K),
    log.sigma = rep(0, 2)))
pos.beta  <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta  <- rnorm(1 + Data$D)
    gamma <- rnorm(Data$K)
    sigma <- runif(2)
    return(c(beta, gamma, sigma))
    }
Data <- list(D = D, K = K, N = N, PGF = PGF, Z = Z, X = X,
    mon.names = mon.names, parm.names = parm.names,
    pos.beta = pos.beta, pos.gamma = pos.gamma,
    pos.sigma = pos.sigma, y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma[2], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    S <- as.vector(tcrossprod(Data$Z, t(gamma)))
    mu <- as.vector(tcrossprod(Data$X, t(beta))) + S
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,S),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,1+D), rep(0,K), c(1,1))
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "sigma"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Penalized Spline, 2D Tensor Product **\[NEW\]**

Two-dimensional penalized splines extend univariate P-splines to
surfaces by forming a tensor product of two marginal basis matrices, one
for each covariate, so that the fitted surface \\\mu_i = \alpha +
\textbf{B}^{\otimes}\_i \gamma\\ flexibly captures interaction effects
across the entire domain. Each marginal basis uses truncated quadratic
functions placed at quantile knots, and the row-wise Kronecker product
of the two evaluations yields an \\N \times K_1 K_2\\ design matrix
whose coefficients \\\gamma\\ are penalised toward zero through a normal
prior with shared standard deviation \\\sigma\_\gamma\\. Tensor product
splines were developed within the P-spline framework of Eilers and Marx
[\[34\]](#ref34) and extended to multiple dimensions by Wood
[\[158\]](#ref158). Two-dimensional P-splines are used in environmental
science to model spatial fields such as temperature as a smooth function
of latitude and longitude.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_i = \alpha +
\textbf{B}^{\otimes}\_i \gamma, \quad i=1,\dots,N\\
\\\textbf{B}^{\otimes}\_{i,\cdot} = \textbf{B}^1\_{i,\cdot} \otimes
\textbf{B}^2\_{i,\cdot}\\ \\\textbf{B}^1\_{i,j} = \max(0, (x\_{1,i} -
k^1_j))^2, \quad j=1,\dots,K_1\\ \\\textbf{B}^2\_{i,k} = \max(0,
(x\_{2,i} - k^2_k))^2, \quad k=1,\dots,K_2\\ \\\alpha \sim
\mathcal{N}(0, 1000)\\ \\\gamma_m \sim \mathcal{N}(0, \sigma^2\_\gamma),
\quad m=1,\dots,K_1 K_2\\ \\\sigma \sim \mathcal{HC}(25)\\
\\\sigma\_\gamma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The tensor product P-spline requires constructing a custom basis matrix
from two marginal spline bases via a row-wise Kronecker product, which
falls outside the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL. The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42227)
N <- 200
K1 <- 5
K2 <- 5
Ktp <- K1 * K2

#### True surface: f(x1,x2) = sin(2*pi*x1) * cos(2*pi*x2)
true.sigma <- 0.3
x1 <- runif(N)
x2 <- runif(N)
f.true <- sin(2 * pi * x1) * cos(2 * pi * x2)
y <- f.true + rnorm(N, 0, true.sigma)

#### Build marginal basis matrices
k1 <- as.vector(quantile(x1, probs = 1:K1 / (K1 + 1)))
k2 <- as.vector(quantile(x2, probs = 1:K2 / (K2 + 1)))
B1 <- outer(x1, k1, function(x, k) pmax(0, (x - k)^2))
B2 <- outer(x2, k2, function(x, k) pmax(0, (x - k)^2))

#### Tensor product: row-wise Kronecker product
B_tp <- matrix(0, N, Ktp)
for (i in seq_len(N)) {
    B_tp[i, ] <- kronecker(B1[i, ], B2[i, ])
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, gamma=rep(0,Ktp),
    sigma=0, sigma.gamma=0))
pos.alpha <- grep("alpha", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.sigma <- grep("^sigma$", parm.names)
pos.sigma.gamma <- grep("sigma.gamma", parm.names)
PGF <- function(Data) {
    alpha <- rnorm(1)
    gamma <- rnorm(Data$Ktp, 0, 0.1)
    sigma <- runif(1)
    sigma.gamma <- runif(1)
    return(c(alpha, gamma, sigma, sigma.gamma))
    }
Data <- list(N=N, Ktp=Ktp, PGF=PGF, B_tp=B_tp, mon.names=mon.names,
    parm.names=parm.names, pos.alpha=pos.alpha, pos.gamma=pos.gamma,
    pos.sigma=pos.sigma, pos.sigma.gamma=pos.sigma.gamma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- parm[Data$pos.alpha]
    gamma <- parm[Data$pos.gamma]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    sigma.gamma <- interval(parm[Data$pos.sigma.gamma], 1e-100, Inf)
    parm[Data$pos.sigma.gamma] <- sigma.gamma
    ### Log-Prior
    alpha.prior <- dnormv(alpha, 0, 1000, log=TRUE)
    gamma.prior <- sum(dnorm(gamma, 0, sigma.gamma, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    sigma.gamma.prior <- dhalfcauchy(sigma.gamma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- alpha + as.vector(Data$B_tp %*% gamma)
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + alpha.prior + gamma.prior + sigma.prior + sigma.gamma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, rep(0, Ktp), 1, 1)
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Penalized Spline, Heteroskedastic **\[NEW\]**

The heteroskedastic penalized spline fits dual smoothing models for both
the mean and the log-variance of a response, using truncated polynomial
bases with Bayesian penalisation through random-effect priors on the
spline coefficients. The mean function \\\mu_i = \textbf{X}\_\mu\beta +
\textbf{Z}\_\mu\gamma\\ captures the smooth trend, while
\\\log(\sigma_i) = \textbf{X}\_v\delta + \textbf{Z}\_v\eta\\ allows the
residual variance to vary smoothly across the covariate domain. The
penalty standard deviations \\\sigma\_\gamma\\ and \\\sigma\_\eta\\
control the roughness of each spline, with larger values permitting more
flexibility. The P-spline framework was introduced by Eilers and Marx
[\[34\]](#ref34), and the Bayesian extension to heteroskedastic
smoothing was developed by Lang and Brezger [\[157\]](#ref157).
Heteroskedastic spline models are used in income distribution modelling
where both the conditional mean and variance of log-earnings depend on
age.

#### Form

\\y_i \sim \mathcal{N}(\mu_i, \sigma_i^2), \quad i=1,\dots,N\\ \\\mu_i =
\textbf{X}\_\mu \beta + \textbf{Z}\_\mu \gamma\\ \\\log(\sigma_i) =
\textbf{X}\_v \delta + \textbf{Z}\_v \eta\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,D+1\\ \\\gamma_k \sim
\mathcal{N}(0, \sigma\_\gamma^2), \quad k=1,\dots,K\_\mu\\ \\\delta_j
\sim \mathcal{N}(0, 1000), \quad j=1,\dots,D+1\\ \\\eta_k \sim
\mathcal{N}(0, \sigma\_\eta^2), \quad k=1,\dots,K_v\\ \\\sigma\_\gamma
\sim \mathcal{HC}(25)\\ \\\sigma\_\eta \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the dual spline structure with separate penalised bases for the
mean and the log-variance introduces a heteroskedastic likelihood with
design matrices constructed from truncated polynomial knots, which falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42226)
N <- 200
D <- 2
K.mu <- 8
K.var <- 5

#### Covariate and true functions
x <- seq(0, 1, length.out=N)
f.mu <- sin(2 * pi * x)
f.log.sigma <- -0.5 + x
y <- rnorm(N, f.mu, exp(f.log.sigma))

#### Build truncated polynomial bases
knots.mu <- seq(0, 1, length.out=K.mu + 2)[-c(1, K.mu + 2)]
knots.var <- seq(0, 1, length.out=K.var + 2)[-c(1, K.var + 2)]
X.mu <- cbind(1, x, x^2)
Z.mu <- sapply(knots.mu, function(k) pmax(x - k, 0)^D)
X.var <- cbind(1, x, x^2)
Z.var <- sapply(knots.var, function(k) pmax(x - k, 0)^D)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0, D+1), gamma=rep(0, K.mu),
    delta=rep(0, D+1), eta=rep(0, K.var), sigma.gamma=0, sigma.eta=0))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.eta <- grep("eta", parm.names)
pos.sigma.gamma <- grep("sigma.gamma", parm.names)
pos.sigma.eta <- grep("sigma.eta", parm.names)
PGF <- function(Data) {
    return(c(rnorm(Data$D + 1, 0, 0.1), rnorm(Data$K.mu, 0, 0.1),
        rnorm(Data$D + 1, 0, 0.1), rnorm(Data$K.var, 0, 0.1),
        runif(1), runif(1)))
    }
Data <- list(D=D, K.mu=K.mu, K.var=K.var, N=N, PGF=PGF, X.mu=X.mu,
    X.var=X.var, Z.mu=Z.mu, Z.var=Z.var, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.delta=pos.delta,
    pos.eta=pos.eta, pos.gamma=pos.gamma, pos.sigma.eta=pos.sigma.eta,
    pos.sigma.gamma=pos.sigma.gamma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    delta <- parm[Data$pos.delta]
    eta <- parm[Data$pos.eta]
    sigma.gamma <- interval(parm[Data$pos.sigma.gamma], 1e-100, Inf)
    parm[Data$pos.sigma.gamma] <- sigma.gamma
    sigma.eta <- interval(parm[Data$pos.sigma.eta], 1e-100, Inf)
    parm[Data$pos.sigma.eta] <- sigma.eta
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnorm(gamma, 0, sigma.gamma, log=TRUE))
    delta.prior <- sum(dnormv(delta, 0, 1000, log=TRUE))
    eta.prior <- sum(dnorm(eta, 0, sigma.eta, log=TRUE))
    sigma.gamma.prior <- dhalfcauchy(sigma.gamma, 25, log=TRUE)
    sigma.eta.prior <- dhalfcauchy(sigma.eta, 25, log=TRUE)
    ### Log-Likelihood
    mu <- as.vector(Data$X.mu %*% beta + Data$Z.mu %*% gamma)
    log.sigma <- as.vector(Data$X.var %*% delta + Data$Z.var %*% eta)
    sigma.i <- exp(log.sigma)
    LL <- sum(dnorm(Data$y, mu, sigma.i, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + delta.prior + eta.prior +
        sigma.gamma.prior + sigma.eta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu, sigma.i), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, D+1), rep(0, K.mu), rep(0, D+1), rep(0, K.var),
    1, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Spatial and spatio-temporal models

### Conditional Autoregression (CAR), Poisson

The conditional autoregressive (CAR) Poisson model extends Poisson
regression with spatially structured random effects, where each area’s
error term is conditionally normal given its neighbors’ errors weighted
by an adjacency matrix. The spatial autocorrelation parameter \\\rho\\
controls the strength of borrowing across contiguous areas, making this
formulation appropriate for disease mapping and other areal count data
where neighboring regions are expected to share similar risk profiles.
The conditional autoregressive (CAR) model was introduced by Besag
[\[10\]](#ref10) for modelling spatial dependence on irregular lattices.
CAR models are central to disease mapping in public health, where
county-level disease counts are spatially smoothed to identify
geographic clusters.

#### Form

\\\textbf{y} \sim \mathcal{P}(\lambda)\\ \\\lambda =
\exp(\log(\textbf{E}) + \beta_1 + \beta_2 \textbf{x} + \epsilon)\\
\\\epsilon \sim \mathcal{N}(\epsilon\_\mu, \sigma^2)\\
\\\epsilon\_{\mu\[i\]} = \rho \sum^J\_{j=1} \textbf{A}\_{i,j}
\epsilon_i, \quad i=1,\dots,N\\ \\\beta_j \sim \mathcal{N}(0, 1000),
\quad j=1,\dots,J\\ \\\rho \sim \mathcal{U}(-1,1)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The CAR Poisson model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatially structured error term \\\epsilon\_\mu\\
depends on the adjacency matrix \\\textbf{A}\\ and a row-wise weighted
sum of neighboring errors, which requires imperative constructs that
fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42125)
N <- 56

#### Build adjacency matrix (irregular lattice)
A <- matrix(0, N, N)
A[1,c(5,9,11,19)] <- 1; A[2,c(7,10)] <- 1
A[3,c(6,12)] <- 1; A[4,c(18,20,28)] <- 1; A[5,c(1,11,12,13,19)] <- 1
A[6,c(3,8)] <- 1; A[7,c(2,10,13,16,17)] <- 1; A[8,6] <- 1
A[9,c(1,11,17,19,23,29)] <- 1; A[10,c(2,7,16,22)] <- 1
A[11,c(1,5,9,12)] <- 1; A[12,c(3,5,11)] <- 1; A[13,c(5,7,17,19)] <- 1
A[14,c(31,32,35)] <- 1; A[15,c(25,29,50)] <- 1
A[16,c(7,10,17,21,22,29)] <- 1; A[17,c(7,9,13,16,19,29)] <- 1
A[18,c(4,20,28,33,55,56)] <- 1; A[19,c(1,5,9,13,17)] <- 1
A[20,c(4,18,55)] <- 1; A[21,c(16,29,50)] <- 1; A[22,c(10,16)] <- 1
A[23,c(9,29,34,36,37,39)] <- 1; A[24,c(27,30,31,44,47,48,55,56)] <- 1
A[25,c(15,26,29)] <- 1; A[26,c(25,29,42,43)] <- 1
A[27,c(24,31,32,55)] <- 1; A[28,c(4,18,33,45)] <- 1
A[29,c(9,15,16,17,21,23,25,26,34,43,50)] <- 1
A[30,c(24,38,42,44,45,56)] <- 1; A[31,c(14,24,27,32,35,46,47)] <- 1
A[32,c(14,27,31,35)] <- 1; A[33,c(18,28,45,56)] <- 1
A[34,c(23,29,39,40,42,43,51,52,54)] <- 1; A[35,c(14,31,32,37,46)] <- 1
A[36,c(23,37,39,41)] <- 1; A[37,c(23,35,36,41,46)] <- 1
A[38,c(30,42,44,49,51,54)] <- 1; A[39,c(23,34,36,40,41)] <- 1
A[40,c(34,39,41,49,52)] <- 1; A[41,c(36,37,39,40,46,49,53)] <- 1
A[42,c(26,30,34,38,43,51)] <- 1; A[43,c(26,29,34,42)] <- 1
A[44,c(24,30,38,48,49)] <- 1; A[45,c(28,30,33,56)] <- 1
A[46,c(31,35,37,41,47,53)] <- 1; A[47,c(24,31,46,48,49,53)] <- 1
A[48,c(24,44,47,49)] <- 1; A[49,c(38,40,41,44,47,48,52,53,54)] <- 1
A[50,c(15,21,29)] <- 1; A[51,c(34,38,42,54)] <- 1
A[52,c(34,40,49,54)] <- 1; A[53,c(41,46,47,49)] <- 1
A[54,c(34,38,49,51,52)] <- 1; A[55,c(18,20,24,27,56)] <- 1
A[56,c(18,24,30,33,45,55)] <- 1

#### True parameter values
true.beta <- c(0.3, 0.15)
true.rho <- 0.6
true.sigma <- 0.4

#### Simulate expected counts and predictor
E <- rgamma(N, shape = 5, rate = 0.5)
x <- sample(c(0, 1, 7, 10, 16, 24), N, replace = TRUE)

#### Generate spatially correlated errors via CAR structure
epsilon <- rnorm(N, 0, true.sigma)
for (iter in 1:50) {
    epsilon.mu <- true.rho * as.vector(A %*% epsilon) /
        pmax(rowSums(A), 1)
    epsilon <- rnorm(N, epsilon.mu, true.sigma)
}

#### Generate Poisson counts
lambda.true <- exp(log(E) + true.beta[1] + true.beta[2] * x / 10 + epsilon)
y <- rpois(N, lambda.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,2), epsilon=rep(0,N), rho=0,
    sigma=0))
pos.beta <- grep("beta", parm.names)
pos.epsilon <- grep("epsilon", parm.names)
pos.rho <- grep("rho", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(2)
    epsilon <- rnorm(Data$N)
    rho <- runif(1,-1,1)
    sigma <- runif(1)
    return(c(beta, epsilon, rho, sigma))
    }
Data <- list(A=A, E=E, N=N, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.epsilon=pos.epsilon,
    pos.rho=pos.rho, pos.sigma=pos.sigma, x=x, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    epsilon <- parm[Data$pos.epsilon]
    parm[Data$pos.rho] <- rho <- interval(parm[Data$pos.rho], -1, 1)
    epsilon.mu <- rho * rowSums(epsilon * Data$A)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    epsilon.prior <- sum(dnorm(epsilon, epsilon.mu, sigma, log=TRUE))
    rho.prior <- dunif(rho, -1, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    lambda <- exp(log(Data$E) + beta[1] + beta[2]*Data$x/10 + epsilon)
    LL <- sum(dpois(Data$y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + epsilon.prior + rho.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(length(lambda), lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,2), rep(0,N), 0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("beta[1] -- true:", true.beta[1],
    " post. mean:", round(fit$Summary2[1, "Mean"], 3), "\n")
cat("beta[2] -- true:", true.beta[2],
    " post. mean:", round(fit$Summary2[2, "Mean"], 3), "\n")
cat("rho    -- true:", true.rho,
    " post. mean:", round(fit$Summary2[pos.rho, "Mean"], 3), "\n")
cat("sigma  -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.rho]   <- true.rho[seq_along(pos.rho)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Geographically Weighted Regression

Geographically weighted regression (GWR) allows regression coefficients
to vary continuously over space by weighting observations according to
their distance from each location, with bandwidth \\h\\ controlling the
rate of spatial decay. The model incorporates location-specific
heteroscedasticity through gamma-distributed weights \\\nu\\ and a
Student-t-like error structure controlled by shape parameter \\\alpha\\,
making it robust to spatial outliers while capturing non-stationary
relationships. Geographically weighted regression was introduced by
Brunsdon et al. [\[17\]](#ref17) to allow regression coefficients to
vary spatially. GWR is applied in urban studies to reveal how the
relationship between house prices and property features varies across a
metropolitan area.

#### Form

\\\textbf{y}\_{i,k} \sim \mathcal{N}(\mu\_{i,k}, \tau^{-1}\_{i,k}),
\quad i=1,\dots,N, \quad k=1,\dots,N\\ \\\mu\_{i,1:N} = \textbf{X}
\beta\_{i,1:J}\\ \\\tau = \frac{1}{\sigma^2} \textbf{w} \nu\\
\\\textbf{w} = \frac{\exp(-0.5 \textbf{Z}^2)}{\textbf{h}}\\ \\\alpha
\sim \mathcal{U}(1.5, 100)\\ \\\beta\_{i,j} \sim \mathcal{N}(0, 1000),
\quad i=1,\dots,N, \quad j=1,\dots,J\\ \\\textbf{h} \sim
\mathcal{N}(0.1, 1000) \in \[0.1, \infty\]\\ \\\nu\_{i,k} \sim
\mathcal{G}(\alpha, 2), \quad i=1,\dots,N, \quad k=1,\dots,N\\
\\\sigma_i \sim \mathcal{HC}(25), \quad i=1,\dots,N\\

#### model_spec() notation

The GWR model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatially weighted precision structure requires
computing pairwise distance-based weights, location-specific regression
coefficients, and gamma-distributed heteroscedasticity terms, all of
which involve imperative array operations that fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42131)
N <- 49
J <- 3

#### Simulate spatial coordinates
easting <- runif(N, 24, 52)
northing <- runif(N, 25, 45)

#### True spatially varying coefficients (smooth surface)
true.beta0 <- 20 + 0.5 * (easting - 38)
true.beta1 <- -0.8 + 0.02 * (northing - 35)
true.beta2 <- 0.3 - 0.01 * (easting - 38)
true.alpha <- 10
true.sigma <- rep(3.0, N)

#### Generate predictors and response
income <- rnorm(N, 15, 5)
housing <- rnorm(N, 35, 15)
X <- cbind(1, income, housing)

#### Location-specific mean
mu.true <- numeric(N)
for (i in 1:N) {
    mu.true[i] <- true.beta0[i] + true.beta1[i] * income[i] +
        true.beta2[i] * housing[i]
}
y.vec <- rnorm(N, mu.true, true.sigma)

#### Build distance and weight matrices
D <- as.matrix(dist(cbind(northing, easting), diag=TRUE, upper=TRUE))
Z <- D / sd(as.vector(D))
y <- matrix(0, N, N)
for (i in 1:N) { for (k in 1:N) { y[i,k] <- y.vec[k] } }
mon.names <- "LP"
parm.names <- as.parm.names(list(alpha=0, beta=matrix(0,N,J), H=0,
    nu=matrix(0,N,N), sigma=rep(0,N)))
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.H <- grep("H", parm.names)
pos.nu <- grep("nu", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    alpha <- runif(1,1.5,100)
    beta <- rnorm(Data$N*Data$J)
    H <- runif(1,0.1,1000)
    nu <- rgamma(Data$N*Data$N,alpha,2)
    sigma <- runif(Data$N)
    return(c(alpha, beta, H, nu, sigma))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, Z=Z, latitude=northing,
    longitude=easting, mon.names=mon.names, parm.names=parm.names,
    pos.alpha=pos.alpha, pos.beta=pos.beta, pos.H=pos.H, pos.nu=pos.nu,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    alpha <- interval(parm[Data$pos.alpha], 1.5, 100)
    parm[Data$pos.alpha] <- alpha
    beta <- matrix(parm[Data$pos.beta], Data$N, Data$J)
    parm[Data$pos.H] <- H <- interval(parm[Data$pos.H], 0.1, Inf)
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    nu <- matrix(nu, Data$N, Data$N)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    alpha.prior <- dunif(alpha, 1.5, 100, log=TRUE)
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    h.prior <- dhalfnorm(H-0.1, 1000, log=TRUE)
    nu.prior <- sum(dgamma(nu, alpha, 2, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    w <- exp(-0.5 * Data$Z^ 2) / H
    tau <- (1/sigma^ 2) * w * nu
    mu <- tcrossprod(Data$X, beta)
    LL <- sum(dnormp(Data$y, mu, tau, log=TRUE))
    #WSE <- w * nu * (Data$y - mu)^ 2; w.y <- w * nu * Data$y
    #WMSE <- rowMeans(WSE); y.w <- rowSums(w.y) / rowSums(w)
    #LAR2 <- 1 - WMSE / sd(y.w)^ 2
    ### Log-Posterior
    LP <- LL + alpha.prior + beta.prior + h.prior + nu.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnormp(prod(dim(mu)), mu, tau), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(runif(1,1.5,100), rep(0,N*J), 1, rep(1,N*N), rep(1,N))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.alpha] <- true.alpha[seq_along(pos.alpha)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[grep("beta0", parm.names)] <- true.beta0[seq_along(grep("beta0", parm.names))]
ground_truth[grep("beta1", parm.names)] <- true.beta1[seq_along(grep("beta1", parm.names))]
ground_truth[grep("beta2", parm.names)] <- true.beta2[seq_along(grep("beta2", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Kriging

This is an example of universal kriging of \\\textbf{y}\\ given
\\\textbf{X}\\, regression effects \\\beta\\, and spatial effects
\\\zeta\\. Euclidean distance between spatial coordinates (longitude and
latitude) is used for each of \\i=1,\dots,N\\ records of \\\textbf{y}\\.
An additional record is created from the same data-generating process to
compare the accuracy of interpolation. For the spatial component,
\\\phi\\ is the rate of spatial decay and \\\kappa\\ is the scale.
\\\kappa\\ is often difficult to identify, so it is set to 1 (Gaussian),
but may be allowed to vary up to 2 (Exponential). In practice, \\\phi\\
is also often difficult to identify. While \\\Sigma\\ is spatial
covariance, spatial correlation is \\\rho = \exp(-\phi \textbf{D})\\. To
extend this to a large data set, consider the predictive process kriging
example in section [link](#kriging.pp). Kriging was developed by
Matheron [\[78\]](#ref78), building on the pioneering work of Krige
[\[68\]](#ref68) for spatial interpolation in geostatistics. Kriging is
the standard tool in mining for estimating ore grades at unsampled
locations from borehole assay data.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2_1)\\ \\ \mu = \textbf{X}
\beta + \zeta\\ \\ \textbf{y}^{new} = \textbf{X} \beta + \sum^N\_{i=1}
\left ( \frac{\rho_i}{\sum \rho} \zeta_i \right )\\ \\ \rho = \exp(-\phi
\textbf{D}^{new})^\kappa\\ \\ \zeta \sim \mathcal{N}\_N(\zeta\_\mu,
\Sigma)\\ \\ \Sigma = \sigma^2_2 \exp(-\phi \textbf{D})^\kappa\\ \\
\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,2\\ \\ \sigma_j \sim
\mathcal{HC}(25) \in \[0.1,10\], \quad j=1,\dots,2\\ \\ \phi \sim
\mathcal{U}(1, 5)\\ \\ \zeta\_\mu = 0\\ \\ \kappa = 1\\

#### model_spec() notation

The universal kriging model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatial effects \\\zeta\\ follow a multivariate normal
with a covariance matrix that depends on the distance matrix and spatial
decay parameter \\\phi\\, requiring imperative matrix construction and
interpolation logic that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42133)
N <- 20

#### True parameter values
true.beta <- c(50, 2)
true.sigma <- c(2.0, 5.0)  # observation noise, spatial scale
true.phi <- 2.0
true.kappa <- 1

#### Simulate spatial coordinates (N training + 1 test)
longitude <- runif(N + 1, 0, 100)
latitude <- runif(N + 1, 0, 100)
#### Generate spatial effects from the true covariance
D.all <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))
Sigma.true <- true.sigma[2]^2 * exp(-true.phi * D.all)^true.kappa
zeta.all <- as.vector(rmvn(1, rep(0, N + 1), Sigma.true))

#### Generate response
X <- cbind(1, runif(N + 1, -2, 2))
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- rnorm(N + 1, mu.true + zeta.all, true.sigma[1])

#### Split into training and test
longitude.new <- longitude[N + 1]; latitude.new <- latitude[N + 1]
Xnew <- X[N + 1, ]; ynew <- y[N + 1]
longitude <- longitude[1:N]; latitude <- latitude[1:N]
X <- X[1:N, ]; y <- y[1:N]
D <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))
D.new <- sqrt((longitude - longitude.new)^2 + (latitude - latitude.new)^2)

#### Assemble Data list
mon.names <- c("LP","ynew")
parm.names <- as.parm.names(list(zeta=rep(0,N), beta=rep(0,2),
    sigma=rep(0,2), phi=0))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.phi <- grep("phi", parm.names)
PGF <- function(Data) {
    beta <- rnorm(2)
    sigma <- runif(2,0.1,10)
    phi <- runif(1,1,5)
    kappa <- 1
    zeta <- rmvn(1, rep(0,Data$N),
        sigma[2]*sigma[2]*exp(-phi*Data$D)^ kappa)
    return(c(zeta, beta, sigma, phi))
    }
Data <- list(D=D, D.new=D.new, latitude=latitude, longitude=longitude,
    N=N, PGF=PGF, X=X, Xnew=Xnew, mon.names=mon.names,
    parm.names=parm.names, pos.zeta=pos.zeta, pos.beta=pos.beta,
    pos.sigma=pos.sigma, pos.phi=pos.phi, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    zeta <- parm[Data$pos.zeta]
    kappa <- 1
    sigma <- interval(parm[Data$pos.sigma], 0.1, 10)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
    Sigma <- sigma[2]*sigma[2] * exp(-phi * Data$D)^ kappa
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    zeta.prior <- dmvn(zeta, rep(0, Data$N), Sigma, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
    phi.prior <- dunif(phi, 1, 5, log=TRUE)
    ### Interpolation
    rho <- exp(-phi * Data$D.new)^ kappa
    ynew <- rnorm(1, sum(beta * Data$Xnew) + sum(rho / sum(rho) * zeta),
        sigma[1])
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta)) + zeta
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + zeta.prior + sigma.prior + phi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,N), rep(0,2), rep(1,2), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("beta    -- true:", true.beta, "\n")
cat("             est:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("sigma   -- true:", true.sigma, "\n")
cat("             est:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
cat("phi     -- true:", true.phi, "\n")
cat("             est:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")
cat("y.new   -- true:", round(ynew, 3),
    " post. mean:", round(mean(fit$Monitor[, 2]), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[grep("kappa", parm.names)] <- true.kappa[seq_along(grep("kappa", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Kriging, Predictive Process

The predictive process approximation to kriging uses the first \\K\\ of
\\N\\ records as knots for the parent spatial process and projects the
remaining locations through the knot-based covariance structure,
substantially reducing the computational cost of the full \\N \times N\\
kriging covariance from \\O(N^3)\\ to \\O(K^3)\\. This makes it
practical for moderately large spatial datasets while retaining the
geostatistical interpretation of the full kriging model. The predictive
process approximation was proposed by Banerjee et al. [\[3\]](#ref3) to
make kriging computationally feasible for large spatial datasets.
Predictive process kriging is used in air quality monitoring to
interpolate pollutant concentrations across a country from a sparse
network of monitoring stations.

#### Form

\\ \textbf{y} \sim \mathcal{N}(\mu, \sigma^2_1)\\ \\ \mu\_{1:K} =
\textbf{X}\_{1:K,1:J} \beta + \zeta\\ \\ \mu\_{(K+1):N} =
\textbf{X}\_{(K+1):N,1:J} \beta + \sum^{N-K}\_{p=1}
\frac{\lambda\_{p,1:K}}{\sum^{N-K}\_{q=1} \lambda\_{q,1:K}} \zeta^T\\ \\
\lambda = \exp(-\phi \textbf{D}\_P)^\kappa\\ \\ \textbf{y}^{new} =
\textbf{X} \beta + \sum^K\_{k=1} (\frac{\rho_k}{\sum \rho} \zeta_k)\\ \\
\rho = \exp(-\phi \textbf{D}^{new})^\kappa\\ \\ \zeta \sim
\mathcal{N}\_K(0, \Sigma)\\ \\ \Sigma = \sigma^2_2 \exp(-\phi
\textbf{D})^\kappa\\ \\ \beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,2\\ \\ \sigma_j \sim \mathcal{HC}(25), \quad j=1,\dots,2\\ \\
\phi \sim \mathrm{N}(0, 1000) \in \[1, 5\]\\ \\ \kappa = 1\\

#### model_spec() notation

The predictive process kriging model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the knot-based spatial projection requires computing
distance matrices between knots and prediction locations, constructing
the predictive process covariance, and performing matrix operations that
fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42134)
N <- 100
K <- 30
longitude <- runif(N+1,0,100)
latitude <- runif(N+1,0,100)
D <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
Sigma <- 10000 * exp(-1.5 * D)
zeta <- colMeans(rmvn(1000, rep(0,N+1), Sigma))
beta <- c(50,2)
X <- matrix(runif((N+1)*2,-2,2),(N+1),2); X[,1] <- 1
mu <- as.vector(tcrossprod(X, t(beta)))
y <- mu + zeta
longitude.new <- longitude[N+1]; latitude.new <- latitude[N+1]
Xnew <- X[N+1,]; ynew <- y[N+1]
longitude <- longitude[1:N]; latitude <- latitude[1:N]
X <- X[1:N,]; y <- y[1:N]
D <- as.matrix(dist(cbind(longitude[1:K],latitude[1:K]), diag=TRUE,
    upper=TRUE))
D.P <- matrix(0, N-K, K)
for (i in (K+1):N) {
    D.P[K+1-i,] <- sqrt((longitude[1:K] - longitude[i])^ 2 +
        (latitude[1:K] - latitude[i])^ 2)}
D.new <- sqrt((longitude[1:K] - longitude.new)^ 2 +
    (latitude[1:K] - latitude.new)^ 2)
mon.names <- c("LP","ynew")
parm.names <- as.parm.names(list(zeta=rep(0,K), beta=rep(0,2),
    sigma=rep(0,2), phi=0))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.phi <- grep("phi", parm.names)
PGF <- function(Data) {
    beta <- rnorm(2)
    sigma <- runif(2,0.1,10)
    phi <- runif(1,1,5)
    kappa <- 1
    zeta <- rmvn(1, rep(0,Data$K),
        sigma[2]*sigma[2]*exp(-phi*Data$D)^ kappa)
    return(c(zeta, beta, sigma, phi))
    }
Data <- list(D=D, D.new=D.new, D.P=D.P, K=K, N=N, PGF=PGF, X=X,
    Xnew=Xnew, latitude=latitude, longitude=longitude,
    mon.names=mon.names, parm.names=parm.names, pos.zeta=pos.zeta,
    pos.beta=pos.beta, pos.sigma=pos.sigma, pos.phi=pos.phi, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    zeta <- parm[Data$pos.zeta]
    kappa <- 1
    sigma <- interval(parm[Data$pos.sigma], 1, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
    Sigma <- sigma[2]*sigma[2] * exp(-phi * Data$D)^ kappa
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    zeta.prior <- dmvn(zeta, rep(0, Data$K), Sigma, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
    phi.prior <- dunif(phi, 1, 5, log=TRUE)
    ### Interpolation
    rho <- exp(-phi * Data$D.new)^ kappa
    ynew <- rnorm(1, sum(beta * Data$Xnew) + sum(rho / sum(rho) * zeta),
        sigma)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    mu[1:Data$K] <- mu[1:Data$K] + zeta
    lambda <- exp(-phi * Data$D.P)^ kappa
    mu[(Data$K+1):Data$N] <- mu[(Data$K+1):Data$N] +
        rowSums(lambda / rowSums(lambda) *
        matrix(zeta, Data$N - Data$K, Data$K, byrow=TRUE))
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + zeta.prior + sigma.prior + phi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew),
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,K), c(mean(y), 0), rep(1,2), 3)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("beta    -- true:", true.beta, "\n")
cat("             est:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("sigma   -- true:", true.sigma, "\n")
cat("             est:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
cat("y.new   -- true:", round(ynew, 3),
    " post. mean:", round(mean(fit$Monitor[, 2]), 3), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Log-Gaussian Cox Process **\[NEW\]**

The log-Gaussian Cox process (LGCP), formalised by Moller et
al. [\[163\]](#ref163), is a doubly stochastic point process in which
the log-intensity surface is a Gaussian random field. On a discretised
spatial domain of \\G\\ grid cells, counts \\y_g\\ in each cell follow a
Poisson distribution whose rate is the product of cell area and a
spatially correlated intensity \\\lambda_g = \exp(\mu + w_g)\\, where
\\\mu\\ is a global intercept and \\w_g\\ are latent Gaussian process
values with an exponential covariance kernel parameterised by marginal
standard deviation \\\sigma_w\\ and length-scale \\\ell\\. The model
jointly estimates the intensity surface and the covariance
hyperparameters from the observed counts. LGCPs are applied in spatial
ecology for modelling species occurrence intensity and in epidemiology
for disease mapping.

#### Form

\\y_g \sim \text{Poisson}(\lambda_g \cdot A_g), \quad g=1,\dots,G\\
\\\log(\lambda_g) = \mu + w_g\\ \\\textbf{w} \sim
\mathcal{N}\_G(\textbf{0},\\ \textbf{K})\\ \\K\_{ij} = \sigma_w^2
\exp\\\left(-\frac{d\_{ij}}{\ell}\right)\\ \\\mu \sim \mathcal{N}(0,
100)\\ \\\sigma_w \sim \mathcal{HC}(5)\\ \\\ell \sim \mathcal{HC}(5)\\

#### model_spec() notation

The LGCP cannot be expressed in
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the latent Gaussian process with a spatial covariance matrix
that depends on estimated hyperparameters \\\sigma_w\\ and \\\ell\\
requires imperative matrix construction inside the model function, which
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42219)
G <- 25

#### True parameter values
true.mu <- 1
true.sigma.w <- 1
true.ell <- 2

#### Generate 5x5 grid coordinates
grid.x <- rep(1:5, each = 5)
grid.y <- rep(1:5, times = 5)
coords <- cbind(grid.x, grid.y)
d.mat <- as.matrix(dist(coords))
A <- 1  # unit cell area

#### Generate latent GP and Poisson counts
K.true <- true.sigma.w^2 * exp(-d.mat / true.ell)
w.true <- as.vector(rmvn(1, rep(0, G), K.true))
lambda.true <- exp(true.mu + w.true) * A
y <- rpois(G, lambda.true)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(mu=0, sigma.w=0, ell=0, w=rep(0, G)))
pos.mu <- grep("^mu$", parm.names)
pos.sigma.w <- grep("sigma.w", parm.names)
pos.ell <- grep("ell", parm.names)
pos.w <- grep("^w", parm.names)
PGF <- function(Data) {
    mu <- rnorm(1, 0, 2)
    sigma.w <- runif(1, 0.5, 2)
    ell <- runif(1, 0.5, 4)
    w <- rnorm(Data$G)
    return(c(mu, sigma.w, ell, w))
    }
Data <- list(A=A, G=G, PGF=PGF, d.mat=d.mat, mon.names=mon.names,
    parm.names=parm.names, pos.mu=pos.mu, pos.sigma.w=pos.sigma.w,
    pos.ell=pos.ell, pos.w=pos.w, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    mu <- parm[Data$pos.mu]
    sigma.w <- interval(parm[Data$pos.sigma.w], 1e-100, Inf)
    parm[Data$pos.sigma.w] <- sigma.w
    ell <- interval(parm[Data$pos.ell], 1e-100, Inf)
    parm[Data$pos.ell] <- ell
    w <- parm[Data$pos.w]
    ### Log-Prior
    mu.prior <- dnormv(mu, 0, 100, log=TRUE)
    sigma.w.prior <- dhalfcauchy(sigma.w, 5, log=TRUE)
    ell.prior <- dhalfcauchy(ell, 5, log=TRUE)
    ## GP prior: rebuild covariance from current hyperparameters
    K <- sigma.w^2 * exp(-Data$d.mat / ell)
    w.prior <- dmvn(w, rep(0, Data$G), K, log=TRUE)
    ### Log-Likelihood
    lambda <- exp(mu + w) * Data$A
    LL <- sum(dpois(Data$y, lambda, log=TRUE))
    ### Log-Posterior
    LP <- LL + mu.prior + sigma.w.prior + ell.prior + w.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rpois(Data$G, lambda), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 1, 1, rep(0, G))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("mu      -- true:", true.mu,
    " post. mean:", round(fit$Summary2[pos.mu, "Mean"], 3), "\n")
cat("sigma.w -- true:", true.sigma.w,
    " post. mean:", round(fit$Summary2[pos.sigma.w, "Mean"], 3), "\n")
cat("ell     -- true:", true.ell,
    " post. mean:", round(fit$Summary2[pos.ell, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.ell] <- true.ell[seq_along(pos.ell)]
ground_truth[pos.mu]  <- true.mu[seq_along(pos.mu)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Network Autocorrelation **\[NEW\]**

The network autocorrelation model extends linear regression by
incorporating spatial dependence through a row-normalised adjacency
matrix \\\textbf{W}\\. The autoregressive parameter \\\rho\\ captures
the degree to which an observation is influenced by the weighted average
of its network neighbours, so that the response vector is determined
jointly by the covariates and the network-filtered feedback. The
reduced-form likelihood involves the determinant of \\(\textbf{I} -
\rho\textbf{W})\\, which encodes the global connectivity structure. The
spatial autoregressive model was introduced by Ord [\[91\]](#ref91)
within a maximum likelihood framework and later extended to Bayesian
estimation. Network autocorrelation models are used in social science to
estimate peer effects in friendship networks and contagion dynamics.

#### Form

\\\textbf{y} = (\textbf{I} - \rho\textbf{W})^{-1}(\textbf{X}\beta +
\epsilon), \quad \epsilon \sim \mathcal{N}(\textbf{0},
\sigma^2\textbf{I})\\ \\\log L = -\frac{N}{2}\log(\sigma^2) +
\log\|\det(\textbf{I} - \rho\textbf{W})\| -
\frac{1}{2\sigma^2}\textbf{e}'\textbf{e}\\ \\\textbf{e} = (\textbf{I} -
\rho\textbf{W})\textbf{y} - \textbf{X}\beta\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\rho \sim \mathcal{U}(-1,
1)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the likelihood depends on the log-determinant of \\(\textbf{I} -
\rho\textbf{W})\\ and the spatial lag transformation of the response
vector, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42223)
N <- 50
J <- 2
true.beta <- c(3, 1.5)
true.rho <- 0.4
true.sigma <- 1

#### Generate random network adjacency (Erdos-Renyi, row-normalised)
W <- matrix(rbinom(N * N, 1, 0.15), N, N)
diag(W) <- 0
rs <- rowSums(W)
rs[rs == 0] <- 1
W <- W / rs

#### Simulate from reduced form
X <- cbind(1, rnorm(N))
A <- diag(N) - true.rho * W
y <- as.vector(solve(A, X %*% true.beta + rnorm(N, 0, true.sigma)))

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0, J), rho=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.rho <- grep("rho", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    return(c(rnorm(Data$J), runif(1, -1, 1), runif(1)))
    }
Data <- list(J=J, N=N, PGF=PGF, W=W, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.rho=pos.rho,
    pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    rho <- interval(parm[Data$pos.rho], -1, 1)
    parm[Data$pos.rho] <- rho
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    rho.prior <- dunif(rho, -1, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    A <- diag(Data$N) - rho * Data$W
    e <- as.vector(A %*% Data$y - Data$X %*% beta)
    log.det.A <- as.numeric(determinant(A, logarithm=TRUE)$modulus)
    LL <- -Data$N / 2 * log(sigma^2) + log.det.A -
        sum(e^2) / (2 * sigma^2)
    ### Log-Posterior
    LP <- LL + beta.prior + rho.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=as.vector(solve(A, Data$X %*% beta +
            rnorm(Data$N, 0, sigma))),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0, 0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.rho]   <- true.rho[seq_along(pos.rho)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Space-Time, Dynamic

The dynamic space-time model applies kriging to a stationary spatial
component and then allows the spatial parameters to evolve through
state-space dynamics. This approach to space-time or spatiotemporal
modeling applies kriging to a stationary spatial component for points in
space \\s=1,\dots,S\\ first at time \\t=1\\, where space is continuous
and time is discrete. Vector \\\zeta\\ contains these spatial effects.
Next, SSM (State Space Model) or DLM (Dynamic Linear Model) components
are applied to the spatial parameters (\\\phi\\, \\\kappa\\, and
\\\lambda\\) and regression effects (\\\beta\\). These parameters are
allowed to vary dynamically with time \\t=2,\dots,T\\, and the resulting
spatial process is estimated for each of these time-periods. When time
is discrete, a dynamic space-time process can be applied. The matrix
\\\Theta\\ contains the dynamically varying stationary spatial effects,
or space-time effects. Spatial coordinates are given in longitude and
latitude for \\s=1,\dots,S\\ points in space and measurements are taken
across discrete time-periods \\t=1,\dots,T\\ for \\\textbf{Y}\_{s,t}\\.
The dependent variable is also a function of design matrix
\\\textbf{X}\\ (which may also be dynamic, but is static in this
example) and dynamic regression effects matrix \\\beta\_{1:K,1:T}\\. For
more information on kriging, see section [link](#kriging). For more
information on SSMs or DLMs, see section [link](#ssm.lin.reg). To extend
this to a large spatial data set, consider incorporating the predictive
process kriging example in section [link](#kriging.pp). Dynamic
space-time models extend the state space framework to spatiotemporal
processes, following Cressie and Wikle [\[28\]](#ref28). Dynamic
space-time models are used in oceanography to forecast sea surface
temperature fields that evolve both spatially and temporally.

#### Form

\\\textbf{Y}\_{s,t} \sim \mathcal{N}(\mu\_{s,t}, \sigma^2_1), \quad
s=1,\dots,S, \quad t=1,\dots,T\\ \\\mu\_{s,t} = \textbf{X}\_{s,1:K}
\beta\_{1:K,t} + \Theta\_{s,t}\\ \\\Theta\_{s,t} =
\frac{\Sigma\_{s,s,t}}{\sum^S\_{r=1} \Sigma\_{r,s,t}} \Theta\_{s,t-1},
\quad s=1,\dots,S, \quad t=2,\dots,T\\ \\\Theta\_{s,1} = \zeta_s\\
\\\zeta \sim \mathcal{N}\_S(0, \Sigma\_{1:S,1:S,1})\\
\\\Sigma\_{1:S,1:S,t} = \lambda^2_t \exp(-\phi_t
\textbf{D})^{\kappa\[t\]}\\ \\\sigma_j \sim \mathcal{HC}(25), \quad
j=1,\dots,4\\ \\\beta\_{k,1} \sim \mathcal{N}(0, 1000), \quad
k=1,\dots,K\\ \\\beta\_{k,t} \sim \mathcal{N}(\beta\_{k,t-1}, \tau^2_k),
\quad k=1,\dots,K, \quad t=2,\dots,T\\ \\\phi_1 \sim
\mathcal{HN}(1000)\\ \\\phi_t \sim \mathcal{N}(\phi\_{t-1}, \sigma^2_2)
\in \[0,\infty\], \quad t=2,\dots,T\\ \\\kappa_1 \sim
\mathcal{HN}(1000)\\ \\\kappa_t \sim \mathcal{N}(\kappa\_{t-1},
\sigma^2_3) \in \[0,\infty\], \quad t=2,\dots,T\\ \\\lambda_1 \sim
\mathcal{HN}(1000)\\ \\\lambda_t \sim \mathcal{N}(\lambda\_{t-1},
\sigma^2_4) \in \[0,\infty\], \quad t=2,\dots,T\\

#### model_spec() notation

The dynamic space-time model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the time-varying spatial covariance parameters (\\\phi_t\\,
\\\kappa_t\\, \\\lambda_t\\) and the recursive space-time effects
\\\Theta\\ require imperative loops and dynamic matrix construction that
fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42189)

#### Simulation dimensions
S <- 20  # Number of sites
T <- 13  # Number of time-periods
K <- 3   # Number of predictors including intercept

#### True parameter values
true.sigma <- c(2.0, 0.5, 0.5, 0.5)

#### Simulate spatial coordinates and distance matrix
longitude <- runif(S, -100, -95)
latitude <- runif(S, 30, 35)
D <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))

#### Static predictors and response
X <- cbind(1, rnorm(S), rnorm(S))
phi.true <- rep(1.5, T); kappa.true <- rep(1.0, T); lambda.true <- rep(3.0, T)
Sigma.1 <- lambda.true[1]^2 * exp(-phi.true[1] * D)^kappa.true[1]
zeta <- as.vector(rmvn(1, rep(0, S), Sigma.1))
Theta <- matrix(zeta, S, T)
for (t in 2:T) {
    Sigma.t <- lambda.true[t]^2 * exp(-phi.true[t] * D)^kappa.true[t]
    for (s in 1:S) {
        Theta[s, t] <- sum(Sigma.t[, s] / sum(Sigma.t[, s]) * Theta[, t-1])
    }
}
beta.true <- matrix(c(rep(50, T), rep(0.5, T), rep(-0.3, T)), K, T, byrow=TRUE)
mu.true <- tcrossprod(X, t(beta.true)) + Theta
Y <- matrix(rnorm(S * T, mu.true, true.sigma[1]), S, T)
mon.names <- "LP"
parm.names <- as.parm.names(list(zeta=rep(0,S), beta=matrix(0,K,T),
    phi=rep(0,T), kappa=rep(0,T), lambda=rep(0,T), sigma=rep(0,4),
    tau=rep(0,K)))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.kappa <- grep("kappa", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K*Data$T, rbind(mean(Data$Y),
        matrix(0, Data$K-1, Data$T)), 1)
    phi <- rhalfnorm(Data$T, 1)
    kappa <- rhalfnorm(Data$T, 1)
    lambda <- rhalfnorm(Data$T, 1)
    Sigma <- lambda[1]*lambda[1]*exp(-phi[1]*Data$D)^ kappa[1]
    zeta <- as.vector(rmvn(1, rep(0,Data$S), Sigma))
    sigma <- runif(4)
    tau <- runif(Data$K)
    return(c(zeta, beta, phi, kappa, lambda, sigma, tau))
    }
Data <- list(D=D, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y, latitude=latitude,
    longitude=longitude, mon.names=mon.names, parm.names=parm.names,
    pos.zeta=pos.zeta, pos.beta=pos.beta, pos.phi=pos.phi,
    pos.kappa=pos.kappa, pos.lambda=pos.lambda, pos.sigma=pos.sigma,
    pos.tau=pos.tau)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$K, Data$T)
    zeta <- parm[Data$pos.zeta]
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
    kappa <- interval(parm[Data$pos.kappa], 1e-100, Inf)
    parm[Data$pos.kappa] <- kappa
    lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    Sigma <- array(0, dim=c(Data$S, Data$S, Data$T))
    for (t in 1:Data$T) {
        Sigma[ , ,t] <- lambda[t]^ 2 * exp(-phi[t] * Data$D)^ kappa[t]}
    ### Log-Prior
    beta.prior <- sum(dnormv(beta[,1], 0, 1000, log=TRUE),
        dnorm(beta[,-1], beta[,-Data$T], matrix(tau, Data$K,
        Data$T-1), log=TRUE))
    zeta.prior <- dmvn(zeta, rep(0,Data$S), Sigma[ , , 1], log=TRUE)
    phi.prior <- sum(dhalfnorm(phi[1], sqrt(1000), log=TRUE),
        dtrunc(phi[-1], "norm", a=0, b=Inf, mean=phi[-Data$T],
        sd=sigma[2], log=TRUE))
    kappa.prior <- sum(dhalfnorm(kappa[1], sqrt(1000), log=TRUE),
        dtrunc(kappa[-1], "norm", a=0, b=Inf, mean=kappa[-Data$T],
        sd=sigma[3], log=TRUE))
    lambda.prior <- sum(dhalfnorm(lambda[1], sqrt(1000), log=TRUE),
        dtrunc(lambda[-1], "norm", a=0, b=Inf, mean=lambda[-Data$T],
        sd=sigma[4], log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    tau.prior <- sum(dhalfcauchy(tau, 25, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    Theta <- matrix(zeta, Data$S, Data$T)
    for (t in 2:Data$T) {
        for (s in 1:Data$S) {
    Theta[s,t] <- sum(Sigma[,s,t] / sum(Sigma[,s,t]) * Theta[,t-1])}}
    mu <- mu + Theta
    LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + zeta.prior + sum(phi.prior) +
        sum(kappa.prior) + sum(lambda.prior) + sigma.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,S), rep(c(mean(Y),rep(0,K-1)),T), rep(1,T),
    rep(1,T), rep(1,T), rep(1,4), rep(1,K))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Space-Time, Nonseparable

This approach to space-time or spatiotemporal modeling applies kriging
both to the stationary spatial and temporal components, where space is
continuous and time is discrete. Matrix \\\Xi\\ contains the space-time
effects. Spatial coordinates are given in longitude and latitude for
\\s=1,\dots,S\\ points in space and measurements are taken across
time-periods \\t=1,\dots,T\\ for \\\textbf{Y}\_{s,t}\\. The dependent
variable is also a function of design matrix \\\textbf{X}\\ and
regression effects vector \\\beta\\. For more information on kriging,
see section [link](#kriging). This example uses a nonseparable,
stationary covariance function in which space and time are separable
only when \\\psi=0\\. To extend this to a large space-time data set,
consider incorporating the predictive process kriging example in section
[link](#kriging.pp). Nonseparable space-time covariance functions were
developed by Gneiting [\[52\]](#ref52) to capture interactions between
spatial and temporal dependence. Nonseparable models are applied in air
pollution monitoring where the spatial correlation of pollutant
concentrations depends on the time of day.

#### Form

\\\textbf{Y}\_{s,t} \sim \mathcal{N}(\mu\_{s,t}, \sigma^2_1), \quad
s=1,\dots,S, \quad t=1,\dots,T\\ \\\mu = \textbf{X} \beta + \Xi\\ \\\Xi
\sim \mathcal{N}\_{ST}(\Xi\_\mu, \Sigma)\\ \\ \Sigma = \sigma^2_2 \exp
\left (-\frac{\textbf{D}\_S}{\phi_1}^\kappa -
\frac{\textbf{D}\_T}{\phi_2}^\lambda - \psi
\frac{\textbf{D}\_S}{\phi_1}^\kappa \frac{\textbf{D}\_T}{\phi_2}^\lambda
\right )\\ \\\beta_k \sim \mathcal{N}(0, 1000), \quad k=1,\dots,K\\
\\\phi_j \sim \mathcal{U}(1, 5), \quad j=1,\dots,2\\ \\\sigma_j \sim
\mathcal{HC}(25), \quad j=1,\dots,2\\ \\\psi \sim \mathcal{HC}(25)\\
\\\Xi\_\mu = 0\\ \\\kappa = 1, \quad \lambda = 1\\

#### model_spec() notation

The nonseparable space-time model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the joint space-time covariance function involves a
nonseparable interaction term \\\psi\\ coupling spatial and temporal
distance matrices, requiring imperative matrix operations that fall
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42190)

#### Simulation dimensions
S <- 10  # Number of sites
T <- 13  # Number of time-periods
K <- 3   # Number of predictors including intercept

#### True parameter values
true.beta <- c(50, 0.5, -0.3)
true.sigma <- c(2.0, 5.0)
true.phi <- c(2.0, 2.0)
true.psi <- 0.5

#### Simulate spatial coordinates
longitude <- runif(S, -100, -95)
latitude <- runif(S, 30, 35)

#### Build expanded distance matrices
D.S <- as.matrix(dist(cbind(rep(longitude, T), rep(latitude, T)),
    diag=TRUE, upper=TRUE))
D.T <- as.matrix(dist(cbind(rep(1:T, each=S), rep(1:T, each=S)),
    diag=TRUE, upper=TRUE))

#### Generate space-time effects from the true nonseparable covariance
kappa <- 1; lambda <- 1
Sigma.true <- true.sigma[2]^2 * exp(-(D.S / true.phi[1])^kappa -
    (D.T / true.phi[2])^lambda -
    true.psi * (D.S / true.phi[1])^kappa * (D.T / true.phi[2])^lambda)
Xi.true <- as.vector(rmvn(1, rep(0, S * T), Sigma.true))

#### Static predictors and response
X <- cbind(1, rnorm(S), rnorm(S))
mu.true <- as.vector(tcrossprod(X, t(true.beta))) + matrix(Xi.true, S, T)
Y <- matrix(rnorm(S * T, mu.true, true.sigma[1]), S, T)
mon.names <- "LP"
parm.names <- as.parm.names(list(Xi=matrix(0,S,T), beta=rep(0,K),
    phi=rep(0,2), sigma=rep(0,2), psi=0))
pos.Xi <- grep("Xi", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.psi <- grep("psi", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
    phi <- runif(2,1,5)
    sigma <- runif(2)
    psi <- runif(1)
    kappa <- 1
    lambda <- 1
    Sigma <- sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^ kappa -
        (Data$D.T / phi[2])^ lambda -
        psi*(Data$D.S / phi[1])^ kappa * (Data$D.T / phi[2])^ lambda)
    Xi <- as.vector(rmvn(1, rep(0,Data$S*Data$T), Sigma))
    return(c(Xi, beta, phi, sigma, psi))
    }
Data <- list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
    latitude=latitude, longitude=longitude, mon.names=mon.names,
    parm.names=parm.names, pos.Xi=pos.Xi, pos.beta=pos.beta,
    pos.phi=pos.phi, pos.sigma=pos.sigma, pos.psi=pos.psi)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    Xi.mu <- rep(0,Data$S*Data$T)
    ### Parameters
    beta <- parm[Data$pos.beta]
    Xi <- parm[Data$pos.Xi]
    kappa <- 1; lambda <- 1
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
    parm[Data$pos.psi] <- psi <- interval(parm[Data$pos.psi], 1e-100, Inf)
    Sigma <- sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^ kappa -
        (Data$D.T / phi[2])^ lambda -
        psi*(Data$D.S / phi[1])^ kappa * (Data$D.T / phi[2])^ lambda)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    Xi.prior <- dmvn(Xi, Xi.mu, Sigma, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    phi.prior <- sum(dunif(phi, 1, 5, log=TRUE))
    psi.prior <- dhalfcauchy(psi, 25, log=TRUE)
    ### Log-Likelihood
    Xi <- matrix(Xi, Data$S, Data$T)
    mu <- as.vector(tcrossprod(Data$X, t(beta))) + Xi
    LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + Xi.prior + sigma.prior + phi.prior + psi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,S*T), c(mean(Y),rep(0,K-1)), rep(1,2), rep(1,2),
    1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.psi]   <- true.psi[seq_along(pos.psi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Space-Time, Separable

This introductory approach to space-time or spatiotemporal modeling
applies kriging both to the stationary spatial and temporal components,
where space is continuous and time is discrete. Vector \\\zeta\\
contains the spatial effects and vector \\\theta\\ contains the temporal
effects. Spatial coordinates are given in longitude and latitude for
\\s=1,\dots,S\\ points in space and measurements are taken across
time-periods \\t=1,\dots,T\\ for \\\textbf{Y}\_{s,t}\\. The dependent
variable is also a function of design matrix \\\textbf{X}\\ and
regression effects vector \\\beta\\. For more information on kriging,
see section [link](#kriging). This example uses separable space-time
covariances, which is more convenient but usually less appropriate than
a nonseparable covariance function. To extend this to a large space-time
data set, consider incorporating the predictive process kriging example
in section [link](#kriging.pp). Separable space-time models assume the
covariance factorises into spatial and temporal components, following
the framework of Cressie and Wikle [\[28\]](#ref28). Separable models
are used in precision agriculture to map crop yield variation across
fields and growing seasons.

#### Form

\\\textbf{Y}\_{s,t} \sim \mathcal{N}(\mu\_{s,t}, \sigma^2_1), \quad
s=1,\dots,S, \quad t=1,\dots,T\\ \\\mu\_{s,t} = \textbf{X}\_{s,1:J}
\beta + \zeta_s + \Theta\_{s,t}\\ \\\Theta\_{s,1:T} = \theta\\ \\\theta
\sim \mathcal{N}\_N(\theta\_\mu, \Sigma_T)\\ \\\Sigma_T = \sigma^2_3
\exp(-\phi_2 \textbf{D}\_T)^\lambda\\ \\\zeta \sim
\mathcal{N}\_N(\zeta\_\mu, \Sigma_S)\\ \\\Sigma_S = \sigma^2_2
\exp(-\phi_1 \textbf{D}\_S)^\kappa\\ \\\beta_j \sim \mathcal{N}(0,
1000), \quad j=1,\dots,2\\ \\\sigma_k \sim \mathcal{HC}(25), \quad
k=1,\dots,3\\ \\\phi_k \sim \mathcal{U}(1, 5), \quad k=1,\dots,2\\
\\\zeta\_\mu = 0\\ \\\theta\_\mu = 0\\ \\\kappa = 1, \quad \lambda = 1\\

#### model_spec() notation

The separable space-time model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatial and temporal random effects require multivariate
normal priors with covariance matrices constructed from distance
matrices and decay parameters, involving imperative matrix operations
that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42191)

#### Simulation dimensions
S <- 20  # Number of sites
T <- 13  # Number of time-periods
K <- 3   # Number of predictors including intercept

#### True parameter values
true.beta <- c(50, 0.5, -0.3)
true.sigma <- c(2.0, 4.0, 3.0)
true.phi <- c(2.0, 2.0)

#### Simulate spatial coordinates
longitude <- runif(S, -100, -95)
latitude <- runif(S, 30, 35)

#### Build distance matrices
D.S <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))
D.T <- as.matrix(dist(cbind(c(1:T), c(1:T)), diag=TRUE, upper=TRUE))

#### Generate spatial and temporal effects
kappa.true <- 1; lambda.true <- 1
Sigma.S.true <- true.sigma[2]^2 * exp(-true.phi[1] * D.S)^kappa.true
Sigma.T.true <- true.sigma[3]^2 * exp(-true.phi[2] * D.T)^lambda.true
zeta.true <- as.vector(rmvn(1, rep(0, S), Sigma.S.true))
theta.true <- as.vector(rmvn(1, rep(0, T), Sigma.T.true))

#### Static predictors and response
X <- cbind(1, rnorm(S), rnorm(S))
Theta.true <- matrix(theta.true, S, T, byrow=TRUE)
mu.true <- as.vector(tcrossprod(X, t(true.beta))) + zeta.true + Theta.true
Y <- matrix(rnorm(S * T, mu.true, true.sigma[1]), S, T)
mon.names <- "LP"
parm.names <- as.parm.names(list(zeta=rep(0,S), theta=rep(0,T),
    beta=rep(0,K), phi=rep(0,2), sigma=rep(0,3)))
pos.zeta <- grep("zeta", parm.names)
pos.theta <- grep("theta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
    phi <- runif(2,1,5)
    sigma <- runif(3)
    kappa <- 1
    lambda <- 1
    Sigma.S <- sigma[2]^ 2 * exp(-phi[1] * Data$D.S)^ kappa
    Sigma.T <- sigma[3]^ 2 * exp(-phi[2] * Data$D.T)^ lambda
    zeta <- as.vector(rmvn(1, rep(0,Data$S), Sigma.S))
    theta <- as.vector(rmvn(1, rep(0,Data$T), Sigma.T))
    return(c(zeta, theta, beta, phi, sigma))
    }
Data <- list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
    latitude=latitude, longitude=longitude, mon.names=mon.names,
    parm.names=parm.names, pos.zeta=pos.zeta, pos.theta=pos.theta,
    pos.beta=pos.beta, pos.phi=pos.phi, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    zeta.mu <- rep(0,Data$S)
    theta.mu <- rep(0,Data$T)
    ### Parameters
    beta <- parm[Data$pos.beta]
    zeta <- parm[Data$pos.zeta]
    theta <- parm[Data$pos.theta]
    kappa <- 1; lambda <- 1
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
    Sigma.S <- sigma[2]^ 2 * exp(-phi[1] * Data$D.S)^ kappa
    Sigma.T <- sigma[3]^ 2 * exp(-phi[2] * Data$D.T)^ lambda
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    zeta.prior <- dmvn(zeta, zeta.mu, Sigma.S, log=TRUE)
    theta.prior <- dmvn(theta, theta.mu, Sigma.T, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    phi.prior <- sum(dunif(phi, 1, 5, log=TRUE))
    ### Log-Likelihood
    Theta <- matrix(theta, Data$S, Data$T, byrow=TRUE)
    mu <- as.vector(tcrossprod(Data$X, t(beta))) + zeta + Theta
    LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + zeta.prior + theta.prior + sigma.prior +
        phi.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,S), rep(0,T), rep(0,K), rep(1,2), rep(1,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Spatial Autoregression (SAR)

The spatial autoregressive (SAR) model captures spatial dependence by
including a spatially lagged dependent variable \\\textbf{z}\\ as a
predictor, where the spatial lag is a weighted average of neighboring
observations constructed from an inverse-distance weight matrix. The
autoregressive parameter \\\phi\\ controls the strength and direction of
spatial spillover effects, making this formulation appropriate for
modeling contagion, diffusion, or network effects across spatial units.
The spatial autoregressive model was introduced by Whittle
[\[113\]](#ref113) and developed for areal data by Ord [\[91\]](#ref91).
SAR models are applied in regional economics to estimate how
unemployment rates in neighbouring counties influence each other through
labour market spillovers.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta +
\phi \textbf{z}\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\phi \sim \mathcal{U}(-1, 1)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The SAR model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatially lagged predictor \\\textbf{z}\\ is constructed
from the weight matrix and the response variable itself, creating a
simultaneous dependence structure that requires imperative matrix
operations outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42188)
N <- 100
J <- 3

#### True parameter values
true.beta <- c(2.0, 1.5, -0.8)
true.phi <- 0.4
true.sigma <- 1.0

#### Simulate spatial coordinates
latitude <- runif(N, 0, 100)
longitude <- runif(N, 0, 100)

#### Build weight matrix from inverse distance
D <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))
W <- exp(-D)
W <- ifelse(D == 0, 0, W)

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
epsilon <- rnorm(N, 0, true.sigma)
y <- mu.true + epsilon

#### Compute spatial lag
Z <- W / matrix(rowSums(W), N, N) * matrix(y, N, N, byrow=TRUE)
z <- rowSums(Z)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), phi=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    phi <- runif(1,-1,1)
    sigma <- runif(1)
    return(c(beta, phi, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, latitude=latitude, longitude=longitude,
    mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta,
    pos.phi=pos.phi, pos.sigma=pos.sigma, y=y, z=z)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], -1, 1)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    phi.prior <- dunif(phi, -1, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta)) + phi*Data$z
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + phi.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 0.5, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:", true.beta, "\n")
cat("Post. mean:", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True phi: ", true.phi, "\n")
cat("Post. mean:", round(fit$Summary2[pos.phi, "Mean"], 3), "\n")
cat("True sigma:", true.sigma, "\n")
cat("Post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Spatial Durbin Error Model **\[NEW\]**

The spatial Durbin error model (SDEM) augments a linear regression with
two forms of spatial dependence: exogenous spatial lag effects on the
covariates through \\\textbf{W}\textbf{X}\gamma\\, and endogenous
spatial autocorrelation in the disturbances through the parameter
\\\lambda\\. The reduced-form error covariance \\\sigma^2\[(\textbf{I} -
\lambda\textbf{W})'(\textbf{I} - \lambda\textbf{W})\]^{-1}\\ produces
heteroskedastic, spatially correlated residuals even after accounting
for the exogenous spillovers. The spatial autoregressive error structure
originates with Ord [\[91\]](#ref91), while the Durbin extension that
adds spatially lagged covariates was formalised by LeSage and Pace
[\[162\]](#ref162). The SDEM is applied in regional economics to model
housing prices where both neighbourhood spillovers and spatially
correlated unobserved factors influence outcomes.

#### Form

\\\textbf{y} = \textbf{X}\beta + \textbf{W}\textbf{X}\gamma +
\textbf{u}\\ \\\textbf{u} = \lambda \textbf{W}\textbf{u} + \epsilon\\
\\\epsilon \sim \mathcal{N}(\textbf{0}, \sigma^2 \textbf{I})\\ \\\beta_j
\sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\gamma_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\lambda \sim \mathcal{U}(-1,
1)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The spatial error autocorrelation structure with simultaneous covariate
spatial lags falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42230)
N <- 50
J <- 2

#### True parameter values
true.beta <- c(5, 2)
true.gamma <- c(0, -1)
true.lambda <- 0.5
true.sigma <- 1.0

#### Simulate spatial coordinates on a grid
coords <- cbind(rep(1:ceiling(sqrt(N)), each=ceiling(sqrt(N))),
    rep(1:ceiling(sqrt(N)), ceiling(sqrt(N))))[1:N, ]
D <- as.matrix(dist(coords))

#### Build k-nearest-neighbor weight matrix (k=4, row-normalised)
k <- 4
W <- matrix(0, N, N)
for (i in seq_len(N)) {
    neighbors <- order(D[i, ])[2:(k + 1)]
    W[i, neighbors] <- 1
}
W <- W / rowSums(W)

#### Design matrix and spatial lag of X
X <- cbind(1, rnorm(N))
WX <- W %*% X

#### Generate spatially correlated errors
A <- diag(N) - true.lambda * W
epsilon <- rnorm(N, 0, true.sigma)
u <- solve(A, epsilon)

#### Response
mu.true <- X %*% true.beta + WX %*% true.gamma
y <- as.vector(mu.true + u)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), gamma=rep(0,J),
    lambda=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    gamma <- rnorm(Data$J)
    lambda <- runif(1, -1, 1)
    sigma <- runif(1)
    return(c(beta, gamma, lambda, sigma))
    }
Data <- list(N=N, J=J, PGF=PGF, X=X, W=W, WX=WX,
    mon.names=mon.names, parm.names=parm.names,
    pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.lambda=pos.lambda, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    lambda <- interval(parm[Data$pos.lambda], -0.999, 0.999)
    parm[Data$pos.lambda] <- lambda
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    gamma.prior <- sum(dnormv(gamma, 0, 1000, log=TRUE))
    lambda.prior <- dunif(lambda, -1, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    A <- diag(Data$N) - lambda * Data$W
    mu <- Data$X %*% beta + Data$WX %*% gamma
    e <- A %*% (Data$y - mu)
    log.det.A <- as.numeric(determinant(A, logarithm=TRUE)$modulus)
    LL <- -Data$N / 2 * log(sigma^2) + log.det.A -
        sum(e^2) / (2 * sigma^2)
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + lambda.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, J), rep(0, J), 0, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:    ", true.beta, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True gamma:   ", true.gamma, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.gamma, "Mean"], 3), "\n")
cat("True lambda:  ", true.lambda, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")
cat("True sigma:   ", true.sigma, "\n")
cat("Post. mean:   ", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma]  <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### STARMA(p,q)

The space-time autoregressive moving average (STARMA) model generalizes
classical ARMA to spatiotemporal data by incorporating spatial lags of
both the response and the residuals through adjacency arrays constructed
from \\K\\ nearest neighbors. The autoregressive parameters
\\\phi\_{l,p}\\ capture how past values at neighboring locations
influence the current response, while the moving average parameters
\\\theta\_{m,q}\\ smooth residual spatial dependence, making this a
flexible framework for regularly observed spatiotemporal processes. The
space-time ARMA model was proposed by Pfeifer and Deutsch
[\[93\]](#ref93) to extend ARMA dynamics to spatially referenced time
series. STARMA models are used in traffic engineering to forecast
vehicle flow at multiple road sensors that exhibit both temporal
persistence and spatial propagation.

#### Form

\\\textbf{Y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu\_{s,t} =
\sum^J\_{j=1} \textbf{X}\_{s,t,j} \beta_j + \sum^L\_{l=1} \sum^P\_{p=1}
\phi\_{l,p} \textbf{W}^1\_{s,t-p,l} + \sum^M\_{m=1} \sum^Q\_{q=1}
\theta\_{m,q} \textbf{W}^2\_{s,t-q,m}, \quad j=1,\dots,J, \quad
s=1,\dots,S, \quad t=p,\dots,T\\ \\\textbf{W}^1\_{1:S,1:T,l} =
\textbf{V}\_{1:S,1:S,l} \textbf{Y}, \quad l=1,\dots,L\\
\\\textbf{W}^2\_{1:S,1:T,m} = \textbf{V}\_{1:S,1:S,m} \epsilon, \quad
m=1,\dots,M\\ \\\epsilon = \textbf{Y} - \mu\\ \\\beta_j \sim
\mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\phi\_{l,p} \sim
\mathcal{U}(-1, 1), \quad l=1,\dots,L, \quad p=1,\dots,P\\ \\\sigma \sim
\mathcal{HC}(25)\\ \\\theta\_{m,q} \sim \mathcal{N}(0, 1000), \quad
m=1,\dots,M, \quad q=1,\dots,Q\\ where **V** is an adjacency array that
is scaled so that each row sums to one, \\\beta\\ is a vector of
regression effects, \\\phi\\ is a matrixr of autoregressive space-time
parameters, \\\sigma\\ is the residual variance, and \\\theta\\ is a
matrix of moving average space-time parameters.

#### model_spec() notation

The STARMA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the spatiotemporal lag structure requires constructing
adjacency arrays from nearest neighbors, computing spatial lags of both
the response and residuals, and iterating over autoregressive and moving
average orders, all of which involve imperative loop constructs that
fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42192)

#### Simulation dimensions
S <- 50   # Number of sites
T.obs <- 12  # Number of observed time-periods
J <- 3    # Number of predictors including intercept

#### True parameter values
true.beta <- c(0.5, 0.3, -0.2)
true.phi <- matrix(c(0.3, 0.2, 0.15, 0.1), 2, 2)  # L=2, P=2
true.sigma <- 1.0

#### Simulate spatial coordinates
longitude <- runif(S, -100, -95)
latitude <- runif(S, 30, 35)
#### Build adjacency arrays from K nearest neighbors
K <- 5; L <- 2; M <- 2; P <- 2; Q <- 2
D <- as.matrix(dist(cbind(longitude, latitude), diag=TRUE, upper=TRUE))
A <- V <- array(0, dim=c(S, S, max(L, M)))

#### Generate spatiotemporal data
X <- array(1, dim=c(S, T.obs, J))
X[, , 2] <- matrix(rnorm(S * T.obs), S, T.obs)
X[, , 3] <- matrix(rnorm(S * T.obs), S, T.obs)

Y <- matrix(rnorm(S * T.obs, 0, true.sigma), S, T.obs)
for (l in 1:max(L, M)) {
    A[, , l] <- exp(-D)
    A[, , l] <- apply(A[, , l], 1, rank)
    A[, , l] <- ifelse(A[, , l] > (l-1)*K & A[, , l] <= l*K, 1, 0)
    V[, , l] <- A[, , l] / rowSums(A[, , l])
    V[, , l] <- ifelse(is.nan(V[, , l]), 1/ncol(V[, , l]), V[, , l])
}
for (t in (max(P,Q)+1):T.obs) {
    mu.t <- true.beta[1]*X[,t,1] + true.beta[2]*X[,t,2] + true.beta[3]*X[,t,3]
    for (l in 1:L) { for (p in 1:P) {
        mu.t <- mu.t + true.phi[l,p] * V[,,l] %*% Y[,t-p]
    }}
    Y[, t] <- rnorm(S, mu.t, true.sigma)
}
T <- T.obs

#### Compute spatial lag of Y
W1 <- array(0, dim=c(S, T, max(L, M)))
for (l in 1:max(L, M)) {
    W1[, , l] <- tcrossprod(V[, , l], t(Y))
}
S <- nrow(Y); T <- ncol(Y); J <- dim(X)[3]
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), phi=matrix(0,L,P), sigma=0,
    theta=matrix(0,M,Q)))
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.theta <- grep("theta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    phi <- runif(Data$L*Data$P,-1,1)
    sigma <- runif(1)
    theta <- rnorm(Data$M*Data$Q)
    return(c(beta, phi, sigma, theta))
    }
Data <- list(J=J, K=K, L=L, M=M, P=P, Q=Q, PGF=PGF, S=S, T=T, V=V, W1=W1,
    X=X, Y=Y, latitude=latitude, longitude=longitude, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.phi=pos.phi,
    pos.sigma=pos.sigma, pos.theta=pos.theta)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    phi <- matrix(interval(parm[Data$pos.phi], -1, 1), Data$L, Data$P)
    parm[Data$pos.phi] <- as.vector(phi)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    theta <- matrix(parm[Data$pos.theta], Data$M, Data$Q)
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    phi.prior <- sum(dunif(phi, -1, 1, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    theta.prior <- sum(dnormv(theta, 0, 1000, log=TRUE))
    ### Log-Likelihood
    mu <- beta[1]*Data$X[, , 1]
    for (j in 2:Data$J) mu <- mu + beta[j]*Data$X[, , j]
    for (l in 1:Data$L) {for (p in 1:Data$P) {
        mu[,-c(1:p)] <- mu[,-c(1:p)] +
    phi[l,p]*Data$W1[, 1:(Data$T - p), l]}}
    epsilon <- Data$Y - mu
    for (m in 1:Data$M) {
        W2 <- tcrossprod(Data$V[, , m], t(epsilon))
        for (q in 1:Data$Q) {
    mu[,-c(1:q)] <- mu[,-c(1:q)] +
    theta[m,q]*W2[,1:(Data$T - q)]}}
    LL <- sum(dnorm(Data$Y[,-c(1:max(Data$P,Data$Q))],
        mu[,-c(1:max(Data$P,Data$Q))], sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + phi.prior + sigma.prior + theta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(0,L*P), 1, rep(0,M*Q))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.phi]   <- true.phi[seq_along(pos.phi)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Factor analysis and latent variable models

### Factor Analysis, Approximate Dynamic

The Approximate Dynamic Factor Analysis (ADFA) model has many names,
including the approximate factor model and approximate dynamic factor
model. An ADFA is a Dynamic Factor Analysis (DFA) in which the factor
scores of the dynamic factors are approximated with principal
components. This is a combination of principal components and common
factor analysis, in which the factor loadings of common factors are
estimated from the data and factor scores are estimated from principal
components. This is a two-stage model: principal components are
estimated in the first stage and a decision is made regarding how many
principal components to retain, and ADFA is estimated in the second
stage. For more information on DFA, see section [link](#dsfm).
Approximate dynamic factor analysis provides a computationally efficient
alternative to full dynamic factor models, following the approach of
Stock and Watson [\[101\]](#ref101). Dynamic factor analysis is used in
neuroscience to extract latent neural population dynamics from
multi-electrode array recordings.

#### Form

\\\textbf{Y}\_{t,j} \sim \mathcal{N}(\mu\_{t,j}, \sigma^2_j), \quad
t=2,\dots,T, \quad j=1,\dots,J\\ \\\mu\_{t,} = \textbf{F}\_{t-1,}
\Lambda\\ \\\Lambda\_{p,j} \sim \mathcal{N}(0, 1), \quad p=1,\dots,P,
\quad j=1,\dots,J\\ \\\sigma_j \sim \mathcal{HC}(25), \quad
j=1,\dots,J\\

#### model_spec() notation

The ADFA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the factor scores are pre-computed via principal components
analysis in a separate first stage, and the second-stage likelihood
involves lagged factor score matrices, requiring a two-stage estimation
procedure that falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42180)

#### Simulation dimensions
T <- 99   # Number of time-periods
J <- 10   # Number of time-series
P <- 3    # Number of approximate factors

#### True parameter values
true.Lambda <- matrix(rnorm(P * J, 0, 0.5), P, J)
true.sigma <- rep(0.5, J)

#### Generate factor scores as AR(1) processes
F.true <- matrix(0, T + 1, P)
for (t in 2:(T + 1)) {
    F.true[t, ] <- 0.8 * F.true[t - 1, ] + rnorm(P, 0, 0.3)
}
F.true <- F.true[-1, ]

#### Generate observations
mu.true <- tcrossprod(rbind(rep(0, P), F.true[-T, ]), t(true.Lambda))
Y <- matrix(rnorm(T * J, mu.true, rep(true.sigma, each = T)), T, J)

#### First stage: principal components
PCA <- prcomp(Y, scale=TRUE)
F <- PCA$x[, 1:P]
mon.names <- "LP"
parm.names <- as.parm.names(list(Lambda=matrix(0,P,J), sigma=rep(0,J)))
pos.Lambda <- grep("Lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    Lambda <- rnorm(Data$P*Data$J)
    sigma <- runif(Data$J)
    return(c(Lambda, sigma))
    }
Data <- list(F=F, J=J, P=P, PGF=PGF, T=T, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.Lambda=pos.Lambda, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    Lambda <- matrix(parm[Data$pos.Lambda], Data$P, Data$J)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    Lambda.prior <- sum(dnorm(Lambda, 0, 1, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(rbind(rep(0,Data$P), F[-Data$T,]), t(Lambda))
    Sigma <- matrix(sigma, Data$T, Data$J, byrow=TRUE)
    LL <- sum(dnorm(Data$Y[-1,], mu[-1,], Sigma[-1,], log=TRUE))
    ### Log-Posterior
    LP <- LL + Lambda.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,P*J), rep(1,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.Lambda] <- true.Lambda[seq_along(pos.Lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Factor Analysis, Confirmatory

Confirmatory factor analysis (CFA) tests a pre-specified measurement
model in which each observed variable loads on exactly one latent factor
as indicated by the assignment vector \\\textbf{f}\\, with factor scores
in matrix \\\textbf{F}\\ and loadings in \\\lambda\\. The Wishart prior
on the factor precision matrix \\\Omega\\ allows the latent factors to
correlate, making this model suitable for testing hypothesized latent
structures in psychometric and survey data. Confirmatory factor analysis
was formalised by Jöreskog [\[61\]](#ref61) and placed in a Bayesian
framework by Muthén and Asparouhov [\[88\]](#ref88). CFA is standard in
psychometrics for validating that a questionnaire’s items load on
hypothesised latent constructs such as anxiety and depression.

#### Form

\\\textbf{Y}\_{i,m} \sim \mathcal{N}(\mu\_{i,m}, \sigma^2_m), \quad
i=1,\dots,N, \quad m=1,\dots,M\\ \\\mu = \textbf{F}\_{1:N,\textbf{f}}
\lambda^T\\ \\\textbf{F}\_{i,1:P} \sim \mathcal{N}\_P(0, \Omega^{-1}),
\quad i=1,\dots,N\\ \\\lambda_m \sim \mathcal{N}(0, 1), \quad
m=1,\dots,M\\ \\\sigma_m \sim \mathcal{HC}(25), \quad m=1,\dots,M\\
\\\Omega \sim \mathcal{W}\_N(\textbf{S}), \quad \textbf{S} =
\textbf{I}\_P\\

#### model_spec() notation

The CFA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the factor structure requires a Wishart prior on the
precision matrix, multivariate normal factor scores with a precision
parameterization, and indicator-based indexing of factors to variables,
all of which involve imperative matrix operations that fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42181)

#### Simulation dimensions
N <- 50   # Number of records
M <- 5    # Number of variables
P <- 3    # Number of factors
f <- c(1, 3, 2, 2, 1)  # Factor assignment for each variable

#### True parameter values
true.lambda <- c(0.8, 0.7, 0.9, 0.6, 0.75)
true.sigma <- rep(0.5, M)
true.Omega <- diag(P)

#### Generate factor scores and observations
F.true <- rmvnpc(N, rep(0, P), true.Omega)
mu.true <- F.true[, f] * matrix(true.lambda, N, M, byrow=TRUE)
Y <- matrix(rnorm(N * M, mu.true, rep(true.sigma, each = N)), N, M)
for (m in 1:M) Y[, m] <- CenterScale(Y[, m])
S <- diag(P)
mon.names <- "LP"
parm.names <- as.parm.names(list(F=matrix(0,N,P), lambda=rep(0,M),
    U=diag(P), sigma=rep(0,M)), uppertri=c(0,0,1,0))
pos.F <- grep("F", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    U <- rwishartc(Data$N, Data$S)
    F <- as.vector(rmvnpc(Data$N, rep(0,Data$P), U))
    U <- U[upper.tri(U, diag=TRUE)]
    lambda <- rnorm(Data$M)
    sigma <- runif(Data$M)
    return(c(F, lambda, U, sigma))
    }
Data <- list(M=M, N=N, P=P, PGF=PGF, S=S, Y=Y, f=f, mon.names=mon.names,
    parm.names=parm.names, pos.F=pos.F, pos.lambda=pos.lambda,
    pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    lambda <- parm[Data$pos.lambda]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    F <- matrix(parm[Data$pos.F], Data$N, Data$P)
    U <- as.parm.matrix(U, Data$P, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    ### Log-Prior
    lambda.prior <- sum(dnorm(lambda, 0, 1, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    U.prior <- dwishartc(U, Data$N, Data$S, log=TRUE)
    F.prior <- sum(dmvnpc(F, rep(0,Data$P), U, log=TRUE))
    ### Log-Likelihood
    mu <- F[,Data$f] * matrix(lambda, Data$N, Data$M, byrow=TRUE)
    Sigma <- matrix(sigma, Data$N, Data$M, byrow=TRUE)
    LL <- sum(dnorm(Data$Y, mu, Sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + lambda.prior + sigma.prior + F.prior + U.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,N*P), rep(0,M), upper.triangle(S, diag=TRUE),
    rep(1,M))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True lambda:", true.lambda, "\n")
cat("Post. mean: ", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
ground_truth[grep("Omega", parm.names)]  <- true.Omega[seq_along(grep("Omega", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Factor Analysis, Exploratory

Exploratory factor analysis (EFA) estimates both the factor scores
\\\textbf{F}\\ and the loading matrix \\\Lambda\\ simultaneously, with
an upper-triangular constraint on \\\Lambda\\ to resolve rotational
indeterminacy. The Wishart prior on the factor precision matrix allows
correlated factors, while the half-Cauchy priors on the residual
standard deviations provide weakly informative regularization.
Exploratory factor analysis originated with Spearman [\[100\]](#ref100)
and was developed by Thurstone [\[106\]](#ref106); the Bayesian
treatment follows Lopes and West [\[76\]](#ref76). EFA is applied in
personality psychology to discover the latent dimensions underlying
responses to large item batteries.

#### Form

\\\textbf{Y}\_{i,m} \sim \mathcal{N}(\mu\_{i,m}, \sigma^2_m), \quad
i=1,\dots,N, \quad m=1,\dots,M\\ \\\mu = \textbf{F} \Lambda\\
\\\textbf{F}\_{i,1:P} \sim \mathcal{N}\_P(0, \Omega^{-1}), \quad
i=1,\dots,N\\ \\\Lambda\_{p,m} \sim \mathcal{N}(0, 1), \quad
p=1,\dots,P, \quad m=(p+1),\dots,M\\ \\\Omega \sim
\mathcal{W}\_N(\textbf{S}), \quad \textbf{S} = \textbf{I}\_P\\
\\\sigma_m \sim \mathcal{HC}(25), \quad m=1,\dots,M\\

#### model_spec() notation

The EFA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the upper-triangular loading matrix constraint, the Wishart
prior on the factor precision, and the multivariate normal factor scores
with precision parameterization require imperative matrix constructions
that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42182)

#### Simulation dimensions
N <- 43   # Number of records
M <- 12   # Number of variables
P <- 3    # Number of factors

#### True parameter values (upper-triangular Lambda)
Lambda.true <- matrix(1, P, M)
Lambda.true[lower.tri(Lambda.true)] <- 0
Lambda.true[upper.tri(Lambda.true)] <- rnorm(P*M - P - P*(P-1)/2, 0, 0.5)
true.sigma <- rep(0.5, M)

#### Generate factor scores and observations
S <- diag(P)
F.true <- rmvnpc(N, rep(0, P), diag(P))
mu.true <- tcrossprod(F.true, t(Lambda.true))
Y <- matrix(rnorm(N * M, mu.true, rep(true.sigma, each = N)), N, M)
for (m in 1:M) Y[, m] <- CenterScale(Y[, m])
Lambda <- matrix(NA, P, M)
Lambda[upper.tri(Lambda)] <- 0
mon.names <- "LP"
parm.names <- as.parm.names(list(F=matrix(0,N,P), Lambda=Lambda, U=S,
    sigma=rep(0,M)), uppertri=c(0,0,1,0))
pos.F <- grep("F", parm.names)
pos.Lambda <- grep("Lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    U <- rwishartc(Data$N, Data$S)
    F <- as.vector(rmvnpc(Data$N, rep(0,Data$P), U))
    Lambda <- rnorm(Data$P*Data$M-Data$P-Data$P*(Data$P-1)/2,0,1)
    sigma <- runif(Data$M)
    return(c(F, Lambda, U[upper.tri(U, diag=TRUE)], sigma))
    }
Data <- list(M=M, N=N, P=P, PGF=PGF, S=S, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.F=pos.F, pos.Lambda=pos.Lambda,
    pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    F <- matrix(parm[Data$pos.F], Data$N, Data$P)
    lambda <- parm[Data$pos.Lambda]
    U <- as.parm.matrix(U, Data$P, parm, Data, chol=TRUE)
    diag(U) <- exp(diag(U))
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    F.prior <- sum(dmvnpc(F, rep(0,Data$P), U, log=TRUE))
    Lambda.prior <- sum(dnorm(lambda, 0, 1, log=TRUE))
    U.prior <- dwishartc(U, Data$N, Data$S, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    Lambda <- matrix(1, Data$P, Data$M)
    Lambda[lower.tri(Lambda)] <- 0
    Lambda[upper.tri(Lambda)] <- lambda
        mu <- tcrossprod(F, t(Lambda))
    Sigma <- matrix(sigma, Data$N, Data$M, byrow=TRUE)
    LL <- sum(dnorm(Data$Y, mu, Sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + F.prior + Lambda.prior + U.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
    yhat=rnorm(prod(dim(mu)), mu, Sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,N*P), rep(0,P*M-P-P*(P-1)/2), rep(0,P*(P-1)/2+P),
    rep(1,M))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Factor Analysis, Exploratory Ordinal

The exploratory ordinal factor analysis (EOFA) model extends EFA to
ordinal or categorical outcomes through a probit link, estimating
cumulative probability thresholds \\\alpha\\ and factor loadings
\\\Lambda\\ simultaneously. This formulation is also suitable for
collaborative filtering, where users rate items on a discrete scale and
the latent factors capture preference dimensions. Exploratory ordinal
factor analysis extends EFA to categorical indicators using threshold
models, following Muthén [\[87\]](#ref87) and Albert and Chib
[\[2\]](#ref2). Ordinal factor analysis is used in educational testing
to extract latent ability dimensions from Likert-scale survey items.

#### Form

\\\textbf{Y}\_{i,m} \sim \mathcal{CAT}(\textbf{P}\_{i,m,1:K}), \quad
i=1,\dots,N, \quad m=1,\dots,M\\ \\\textbf{P}\_{,,K} = 1 -
Q\_{,,(K-1)}\\ \\\textbf{P}\_{,,k} = \|Q\_{,,k} - Q\_{,,(k-1)}\|, \quad
k=2,\dots,(K-1)\\ \\\textbf{P}\_{,,1} = Q\_{,,1}\\ \\Q = \phi(\mu)\\
\\\mu\_{,,k} = \alpha_k - \textbf{F} \Lambda, \quad k=1,\dots,(K-1)\\
\\\textbf{F}\_{i,1:P} \sim \mathcal{N}\_P(0, \Omega^{-1}), \quad
i=1,\dots,N\\ \\\gamma_p = 0, \quad p=1,\dots,P\\ \\\Lambda\_{p,m} \sim
\mathcal{N}(0, 1), \quad p=1,\dots,P, \quad m=(p+1),\dots,M\\ \\\Omega
\sim \mathcal{W}\_N(\textbf{S}), \quad \textbf{S} = \textbf{I}\_P\\
\\\alpha_k \sim \mathcal{N}(0, 1) \in \[(k-1),k\] \in \[-5,5\], \quad
k=2,\dots,(K-1)\\

#### model_spec() notation

The EOFA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the ordinal probit link requires computing cumulative
probabilities from sorted thresholds, constructing a three-dimensional
probability array, and evaluating the categorical likelihood, all of
which involve imperative array operations that fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42183)

#### Simulation dimensions
M <- 10  # Number of variables
N <- 20  # Number of records
K <- 3   # Number of discrete values
P <- 3   # Number of factors

#### True parameter values
true.alpha <- sort(c(-0.5, 0.5))
true.Lambda <- matrix(1, P, M)
true.Lambda[lower.tri(true.Lambda)] <- 0
true.Lambda[upper.tri(true.Lambda)] <- rnorm(P*M-P-P*(P-1)/2, 0, 0.5)
true.Omega <- runif(P, 0.5, 1.5)

#### Generate factor scores and ordinal responses
F.true <- rmvnp(N, rep(0, P), true.Omega)
mu <- aperm(array(true.alpha, dim=c(K-1, M, N)), perm=c(3,2,1))
mu <- mu - array(tcrossprod(F.true, t(true.Lambda)), dim=c(N, M, K-1))
Pr <- Q <- pnorm(mu)
Pr[, , -1] <- abs(Q[, , -1] - Q[, , -(K-1)])
Pr <- array(Pr, dim=c(N, M, K))
Pr[, , K] <- 1 - Q[, , (K-1)]
dim(Pr) <- c(N*M, K)
Y <- matrix(rcat(nrow(Pr), Pr), N, M)
S <- diag(P)
Lambda <- matrix(0, P, M)
Lambda[lower.tri(Lambda, diag=TRUE)] <- NA
mon.names <- "LP"
parm.names <- as.parm.names(list(F=matrix(0,N,P), Omega=rep(0,P),
    Lambda=Lambda, alpha=rep(0,K-1)))
pos.F <- grep("F", parm.names)
pos.Omega <- grep("Omega", parm.names)
pos.Lambda <- grep("Lambda", parm.names)
pos.alpha <- grep("alpha", parm.names)
PGF <- function(Data) {
    Omega <- runif(Data$P)
    F <- as.vector(rmvnp(Data$N, rep(0,Data$P), diag(Omega)))
    Lambda <- rnorm(Data$P*Data$M-Data$P-Data$P*(Data$P-1)/2)
    alpha <- sort(rnorm(Data$K-1))
    return(c(F, Omega, Lambda, alpha))
    }
Data <- list(K=K, M=M, N=N, P=P, PGF=PGF, S=S, Y=Y,
    mon.names=mon.names, parm.names=parm.names, pos.F=pos.F,
    pos.Omega=pos.Omega, pos.Lambda=pos.Lambda, pos.alpha=pos.alpha)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    F <- matrix(parm[Data$pos.F], Data$N, Data$P)
    Omega <- interval(parm[Data$pos.Omega], 1e-100, Inf)
    parm[Data$pos.Omega] <- Omega
    lambda <- parm[Data$pos.Lambda]
    alpha <- sort(interval(parm[Data$pos.alpha], -5, 5))
    parm[Data$pos.alpha] <- alpha
    ### Log-Prior
    F.prior <- sum(dmvnp(F, rep(0,Data$P), diag(Omega), log=TRUE))
    Omega.prior <- dwishart(diag(Omega), Data$N, Data$S, log=TRUE)
    Lambda.prior <- sum(dnorm(lambda, 0, 1, log=TRUE))
    alpha.prior <- sum(dnormv(alpha, 0, 10, log=TRUE))
    ### Log-Likelihood
    Lambda <- matrix(1, Data$P, Data$M)
    Lambda[lower.tri(Lambda)] <- 0
    Lambda[upper.tri(Lambda)] <- lambda
    mu <- aperm(array(alpha, dim=c(Data$K-1, Data$M, Data$N)),
        perm=c(3,2,1))
    mu <- mu - array(tcrossprod(F, t(Lambda)),
        dim=c(Data$N, Data$M, Data$K-1))
    P <- Q <- pnorm(mu)
    P[ , , -1] <- abs(Q[ , , -1] - Q[ , , -(Data$K-1)])
    P <- array(P, dim=c(Data$N, Data$M, Data$K))
    P[ , , Data$K] <- 1 - Q[ , , (Data$K-1)]
    y <- as.vector(Data$Y)
    dim(P) <- c(Data$N*Data$M, Data$K)
    LL <- sum(dcat(y, P, log=TRUE))
    ### Log-Posterior
    LP <- LL + F.prior + Omega.prior + Lambda.prior + alpha.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=matrix(rcat(nrow(P), P), Data$N, Data$M), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,N*P), rep(0,P), rep(0,P*M-P-P*(P-1)/2),
    seq(from=-1, to=1, len=K-1))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True alpha:", true.alpha, "\n")
cat("Post. mean:", round(fit$Summary2[pos.alpha, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.Lambda] <- true.Lambda[seq_along(pos.Lambda)]
ground_truth[pos.Omega]  <- true.Omega[seq_along(pos.Omega)]
ground_truth[pos.alpha]  <- true.alpha[seq_along(pos.alpha)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Factor Regression

Factor regression orthogonalizes the independent variables through
latent factors before regressing the response on those factors,
combining confirmatory factor analysis (section [link](#cfa)) with
linear regression (section [link](#linear.reg)). This formulation is
constrained to the case where the number of factors equals the number of
predictors less the intercept, aiming to remove collinearity rather than
reduce dimensionality. Factor regression combines factor analysis with a
regression response model, following West [\[111\]](#ref111) and
Carvalho et al. [\[20\]](#ref20). Factor regression is applied in
financial economics to explain asset returns through a small number of
latent factors extracted from a panel of macroeconomic indicators.

#### Form

\\\textbf{y} \sim \mathcal{N}(\nu, \sigma^2\_{J+1})\\ \\\nu = \textbf{F}
\beta\\ \\\mu\_{i,1} = 1\\ \\\mu\_{i,j+1} = \mu\_{i,j}, \quad
j=1,\dots,J\\ \\\textbf{X}\_{i,j} \sim \mathcal{N}(\mu\_{i,j},
\sigma^2_j), \quad i=1,\dots,N, \quad j=2,\dots,J\\ \\\mu\_{i,j} =
\lambda_j \textbf{F}\_{i,j}, \quad i=1,\dots,N, \quad j=2,\dots,J\\
\\\textbf{F}\_{i,1:J} \sim \mathcal{N}\_{J-1}(0, \Omega^{-1}), \quad
i=1,\dots,N\\ \\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\
\\\lambda_j \sim \mathcal{N}(0, 1), \quad j=1,\dots,(J-1)\\ \\\sigma_j
\sim \mathcal{HC}(25), \quad j=1,\dots,(J+1)\\ \\\Omega \sim
\mathcal{W}\_N(\textbf{S}), \quad \textbf{S} = \textbf{I}\_J\\

#### model_spec() notation

The factor regression model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a Wishart prior on the factor precision matrix,
multivariate normal factor scores with precision parameterization, and a
two-part likelihood coupling factor analysis of the predictors with
regression on the response, all of which involve imperative matrix
operations that fall outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42184)

#### Simulation dimensions
N <- 30
J <- 3

#### True parameter values
true.beta <- c(2.0, 0.5, -0.3, 0.8)  # intercept + J factor effects
true.lambda <- c(0.9, 0.7, 0.8)
true.sigma <- c(0.5, 0.5, 0.5, 1.0)  # J for X + 1 for y
true.Omega <- rep(1, J)

#### Generate factor scores
F.true <- rmvnp(N, rep(0, J), diag(true.Omega))

#### Generate predictors from factor model
X <- F.true * matrix(true.lambda, N, J, byrow=TRUE) +
    matrix(rnorm(N * J, 0, rep(true.sigma[1:J], each = N)), N, J)
for (j in 1:J) X[, j] <- CenterScale(X[, j])

#### Generate response from factor regression
nu.true <- cbind(1, F.true) %*% true.beta
y <- rnorm(N, nu.true, true.sigma[J + 1])
S <- diag(J)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J+1), lambda=rep(0,J),
    sigma=rep(0,J+1), F=matrix(0,N,J), Omega=rep(0,J)))
pos.beta <- grep("beta", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.F <- grep("F", parm.names)
pos.Omega <- grep("Omega", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J+1)
    lambda <- rnorm(Data$J)
    sigma <- runif(Data$J+1)
    Omega <- runif(Data$J)
    F <- as.vector(rmvnp(Data$N, rep(0,Data$J), diag(Omega)))
    return(c(beta, lambda, sigma, F, Omega))
    }
Data <- list(J=J, N=N, PGF=PGF, S=S, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.lambda=pos.lambda,
    pos.sigma=pos.sigma, pos.F=pos.F, pos.Omega=pos.Omega, y=y)
```

#### Model

``` r
Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    lambda <- parm[Data$pos.lambda]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    F <- matrix(Data$pos.F], Data$N, Data$J)
    Omega <- interval(parm[Data$pos.Omega], 1e-100, Inf)
    parm[Data$pos.Omega] <- Omega
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    lambda.prior <- sum(dnorm(lambda, 0, 1, log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    F.prior <- sum(dmvnp(F, rep(0,Data$J), diag(Omega), log=TRUE))
    Omega.prior <- dwishart(diag(Omega), Data$N, Data$S, log=TRUE)
    ### Log-Likelihood
    mu <- F * matrix(lambda, Data$N, Data$J, byrow=TRUE)
    nu <- tcrossprod(cbind(1,F), t(beta))
    LL <- sum(dnorm(Data$X, mu, matrix(sigma[1:Data$J], Data$N, Data$J,
        byrow=TRUE), log=TRUE))
    LL <- LL + dnorm(Data$y, nu, sigma[Data$J+1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + lambda.prior + sigma.prior + F.prior
        Omega.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(Data$N, nu, sigma[Data$J+1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J+1), rep(0,J), rep(0,J+1), rep(0,N*J), rep(1,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
cat("True beta:  ", true.beta, "\n")
cat("Post. mean: ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")
cat("True lambda:", true.lambda, "\n")
cat("Post. mean: ", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.Omega]  <- true.Omega[seq_along(pos.Omega)]
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### IRT Two-Parameter Logistic (2PL) **\[NEW\]**

The two-parameter logistic (2PL) item response theory model, developed
by Birnbaum [\[153\]](#ref153), is the workhorse of educational and
psychological measurement. Each binary response \\y\_{ij}\\ is modelled
as a Bernoulli trial whose success probability depends on the latent
ability \\\theta_i\\ of person \\i\\ and two item characteristics: the
discrimination \\a_j\\ (how sharply the item differentiates between
ability levels) and the difficulty \\b_j\\ (the ability level at which
the probability of a correct response is 0.5). A standard normal prior
on the abilities provides the location and scale identification that the
model requires. Bayesian IRT estimation was formalised by Albert
[\[154\]](#ref154) using data augmentation. The 2PL model is the
standard in large-scale educational testing programmes for calibrating
examination items.

#### Form

\\y\_{ij} \sim \text{Bernoulli}(p\_{ij}), \quad i=1,\dots,I, \\
j=1,\dots,J\\ \\\text{logit}(p\_{ij}) = a_j\\(\theta_i - b_j)\\
\\\theta_i \sim \mathcal{N}(0, 1)\\ \\a_j \sim \text{LogNormal}(0.5,\\
0.5)\\ \\b_j \sim \mathcal{N}(0, 4)\\

#### model_spec() notation

The 2PL model cannot be expressed in
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the double-indexed likelihood with crossed person and item
parameters, the multiplicative interaction \\a_j(\theta_i - b_j)\\, and
the lognormal constraint on discriminations fall outside the declarative
scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42217)
N.persons <- 30
N.items <- 10
N.obs <- N.persons * N.items

#### True parameter values
theta.true <- rnorm(N.persons)
a.true <- exp(rnorm(N.items, 0.5, 0.3))
b.true <- rnorm(N.items, 0, 1)

#### Generate response matrix and flatten
Y <- matrix(NA, N.persons, N.items)
for (i in 1:N.persons) {
    for (j in 1:N.items) {
        p.ij <- plogis(a.true[j] * (theta.true[i] - b.true[j]))
        Y[i,j] <- rbinom(1, 1, p.ij)
    }
}
y <- as.vector(Y)
person <- rep(1:N.persons, N.items)
item <- rep(1:N.items, each = N.persons)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(theta=rep(0, N.persons),
    a=rep(0, N.items), b=rep(0, N.items)))
pos.theta <- grep("theta", parm.names)
pos.a <- grep("^a", parm.names)
pos.b <- grep("^b", parm.names)
PGF <- function(Data) {
    theta <- rnorm(Data$N.persons)
    a <- exp(rnorm(Data$N.items, 0.5, 0.3))
    b <- rnorm(Data$N.items, 0, 1)
    return(c(theta, a, b))
    }
Data <- list(N.persons=N.persons, N.items=N.items, N.obs=N.obs,
    PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    pos.theta=pos.theta, pos.a=pos.a, pos.b=pos.b,
    person=person, item=item, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    theta <- parm[Data$pos.theta]
    a <- interval(parm[Data$pos.a], 1e-100, Inf)
    parm[Data$pos.a] <- a
    b <- parm[Data$pos.b]
    ### Log-Prior
    theta.prior <- sum(dnorm(theta, 0, 1, log=TRUE))
    a.prior <- sum(dlnorm(a, 0.5, 0.5, log=TRUE))
    b.prior <- sum(dnormv(b, 0, 4, log=TRUE))
    ### Log-Likelihood
    p <- plogis(a[Data$item] * (theta[Data$person] - b[Data$item]))
    LL <- sum(dbinom(Data$y, 1, p, log=TRUE))
    ### Log-Posterior
    LP <- LL + theta.prior + a.prior + b.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rbinom(Data$N.obs, 1, p), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, N.persons), rep(1, N.items), rep(0, N.items))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery: item discrimination and difficulty
cat("Item discrimination (a):\n")
cat("  true:      ", round(a.true, 3), "\n")
cat("  post. mean:", round(fit$Summary2[pos.a, "Mean"], 3), "\n")
cat("Item difficulty (b):\n")
cat("  true:      ", round(b.true, 3), "\n")
cat("  post. mean:", round(fit$Summary2[pos.b, "Mean"], 3), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Latent Dirichlet Allocation

Latent Dirichlet Allocation was introduced by Blei et
al. [\[13\]](#ref13) as a generative probabilistic model for discovering
latent topics in text corpora. LDA is widely used in computational
social science to identify thematic structures in political speeches,
news articles, and social media posts. \### Form

\\\textbf{Y}\_{m,n} \sim \mathcal{CAT}(\phi\[\textbf{Z}\_{m,n},\]),
\quad m=1,\dots,M, \quad n=1,\dots,N\\ \\\textbf{Z}\_{m,n} \sim
\mathcal{CAT}(\theta\_{m,1:K})\\ \\\phi\_{k,v} \sim \mathcal{D}(\beta)\\
\\\theta\_{m,k} \sim \mathcal{D}(\alpha)\\ \\\alpha_k = 1, \quad
k=1,\dots,K\\ \\\beta_v = 1, \quad v=1,\dots,V\\

#### Data

``` r

K <- 2 #Number of (latent) topics
M <- 4 #Number of documents in corpus
N <- 15 #Maximum number of (used) words per document
V <- 5 #Maximum number of occurrences of any word (Vocabulary size)
Y <- matrix(rcat(M*N,rep(1/V,V)), M, N)
rownames(Y) <- paste("doc", 1:nrow(Y), sep="")
colnames(Y) <- paste("word", 1:ncol(Y), sep="")
#Note: Y is usually represented as w, a matrix of word counts.
if(min(Y) == 0) Y <- Y + 1 #A zero cannot occur, Y must be 1,2,...,V.
V <- max(Y) #Maximum number of occurrences of any word (Vocabulary size)
alpha <- rep(1,K) # hyperparameters (constant)
beta <- rep(1,V)
mon.names <- "LP"
parm.names <- as.parm.names(list(phi=matrix(0,K,V), theta=matrix(0,M,K),
    Z=matrix(0,M,N)))
pos.phi <- grep("phi", parm.names)
pos.theta <- grep("theta", parm.names)
pos.Z <- grep("Z", parm.names)
PGF <- function(Data) {
    phi <- matrix(runif(Data$J*Data$V), Data$K, Data$V)
    phi <- phi / rowSums(phi)
    theta <- matrix(runif(Data$M*Data$K), Data$M, Data$K)
    theta <- theta / rowSums(theta)
    z <- rcat(Data$M*Data$N, rep(1/Data$K,Data$K))
    return(c(as.vector(phi), as.vector(theta), z))}
Data <- list(K=K, M=M, N=N, PGF=PGF, V=V, Y=Y, alpha=alpha, beta=beta,
    mon.names=mon.names, parm.names=parm.names, pos.phi=pos.phi,
    pos.theta=pos.theta, pos.Z=pos.Z)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    phi <- matrix(interval(parm[Data$pos.phi], 0, 1), Data$K, Data$V)
    phi <- phi / rowSums(phi)
    parm[Data$pos.phi] <- as.vector(phi)
    theta <- matrix(interval(parm[Data$pos.theta], 0, 1), Data$M, Data$K)
    theta <- theta / rowSums(theta)
    parm[Data$pos.theta] <- as.vector(theta)
    Z <- matrix(parm[Data$pos.Z], Data$M, Data$N)
    ### Log-Prior
    phi.prior <- sum(ddirichlet(phi, beta, log=TRUE))
    theta.prior <- sum(ddirichlet(theta, alpha, log=TRUE))
    ### Log-Likelihood
    LL <- Z.prior <- 0
    Yhat <- Data$Y
    for (m in 1:Data$M) {for (n in 1:Data$N) {
        Z.prior + Z.prior + dcat(Z[m,n], theta[m,], log=TRUE)
        LL <- LL + dcat(Data$Y[m,n], as.vector(phi[Z[m,n],]), log=TRUE)
        Yhat[m,n] <- rcat(1, as.vector(phi[Z[m,n],]))}}
    ### Log-Posterior
    LP <- LL + phi.prior + theta.prior + Z.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=Yhat, parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(1/V,K*V), rep(1/K,M*K), rcat(M*N,rep(1/K,K)))
```

#### model_spec() notation

The LDA model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it involves discrete latent topic assignments \\Z\_{m,n}\\
drawn from document-specific Dirichlet-Multinomial distributions, with
the likelihood coupling topic-word distributions \\\phi\\ to observed
words through these latent indices. This discrete mixture structure with
matrix-indexed categorical draws falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Inspect estimated topic-word distributions
phi.post <- matrix(fit$Summary2[pos.phi, "Mean"], K, V)
cat("Estimated topic-word distributions (phi):\n")
print(round(phi.post, 3))

#### Inspect estimated document-topic proportions
theta.post <- matrix(fit$Summary2[pos.theta, "Mean"], M, K)
cat("Estimated document-topic proportions (theta):\n")
print(round(theta.post, 3))

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

  

## Mixture and clustering models

### Cluster Analysis, Confirmatory (CCA)

This is a parametric, model-based, cluster analysis, also called a
finite mixture model or latent class cluster analysis, where the number
of clusters \\C\\ is fixed. When the number of clusters is unknown,
exploratory cluster analysis should be used. The record-level cluster
membership parameter vector, \\\theta\\, is a vector of discrete
parameters. Discrete parameters are not supported in all algorithms. The
example below is updated with the Slice sampler. Bayesian confirmatory
cluster analysis assigns observations to predetermined cluster
structures, following the model-based clustering framework of Banfield
and Raftery [\[4\]](#ref4). CCA is applied in market segmentation to
confirm whether customer purchase patterns match hypothesised consumer
archetypes.

#### Form

\\\textbf{Y}\_{i,j} \sim \mathcal{N}(\mu\_{\theta\[i\],j},
\sigma^2\_{\theta\[i\]}), \quad i=1,\dots,N, \quad j=1,\dots,J\\
\\\theta_i \sim \mathcal{CAT}(\pi\_{1:C}), \quad i=1,\dots,N\\
\\\pi\_{1:C} \sim \mathcal{D}(\alpha\_{1:C})\\ \\\alpha_c = 1\\
\\\mu\_{c,j} \sim \mathcal{N}(0, \nu^2_c), \quad c=1,\dots,C, \quad
j=1,\dots,J\\ \\\sigma_c \sim \mathcal{HC}(25), \quad c=1,\dots,C\\
\\\nu_c \sim \mathcal{HC}(25), \quad c=1,\dots,C\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires discrete cluster allocation parameters
\\\theta_i\\ sampled from a Categorical distribution, with mixing
proportions \\\pi\\ derived from the empirical cluster frequencies and a
Dirichlet prior. The combination of discrete latent allocation,
per-cluster indexing of means and variances, and the
Dirichlet-Categorical hierarchy falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42119)
N <- 100
J <- 4
C <- 3

#### True parameter values
true.pi <- c(0.4, 0.35, 0.25)
true.mu <- matrix(c(-2, 0, 2, 1, -1, 0.5, 0, 1.5, -1, 0.8, -0.5, 1.2),
    C, J)
true.sigma <- c(0.5, 0.8, 0.6)

#### Generate cluster memberships and data
true.theta <- rcat(N, true.pi)
Y <- matrix(0, N, J)
for (i in 1:N)
    Y[i, ] <- rnorm(J, true.mu[true.theta[i], ], true.sigma[true.theta[i]])

#### Assemble Data list
alpha <- rep(1, C)
mon.names <- c("LP", paste("pi[", 1:C, "]", sep=""))
parm.names <- as.parm.names(list(theta=rep(0,N), mu=matrix(0,C,J),
    nu=rep(0,C), sigma=rep(0,C)))
pos.theta <- grep("theta", parm.names)
pos.mu <- grep("mu", parm.names)
pos.nu <- grep("nu", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    theta <- rcat(Data$N, p=rep(1/Data$C, Data$C))
    mu <- rnorm(Data$C*Data$J)
    nu <- runif(Data$C)
    sigma <- runif(Data$C)
    return(c(theta, mu, nu, sigma))
    }
Data <- list(C=C, J=J, N=N, PGF=PGF, Y=Y, alpha=alpha,
    mon.names=mon.names, parm.names=parm.names, pos.theta=pos.theta,
    pos.mu=pos.mu, pos.nu=pos.nu, pos.sigma=pos.sigma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    theta <- parm[Data$pos.theta]
    mu <- matrix(parm[Data$pos.mu], Data$C, Data$J)
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    pi <- rep(0, Data$C)
    tab <- table(theta)
    pi[as.numeric(names(tab))] <- as.vector(tab)
    pi <- pi / sum(pi)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    theta.prior <- sum(dcat(theta, pi, log=TRUE))
    mu.prior <- sum(dnorm(mu, 0, matrix(nu, Data$C, Data$J), log=TRUE))
    nu.prior <- sum(dhalfcauchy(nu, 25, log=TRUE))
    pi.prior <- ddirichlet(pi, Data$alpha, log=TRUE)
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnorm(Data$Y, mu[theta,], sigma[theta], log=TRUE))
    ### Log-Posterior
    LP <- LL + theta.prior + mu.prior + nu.prior + pi.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,pi),
        yhat=rnorm(prod(dim(mu[theta,])), mu[theta,], sigma[theta]),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rcat(N, rep(1/C, C)), rep(0, C*J), rep(1, C), rep(1, C))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: cluster-level means
post.mu <- matrix(fit$Summary2[pos.mu, "Mean"], C, J)
cat("True mu:\n"); print(round(true.mu, 2))
cat("Posterior mean mu:\n"); print(round(post.mu, 2))

### Parameter recovery: mixing proportions (from Monitor)
post.pi <- colMeans(fit$Monitor[, 2:(C+1)])
cat("True pi:     ", round(true.pi, 3), "\n")
cat("Posterior pi:", round(post.pi, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
ground_truth[grep("pi", parm.names)]    <- true.pi[seq_along(grep("pi", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Cluster Analysis, Exploratory (ECA)

This is a nonparametric, model-based, cluster analysis, also called an
infinite mixture model or latent class cluster analysis, where the
number of clusters \\C\\ is unknown, and a Dirichlet process via
truncated stick-breaking is used to estimated the number of clusters.
The record-level cluster membership parameter vector, \\\theta\\, is a
vector of discrete parameters. Discrete parameters are not supported in
all algorithms. The example below is updated with the Slice sampler.
Bayesian exploratory cluster analysis via Dirichlet process mixtures was
developed by Ferguson [\[39\]](#ref39) and popularised for practical
inference by Escobar and West [\[38\]](#ref38). Dirichlet process
clustering is used in genomics for unsupervised discovery of gene
expression subtypes from microarray data.

#### Form

\\\textbf{Y}\_{i,j} \sim \mathcal{N}(\mu\_{\theta\[i\],j},
\sigma^2\_{\theta\[i\]}), \quad i=1,\dots,N, \quad j=1,\dots,J\\
\\\theta_i \sim \mathcal{CAT}(\pi\_{1:C}), \quad i=1,\dots,N\\
\\\mu\_{c,j} \sim \mathcal{N}(0, \nu^2_c), \quad c=1,\dots,C, \quad
j=1,\dots,J\\ \\\sigma_c \sim \mathcal{HC}(25), \quad c=1,\dots,C\\
\\\pi = \mathrm{Stick}(\delta)\\ \\\delta_c \sim \mathcal{BETA}(1,
\gamma), c=1,\dots,(C-1)\\ \\\gamma \sim \mathcal{G}(\alpha, \beta)\\
\\\alpha \sim \mathcal{HC}(25)\\ \\\beta \sim \mathcal{HC}(25)\\ \\\nu_c
\sim \mathcal{HC}(25), \quad c=1,\dots,C\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a Dirichlet process via truncated stick-breaking
with discrete cluster allocation parameters, hierarchical gamma
hyperpriors on the concentration parameter, and per-cluster indexing of
means and variances. This multi-level nonparametric structure falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42120)
N <- 100
J <- 4
C <- 5  # Maximum number of clusters to explore (true number is 3)

#### True parameter values (3 active clusters out of 5 maximum)
true.pi.active <- c(0.45, 0.35, 0.20)
true.mu <- matrix(c(-2, 0, 2, 1, -1, 0.5, 0.5, 1.5, -1, 0.8, -0.5, 1.2),
    3, J)
true.sigma <- c(0.4, 0.6, 0.5)

#### Generate cluster memberships and data
true.theta <- rcat(N, true.pi.active)
Y <- matrix(0, N, J)
for (i in 1:N)
    Y[i, ] <- rnorm(J, true.mu[true.theta[i], ], true.sigma[true.theta[i]])

#### Assemble Data list
mon.names <- c("LP", paste("pi[", 1:C, "]", sep=""))
parm.names <- as.parm.names(list(theta=rep(0,N), delta=rep(0,C-1),
    mu=matrix(0,C,J), nu=rep(0,C), sigma=rep(0,C), alpha=0, beta=0,
    gamma=0))
pos.theta <- grep("theta", parm.names)
pos.delta <- grep("delta", parm.names)
pos.mu <- grep("mu", parm.names)
pos.nu <- grep("nu", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.alpha <- grep("alpha", parm.names)
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    mu <- rnorm(Data$C*Data$J)
    nu <- runif(Data$C)
    sigma <- runif(Data$C)
    alpha <- runif(1)
    beta <- runif(1)
    gamma <- rgamma(1, alpha, beta)
    delta <- rev(sort(rbeta(Data$C-1, 1, gamma)))
    theta <- rcat(Data$N, Stick(delta))
    return(c(theta, delta, mu, nu, sigma, alpha, beta, gamma))
    }
Data <- list(C=C, J=J, N=N, PGF=PGF, Y=Y, mon.names=mon.names,
    parm.names=parm.names, pos.theta=pos.theta, pos.delta=pos.delta,
    pos.mu=pos.mu, pos.nu=pos.nu, pos.sigma=pos.sigma,
    pos.alpha=pos.alpha, pos.beta=pos.beta, pos.gamma=pos.gamma)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperhyperparameters
    alpha <- interval(parm[Data$pos.alpha], 1e-100, Inf)
    parm[Data$pos.alpha] <- alpha
    beta <- interval(parm[Data$pos.beta], 1e-100, Inf)
    parm[Data$pos.beta] <- beta
    ### Hyperparameters
    delta <- interval(parm[Data$pos.delta], 1e-10, 1-1e-10)
    parm[Data$pos.delta] <- delta
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    ### Parameters
    theta <- parm[Data$pos.theta]
    mu <- matrix(parm[Data$pos.mu], Data$C, Data$J)
    pi <- Stick(delta)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperhyperprior
    alpha.prior <- dhalfcauchy(alpha, 25, log=TRUE)
    beta.prior <- dhalfcauchy(beta, 25, log=TRUE)
    ### Log-Hyperprior
    delta.prior <- dStick(delta, gamma, log=TRUE)
    gamma.prior <- dgamma(gamma, alpha, beta, log=TRUE)
    nu.prior <- sum(dhalfcauchy(nu, 25, log=TRUE))
    ### Log-Prior
    theta.prior <- sum(dcat(theta, pi, log=TRUE))
    mu.prior <- sum(dnorm(mu, 0, matrix(nu, Data$C, Data$J), log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    LL <- sum(dnorm(Data$Y, mu[theta,], sigma[theta], log=TRUE))
    ### Log-Posterior
    LP <- LL + theta.prior + delta.prior + mu.prior + nu.prior +
        alpha.prior + beta.prior + gamma.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,pi),
        yhat=rnorm(prod(dim(mu[theta,])), mu[theta,], sigma[theta]),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rcat(N, rev(sort(rStick(C-1, 1)))), rep(0.5, C-1),
    rep(0, C*J), rep(1, C), rep(1, C), rep(1, 3))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: mixing proportions
post.pi <- colMeans(fit$Monitor[, 2:(C+1)])
cat("Posterior pi:", round(post.pi, 3), "\n")
cat("(True model has 3 active clusters with pi =",
    round(true.pi.active, 3), ")\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.mu]    <- true.mu[seq_along(pos.mu)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Dawid-Skene (Crowdsourcing) **\[NEW\]**

The Dawid-Skene model estimates the unknown true class labels of items
annotated by multiple noisy raters, jointly inferring each annotator’s
accuracy and the prevalence of each class from the observed labels
alone. In this simplified formulation each annotator \\j\\ has a single
accuracy parameter \\\pi_j\\ governing how often they assign the correct
label, with errors distributed uniformly across the remaining
categories, and the class prevalences follow a Dirichlet prior. The
model was introduced by Dawid and Skene [\[156\]](#ref156) as the
foundational approach to consensus estimation in the presence of
heterogeneous annotator quality. Crowdsourcing quality control in
platforms such as Amazon Mechanical Turk relies on Dawid-Skene to
aggregate noisy human labels for training machine learning classifiers.

#### Form

\\y\_{i,j} \sim \mathcal{CAT}(\textbf{p}\_{z_i,j}), \quad i=1,\dots,N,
\quad j=1,\dots,J\\ \\p\_{k,j,l} = \begin{cases} \pi_j & \text{if } l =
k \\ (1 - \pi_j)/(K-1) & \text{if } l \neq k \end{cases}\\ \\z_i \sim
\mathcal{CAT}(\textbf{prev}), \quad i=1,\dots,N\\ \\\pi_j \sim
\mathcal{BETA}(5, 1), \quad j=1,\dots,J\\ \\\textbf{prev} \sim
\mathcal{D}(\textbf{1}\_K)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it involves discrete latent class assignments \\z_i\\ that
must be marginalized or sampled as pseudo-continuous parameters,
combined with a confusion matrix structure linking annotator labels to
true classes. The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42208)
N <- 30   # Items
J <- 5    # Annotators
K <- 3    # Classes

#### True parameter values
true.acc <- c(0.90, 0.80, 0.85, 0.70, 0.75)
true.prev <- c(0.40, 0.35, 0.25)

#### Simulate true labels
true.z <- sample(1:K, N, replace=TRUE, prob=true.prev)

#### Simulate annotations
Y <- matrix(NA, N, J)
for (i in 1:N) {
    for (j in 1:J) {
        if (runif(1) < true.acc[j]) {
            Y[i,j] <- true.z[i]
        } else {
            Y[i,j] <- sample(setdiff(1:K, true.z[i]), 1)
        }
    }
}

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(
    z.raw=rep(0,N), acc=rep(0,J), prev.raw=rep(0,K-1)))
pos.z.raw <- grep("z.raw", parm.names)
pos.acc <- grep("acc", parm.names)
pos.prev.raw <- grep("prev.raw", parm.names)
PGF <- function(Data) {
    z.raw <- runif(Data$N, 0.5, Data$K + 0.5)
    acc <- rbeta(Data$J, 5, 1)
    prev.raw <- rnorm(Data$K - 1)
    return(c(z.raw, acc, prev.raw))
    }
Data <- list(N=N, J=J, K=K, Y=Y, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, pos.z.raw=pos.z.raw, pos.acc=pos.acc,
    pos.prev.raw=pos.prev.raw, true.z=true.z)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    z.raw <- interval(parm[Data$pos.z.raw], 0.5, Data$K + 0.5)
    parm[Data$pos.z.raw] <- z.raw
    z <- floor(z.raw)
    z[z < 1] <- 1L; z[z > Data$K] <- Data$K
    acc <- interval(parm[Data$pos.acc], 1e-100, 1 - 1e-100)
    parm[Data$pos.acc] <- acc
    prev.raw <- parm[Data$pos.prev.raw]
    log.prev <- c(prev.raw, 0) - log(sum(exp(c(prev.raw, 0))))
    prev <- exp(log.prev)
    ### Log-Prior
    acc.prior <- sum(dbeta(acc, 5, 1, log=TRUE))
    prev.prior <- ddirichlet(prev, rep(1, Data$K), log=TRUE)
    ### Log-Likelihood
    LL <- 0
    for (i in 1:Data$N) {
        LL <- LL + log.prev[z[i]]
        for (j in 1:Data$J) {
            if (Data$Y[i,j] == z[i]) {
                LL <- LL + log(acc[j])
            } else {
                LL <- LL + log((1 - acc[j]) / (Data$K - 1))
            }
        }
    }
    ### Log-Posterior
    LP <- LL + acc.prior + prev.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=sample(1:Data$K, Data$N*Data$J, replace=TRUE, prob=prev),
        parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(2, N), rep(0.8, J), rep(0, K-1))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery: annotator accuracies
cat("True accuracies: ", true.acc, "\n")
cat("Posterior mean:  ", round(fit$Summary2[pos.acc, "Mean"], 3), "\n")

#### Classification accuracy
post.z <- floor(fit$Summary2[pos.z.raw, "Mean"])
post.z[post.z < 1] <- 1; post.z[post.z > K] <- K
cat("Classification accuracy:",
    round(mean(post.z == true.z), 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.acc]  <- true.acc[seq_along(pos.acc)]
ground_truth[grep("prev", parm.names)] <- true.prev[seq_along(grep("prev", parm.names))]
ground_truth[grep("z", parm.names)]    <- true.z[seq_along(grep("z", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Mixture Model, Finite

This finite mixture model (FMM) imposes a multilevel structure on each
of the \\J\\ regression effects in \\\beta\\, so that mixture components
share a common residual standard deviation, \\\nu_m\\. Identifiability
is gained at the expense of some shrinkage. The record-level mixture
membership parameter vector, \\\theta\\, is a vector of discrete
parameters. Discrete parameters are not supported in all algorithms. The
example below is updated with the Slice sampler. Finite mixture models
were systematised by McLachlan and Peel [\[85\]](#ref85), with Bayesian
inference via data augmentation following Diebolt and Robert
[\[30\]](#ref30). Finite mixtures are used in astronomy to decompose the
distribution of galaxy recession velocities into distinct cluster
components.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_i =
\textbf{X}\_{i,1:J}\beta\_{\theta\[i\],1:J}, \quad i=1,\dots,N\\
\\\theta_i \sim \mathcal{CAT}(\pi\_{1:M}), \quad i=1,\dots,N\\
\\\beta\_{m,j} \sim \mathcal{N}(0, \nu^2_m), \quad j=1,\dots,J, \quad
m=2,\dots,M\\ \\\beta\_{1,j} \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\nu_m \sim \mathcal{HC}(25), \quad m=1,\dots,M\\
\\\sigma \sim \mathcal{HC}(25)\\ \\\pi\_{1:M} \sim
\mathcal{D}(\alpha\_{1:M})\\ \\\alpha_m = 1\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires discrete mixture allocation parameters
\\\theta_i\\ sampled from a Categorical distribution, with mixing
proportions derived from the empirical allocation frequencies and a
Dirichlet prior. The combination of discrete latent allocation,
per-mixture regression coefficients, and the Dirichlet-Categorical
hierarchy falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42176)
N <- 100
J <- 3   # Number of predictors (including intercept)
M <- 2   # Number of mixtures

#### True parameter values
true.pi <- c(0.6, 0.4)
true.beta <- matrix(c(2.0, 0.5, -0.3, 3.5, -0.8, 0.6), M, J)
true.sigma <- 0.5

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
for (j in 2:J) X[, j] <- CenterScale(X[, j])
true.theta <- rcat(N, true.pi)
mu.true <- rowSums(true.beta[true.theta, ] * X)
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
alpha <- rep(1, M)
mon.names <- "LP"
parm.names <- as.parm.names(list(theta=rep(0,N), beta=matrix(0,M,J),
    nu=rep(0,M), sigma=0))
pos.theta <- grep("theta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.nu <- grep("nu", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    theta <- rcat(Data$N, rep(1/Data$M, Data$M))
    nu <- runif(Data$M)
    beta <- rnormv(Data$M*Data$J, 0,
        cbind(1000, matrix(nu, Data$M, Data$J-1)))
    sigma <- runif(1)
    return(c(theta, beta, nu, sigma))
    }
Data <- list(J=J, M=M, N=N, PGF=PGF, X=X, alpha=alpha,
    mon.names=mon.names, parm.names=parm.names, pos.theta=pos.theta,
    pos.beta=pos.beta, pos.nu=pos.nu, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$M, Data$J)
    theta <- parm[Data$pos.theta]
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    pi <- rep(0, Data$M)
    tab <- table(theta)
    pi[as.numeric(names(tab))] <- as.vector(tab)
    pi <- pi / sum(pi)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0,
        cbind(1000, matrix(nu, Data$M, Data$J-1)), log=TRUE))
    theta.prior <- sum(dcat(theta, p=pi, log=TRUE))
    pi.prior <- ddirichlet(pi, Data$alpha, log=TRUE)
    nu.prior <- sum(dhalfcauchy(nu, 25, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- rowSums(beta[theta,] * Data$X)
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + theta.prior + pi.prior + nu.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rcat(N, rep(1/M, M)), rep(0, M*J), rep(1, M), 1)
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: regression coefficients
post.beta <- matrix(fit$Summary2[pos.beta, "Mean"], M, J)
cat("True beta:\n"); print(round(true.beta, 3))
cat("Posterior mean beta:\n"); print(round(post.beta, 3))

### Parameter recovery: residual standard deviation
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
ground_truth[grep("pi", parm.names)]    <- true.pi[seq_along(grep("pi", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Mixture Model, Infinite

This infinite mixture model (IMM) uses a Dirichlet process via truncated
stick-breaking. The record-level mixture membership parameter vector,
\\\theta\\, is a vector of discrete parameters. Discrete parameters are
not supported in all algorithms. The example below is updated with the
Slice sampler. Dirichlet process mixture models were introduced by
Ferguson [\[39\]](#ref39) and made computationally practical by Escobar
and West [\[38\]](#ref38). Infinite mixtures are applied in proteomics
for clustering mass spectrometry peaks without pre-specifying the number
of protein groups.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu_i =
\textbf{X}\_{i,1:J}\beta\_{\theta\[i\],1:J}, \quad i=1,\dots,N\\
\\\theta_i \sim \mathcal{CAT}(\pi\_{1:M}), \quad i=1,\dots,N\\
\\\beta\_{m,j} \sim \mathcal{N}(0, \nu^2_m), \quad j=1,\dots,J, \quad
m=2,\dots,M\\ \\\beta\_{1,j} \sim \mathcal{N}(0, 1000), \quad
j=1,\dots,J\\ \\\nu_m \sim \mathcal{HC}(25), \quad m=1,\dots,M\\
\\\sigma \sim \mathcal{HC}(25)\\ \\\pi = \mathrm{Stick}(\delta)\\
\\\delta_m \sim \mathcal{BETA}(1, \gamma), m=1,\dots,(M-1)\\ \\\gamma
\sim \mathcal{G}(\alpha, \iota)\\ \\\alpha \sim \mathcal{HC}(25)\\
\\\iota \sim \mathcal{HC}(25)\\

#### model_spec() notation

This model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires a Dirichlet process via truncated stick-breaking
with discrete mixture allocation parameters, hierarchical gamma
hyperpriors on the concentration parameter, and per-mixture regression
coefficients. This multi-level nonparametric mixture structure falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

``` r

set.seed(42177)
N <- 100
J <- 3   # Number of predictors (including intercept)
M <- 3   # Maximum number of mixtures to explore (true = 2)

#### True parameter values (2 active mixtures out of 3 maximum)
true.pi.active <- c(0.65, 0.35)
true.beta <- matrix(c(2.0, 0.5, -0.3, 3.5, -0.8, 0.6), 2, J)
true.sigma <- 0.4

#### Generate design matrix and response
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
for (j in 2:J) X[, j] <- CenterScale(X[, j])
true.theta <- rcat(N, true.pi.active)
mu.true <- rowSums(true.beta[true.theta, ] * X)
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- c("LP", as.parm.names(list(pi=rep(0,M))))
parm.names <- as.parm.names(list(theta=rep(0,N), delta=rep(0,M-1),
    beta=matrix(0,M,J), nu=rep(0,M), sigma=0, alpha=0, iota=0, gamma=0))
pos.theta <- grep("theta", parm.names)
pos.delta <- grep("delta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.nu <- grep("nu", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.alpha <- grep("alpha", parm.names)
pos.iota <- grep("iota", parm.names)
pos.gamma <- grep("gamma", parm.names)
PGF <- function(Data) {
    nu <- runif(Data$M)
    beta <- rnormv(Data$M*Data$J, 0,
        cbind(1000, matrix(nu, Data$M, Data$J-1)))
    sigma <- runif(1)
    alpha <- runif(1)
    iota <- runif(1)
    gamma <- rgamma(1, alpha, iota)
    delta <- rev(sort(rbeta(Data$M-1, 1, gamma)))
    theta <- rcat(Data$N, Stick(delta))
    return(c(theta, delta, beta, nu, sigma, alpha, iota, gamma))
    }
Data <- list(J=J, M=M, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.theta=pos.theta, pos.delta=pos.delta,
    pos.beta=pos.beta, pos.nu=pos.nu, pos.sigma=pos.sigma,
    pos.alpha=pos.alpha, pos.iota=pos.iota, pos.gamma=pos.gamma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperhyperparameters
    alpha <- interval(parm[Data$pos.alpha], 1e-100, Inf)
    parm[Data$pos.alpha] <- alpha
    iota <- interval(parm[Data$pos.iota], 1e-100, Inf)
    parm[Data$pos.iota] <- iota
    ### Hyperparameters
    delta <- interval(parm[Data$pos.delta], 1e-10, 1-1e-10)
    parm[Data$pos.delta] <- delta
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.gamma] <- gamma
    parm[Data$pos.nu] <- nu <- interval(parm[Data$pos.nu], 1e-100, Inf)
    ### Parameters
    beta <- matrix(parm[Data$pos.beta], Data$M, Data$J)
    theta <- parm[Data$pos.theta]
    pi <- Stick(delta)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperhyperprior
    alpha.prior <- dhalfcauchy(alpha, 25, log=TRUE)
    iota.prior <- dhalfcauchy(iota, 25, log=TRUE)
    ### Log-Hyperprior
    delta.prior <- dStick(delta, gamma, log=TRUE)
    gamma.prior <- dgamma(gamma, alpha, iota, log=TRUE)
    nu.prior <- sum(dhalfcauchy(nu, 25, log=TRUE))
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0,
        cbind(1000, matrix(nu, Data$M, Data$J-1)), log=TRUE))
    theta.prior <- sum(dcat(theta, pi, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- rowSums(beta[theta,]*Data$X)
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + delta.prior + theta.prior + nu.prior +
        sigma.prior + alpha.prior + iota.prior + gamma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,pi),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rcat(N, rev(sort(rStick(M-1, 1)))), rep(0.5, M-1),
    rep(0, M*J), rep(1, M), rep(1, 4))
```

#### Fitting and recovery

``` r

### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

### Parameter recovery: mixing proportions
post.pi <- colMeans(fit$Monitor[, 2:(M+1)])
cat("Posterior pi:", round(post.pi, 3), "\n")
cat("(True model has 2 active mixtures with pi =",
    round(true.pi.active, 3), ")\n")

### Parameter recovery: residual standard deviation
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.theta] <- true.theta[seq_along(pos.theta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Variable selection and regularization

### Ridge Regression

Ridge regression adds a shared normal prior \\\mathcal{N}(0, \lambda)\\
on all regression coefficients, shrinking them toward zero in proportion
to the precision \\1/\lambda\\. This is the Bayesian analogue of
L2-penalized least squares and stabilizes estimation when predictors are
collinear. Ridge regression was introduced by Hoerl and Kennard
[\[57\]](#ref57), with the Bayesian interpretation as a Gaussian prior
on coefficients clarified by Lindley and Smith [\[72\]](#ref72). Ridge
regression is used in chemometrics to handle multicollinearity when
predicting substance concentrations from highly correlated spectral
absorbance measurements.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2_1)\\ \\\mu =
\textbf{X}\beta\\ \\\beta_1 \sim \mathcal{N}(0, 1000)\\ \\\beta_j \sim
\mathcal{N}(0, \sigma^2_2), \quad j=2,\dots,J\\ \\\sigma_k \sim
\mathcal{HC}(25), \quad k=1,\dots,2\\

#### model_spec

The ridge structure can be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL by assigning distinct priors to the intercept and slope
coefficients, with the shrinkage variance governed by a half-Cauchy
hyperprior.

``` r

spec <- model_spec("
  y ~ Normal(mu, sigma1)
  mu = X %*% beta
  beta[1] ~ Normal(0, 1000)
  beta[j] ~ Normal(0, sigma2), j = 2,...,J
  sigma1 ~ HalfCauchy(25)
  sigma2 ~ HalfCauchy(25)
")
```

#### Ground truth and data

``` r

set.seed(42154)
N <- 200
J <- 8

#### True coefficients: intercept + 7 slopes, some deliberately small
true.beta  <- c(3.0, 1.2, -0.8, 0.05, 0.0, 0.4, -0.03, 0.1)
true.sigma <- 0.8

#### Correlated design matrix (collinearity)
Sigma.x <- outer(1:(J - 1), 1:(J - 1), function(i, j) 0.5^abs(i - j))
L <- chol(Sigma.x)
X.raw <- matrix(rnorm(N * (J - 1)), N, J - 1) %*% L
X <- cbind(1, X.raw)

#### Generate response
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- mu.true + rnorm(N, 0, true.sigma)

mon.names <- "LP"
parm.names <- as.parm.names(list(beta = rep(0, J), sigma = rep(0, 2)))
pos.beta  <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta  <- rnorm(Data$J)
    sigma <- runif(2)
    return(c(beta, sigma))
    }
Data <- list(J = J, PGF = PGF, X = X, mon.names = mon.names,
    parm.names = parm.names, pos.beta = pos.beta, pos.sigma = pos.sigma,
    y = y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnorm(beta, 0, c(1000, rep(sigma[2], Data$J-1)),
        log=TRUE))
    sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(1,J), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm selection
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended algorithm
eval(parse(text = rx$Recommend$primary$code))

#### Inspect
print(fit)
summary(fit)
caterpillar.plot(fit, Parms = c("beta", "sigma"))

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Variable Selection, BAL

The Bayesian Adaptive Lasso (BAL) performs variable selection by placing
Laplace (double-exponential) shrinkage priors on the regression
coefficients, with coefficient-specific scale parameters \\\gamma_j\\
drawn from an inverse-gamma distribution whose own shape and rate
hyperparameters \\\delta\\ and \\\tau\\ are learned from the data. This
three-level hierarchy allows the penalty to adapt to each coefficient
individually, shrinking irrelevant predictors toward zero while leaving
genuinely active effects relatively unconstrained. The intercept
receives a flat Laplace prior with scale 1000 so that it is effectively
unpenalised. The Bayesian adaptive LASSO was introduced by Leng et
al. [\[71\]](#ref71), extending the adaptive LASSO of Zou
[\[123\]](#ref123) to a fully Bayesian framework. BAL is applied in
genomics for selecting relevant SNPs from genome-wide association
studies while controlling for correlated markers.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_1 \sim \mathcal{L}(0, 1000)\\ \\\beta_j \sim \mathcal{L}(0,
\gamma_j), \quad j=2,\dots,J\\ \\\gamma_j \sim \mathcal{G}^{-1}(\delta,
\tau) \in \[0,\infty\], \quad j=2,\dots,J\\ \\\delta \sim
\mathcal{HC}(25)\\ \\\tau \sim \mathcal{HC}(25)\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The BAL model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the intercept \\\beta_1\\ receives a different prior
(Laplace with fixed scale 1000) from the remaining coefficients
\\\beta\_{2:J}\\ (Laplace with coefficient-specific scale \\\gamma_j\\).
This requires assigning heterogeneous priors to different elements of
the same parameter vector, which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

We simulate a regression problem with \\J=10\\ predictors where most
coefficients are near zero, so the adaptive shrinkage has genuine work
to do. The true intercept is 3.0, three predictors have moderate
effects, and the remaining six are exactly zero. Predictors are
standardised to match the centred-and-scaled convention used in the
original `demonsnacks`-based example.

``` r

set.seed(42148)
N <- 200
J <- 10

#### True parameter values
true.beta <- c(3.0, 1.5, -0.9, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0)
true.gamma <- c(0.8, 0.6, 0.3, 0.3, 0.3, 0.5, 0.3, 0.3, 0.3)
true.delta <- 1.2
true.tau <- 0.8
true.sigma <- 1.0

#### Generate design matrix (intercept + J-1 standardised predictors)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), gamma=rep(0,J-1), delta=0,
    tau=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.delta <- grep("delta", parm.names)
pos.tau <- grep("tau", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    delta <- runif(1)
    tau <- runif(1)
    sigma <- runif(1)
    gamma <- rinvgamma(Data$J-1, delta, tau)
    beta <- rlaplace(Data$J, 0, c(1,gamma))
    return(c(beta, gamma, delta, tau, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.delta=pos.delta, pos.tau=pos.tau, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperhyperparameters
    delta <- interval(parm[Data$pos.delta], 1e-100, Inf)
    parm[Data$pos.delta] <- delta
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    ### Hyperparameters
    gamma <- interval(parm[Data$pos.gamma], 1e-100, Inf)
    parm[Data$pos.gamma] <- gamma
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Hyperhyperprior
    delta.prior <- dhalfcauchy(delta, 25, log=TRUE)
    tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
    ### Log-Hyperprior
    gamma.prior <- sum(dinvgamma(gamma, delta, tau, log=TRUE))
    ### Log-Prior
    beta.prior <- sum(dlaplace(beta, 0, c(1000, gamma), log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + gamma.prior + delta.prior + tau.prior +
        sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(0,J-1), rep(1,3))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: regression coefficients
cat("True beta:      ", true.beta, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Parameter recovery: shrinkage scales gamma
cat("True gamma:      ", true.gamma, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.gamma, "Mean"], 3), "\n")

#### Parameter recovery: hyperparameters and sigma
cat("delta -- true:", true.delta,
    " post. mean:", round(fit$Summary2[pos.delta, "Mean"], 3), "\n")
cat("tau   -- true:", true.tau,
    " post. mean:", round(fit$Summary2[pos.tau, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.delta] <- true.delta[seq_along(pos.delta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.tau]   <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Variable Selection, Horseshoe

The horseshoe prior is a global-local shrinkage prior that provides
stronger shrinkage of noise coefficients toward zero than the ordinary
LASSO, while leaving large signals essentially unpenalised. Each
regression coefficient \\\beta_j\\ (for \\j \ge 2\\) is drawn from a
horseshoe distribution parameterised by a local scale \\\lambda_j\\ and
a global scale \\\tau\\, both given half-Cauchy(1) priors. The intercept
receives a vague normal prior, and the residual standard deviation a
half-Cauchy(25) prior. The horseshoe prior was introduced by Carvalho et
al. [\[21\]](#ref21) as a default prior for sparse signal recovery with
near-minimax properties. Horseshoe regression is used in
high-dimensional neuroscience studies to identify the few brain voxels
associated with a cognitive task from thousands of candidates.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_1 \sim \mathcal{N}(0, 1000)\\ \\\beta_j \sim
\mathcal{HS}(\lambda_j, \tau), \quad j=2,\dots,J\\ \\\lambda_j \sim
\mathcal{HC}(1), \quad j=2,\dots,J\\ \\\tau \sim \mathcal{HC}(1)\\
\\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The horseshoe model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the intercept \\\beta_1\\ uses a normal prior while the
remaining coefficients \\\beta\_{2:J}\\ use the horseshoe distribution
with per-coefficient local scales \\\lambda_j\\. Assigning heterogeneous
priors to different elements of the same parameter vector, each with its
own auxiliary parameter, falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

We generate a sparse regression with \\J=10\\ predictors. Three
coefficients are meaningfully different from zero and the rest are
exactly zero, providing a clear test of the horseshoe’s ability to
separate signal from noise.

``` r

set.seed(42170)
N <- 200
J <- 10

#### True parameter values
true.beta <- c(2.0, 1.8, 0.0, -1.2, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0)
true.lambda <- rep(0.5, J - 1)
true.tau <- 0.3
true.sigma <- 0.8

#### Generate design matrix (intercept + J-1 predictors)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), lambda=rep(0,J-1),
    tau=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.tau <- grep("tau", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    lambda <- runif(Data$J-1)
    tau <- runif(1)
    sigma <- runif(1)
    return(c(beta, lambda, tau, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.lambda=pos.lambda,
    pos.tau=pos.tau, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta[1], 0, 1000, log=TRUE),
        dhs(beta[-1], lambda, tau, log=TRUE))
    lambda.prior <- sum(dhalfcauchy(lambda, 1, log=TRUE))
    tau.prior <- dhalfcauchy(tau, 1, log=TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(1,J-1), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: regression coefficients
cat("True beta:      ", true.beta, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Parameter recovery: local scales lambda
cat("True lambda:     ", true.lambda, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")

#### Parameter recovery: global scale and sigma
cat("tau   -- true:", true.tau,
    " post. mean:", round(fit$Summary2[pos.tau, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.tau]    <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Variable Selection, LASSO

The Bayesian LASSO places a penalised prior on the regression
coefficients through the `dlasso` distribution, which combines a normal
likelihood conditional on coefficient-specific precision parameters
\\\tau_j\\ that are themselves exponentially distributed with rate
governed by a global penalty \\\lambda\\. This formulation is the Park
and Casella [\[92\]](#ref92) representation that integrates the LASSO
penalty into a proper Bayesian hierarchy. The intercept receives a vague
normal prior and is not penalised. The Bayesian LASSO was developed by
Park and Casella [\[92\]](#ref92) using Laplace (double-exponential)
priors on regression coefficients. LASSO is applied in text
classification to select informative features from thousands of word
frequency predictors.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_1 \sim \mathcal{N}(0, 1000)\\ \\\beta_j \sim \mathcal{LASSO}(0,
\sigma, \tau_j, \lambda), \quad j=2,\dots,J\\

#### model_spec() notation

The LASSO model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because the intercept \\\beta_1\\ uses a vague normal prior while
the remaining coefficients \\\beta\_{2:J}\\ use the `dlasso`
distribution with per-coefficient local precisions \\\tau_j\\. Assigning
heterogeneous priors to different elements of the same parameter vector
falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

We simulate a sparse regression with \\J=10\\ predictors. Three
coefficients are active while the remaining six are exactly zero,
providing a clear target for the LASSO penalty.

``` r

set.seed(42171)
N <- 200
J <- 10

#### True parameter values
true.beta <- c(1.5, 2.0, 0.0, 0.0, -1.1, 0.0, 0.0, 0.0, 0.0, 0.8)
true.sigma <- 0.9
true.tau <- rep(0.5, J - 1)
true.lambda <- 1.0

#### Generate design matrix (intercept + J-1 predictors)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, tau=rep(0,J-1),
    lambda=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.tau <- grep("tau", parm.names)
pos.lambda <- grep("lambda", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    tau <- runif(Data$J-1)
    lambda <- runif(1)
    return(c(beta, sigma, tau, lambda))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma,
    pos.tau=pos.tau, pos.lambda=pos.lambda, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
    lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    ### Log-Prior
    beta.prior <- sum(dnormv(beta[1], 0, 1000, log=TRUE),
        dlasso(beta[-1], sigma, tau, lambda, log=TRUE))
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1, rep(1,J-1), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using the recommended algorithm (shown in print output above)
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: regression coefficients
cat("True beta:      ", true.beta, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Parameter recovery: sigma and lambda
cat("sigma  -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")
cat("lambda -- true:", true.lambda,
    " post. mean:", round(fit$Summary2[pos.lambda, "Mean"], 3), "\n")

#### Parameter recovery: local precisions tau
cat("True tau:       ", true.tau, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.tau, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]   <- true.beta[seq_along(pos.beta)]
ground_truth[pos.lambda] <- true.lambda[seq_along(pos.lambda)]
ground_truth[pos.sigma]  <- true.sigma[seq_along(pos.sigma)]
ground_truth[pos.tau]    <- true.tau[seq_along(pos.tau)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Variable Selection, RJ

This example uses the RJ (Reversible-Jump) algorithm of the `lucifer`
function for variable selection and Bayesian Model Averaging (BMA). In a
high-dimensional sparse regression with \\J=100\\ predictors, 90% of
which have zero true coefficients, the reversible-jump mechanism
proposes birth and death moves that add or remove predictors from the
active set. The binomial prior on model size (with mean \\0.4 \times 99
\approx 40\\) regularises complexity, and the algorithm explores the
joint space of model indicators and continuous parameters
simultaneously. Other MCMC algorithms will not perform variable
selection with this example as presented; the RJ algorithm must be used.
Reversible jump MCMC was introduced by Green [\[53\]](#ref53) for
exploring models of varying dimension within a single Markov chain.
Reversible jump methods are used in astrophysics for determining the
number of spectral emission lines in gamma-ray burst data.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

The RJ variable selection model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires reversible-jump specifications (birth/death move
probabilities, model-size prior, and a selectable indicator vector) that
operate on the discrete model-indicator space. These meta-algorithmic
constructs are external to the model likelihood and prior, and fall
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

We simulate a large sparse regression where only 10% of the 100
predictors have nonzero coefficients. The seed is set so that the
sparsity pattern is reproducible and recovery of the active set can be
verified.

``` r

set.seed(42172)
N <- 1000
J <- 100
true.sigma <- 0.5

#### True coefficients: 10% nonzero, rest exactly zero
true.beta <- rep(0, J)
true.beta[1] <- 1.5  # intercept
active <- sort(sample(2:J, round(J * 0.1)))
true.beta[active] <- runif(length(active), -3, 3)

#### Generate design matrix (intercept column + J-1 predictors)
X <- matrix(1, N, J)
for (j in 2:J) X[,j] <- rnorm(N, runif(1, -3, 3), runif(1, 0.1, 1))

#### Generate response from the true model
mu.true <- X %*% true.beta
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

#### Reversible-Jump specifications
bin.n <- J-1 #Maximum allowable model size
bin.p <- 0.4 #Most probable size: bin.p x bin.n is binomial mean and median
parm.p <- rep(1/J,J+1)
selectable <- c(0, rep(1,J-1), 0)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- GIV(Model, Data, PGF=TRUE)
```

#### Fitting and recovery

``` r

#### fit using the RJ algorithm (required for reversible-jump variable selection)
fit <- lucifer(Model, Data, Initial.Values,
    Covar = NULL, Iterations = 100000, Status = 1000, Thinning = 100,
    Algorithm = "RJ",
    Specs = list(bin.n = bin.n, bin.p = bin.p,
                 parm.p = parm.p, selectable = selectable,
                 selected = rep(0, length(Initial.Values))))
print(fit)
summary(fit)

#### Parameter recovery: identify active predictors
post.beta <- fit$Summary2[pos.beta, "Mean"]
cat("True active set:      ", sort(active), "\n")
cat("Recovered active set: ",
    sort(which(abs(post.beta[-1]) > 0.1) + 1), "\n")

#### Parameter recovery: nonzero coefficients
for (j in active)
    cat("beta[", j, "] -- true:", round(true.beta[j], 3),
        " post. mean:", round(post.beta[j], 3), "\n")

#### Parameter recovery: sigma
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Variable Selection, SSVS

This example uses a modified form of the random-effects (or global
adaptation) Stochastic Search Variable Selection (SSVS) algorithm
presented in O’Hara and Sillanpää [\[90\]](#ref90), which selects
variables according to practical significance rather than statistical
significance. Here, SSVS is applied to linear regression, though this
method is widely applicable. For \\J\\ variables, each regression effect
\\\beta_j\\ is conditional on \\\gamma_j\\, a binary inclusion variable.
Each \\\beta_j\\ is a discrete mixture distribution with respect to
\\\gamma_j = 0\\ or \\\gamma_j = 1\\, with precision 100 or
\\\beta\_\sigma = 0.1\\, respectively. As with other representations of
SSVS, these precisions may require tuning. Stochastic search variable
selection was proposed by George and McCulloch [\[49\]](#ref49) using
spike-and-slab priors to identify active predictors. SSVS is applied in
econometrics for selecting macroeconomic predictors of GDP growth from a
large candidate set.

The binary inclusion variables are discrete parameters, and discrete
parameters are not supported in all algorithms. The example below is
updated with the Slice sampler.

When the goal is to select the best model, each \\\textbf{X}\_{1:N,j}\\
is retained for a future run when the posterior mean of \\\gamma_j \ge
0.5\\. When the goal is model-averaging, the results of this model may
be used directly, which would please L. J. Savage, who said that “models
should be as big as an elephant” [\[32\]](#ref32).

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}
\beta\\ \\\beta_1 \sim \mathcal{N}(0, 1000)\\ \\(\beta_j \| \gamma_j)
\sim (1 - \gamma_j)\mathcal{N}(0, 0.01) + \gamma_j \mathcal{N}(0,
\beta^2\_\sigma) \quad j=2,\dots,J\\ \\\beta\_\sigma \sim
\mathcal{HC}(25)\\ \\\gamma_j \sim \mathcal{BERN}(1/(J-1)), \quad
j=1,\dots,(J-1)\\ \\\sigma \sim \mathcal{HC}(25)\\

#### model_spec() notation

The SSVS model cannot be expressed in the
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
DSL because it requires discrete binary inclusion indicators
\\\gamma_j\\ that switch between two mixture components. The discrete
parameters, the conditional mixture prior structure, and the imperative
logic that sets `beta.sigma[gamma == 0] <- 0.1` all fall outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).
The manual `Model` function below handles this directly.

#### Ground truth and data

We simulate a regression with \\J=10\\ predictors where four
coefficients are genuinely active (\\\gamma_j = 1\\) and five are
excluded (\\\gamma_j = 0\\), so the SSVS mechanism has a clear signal to
recover. The intercept is always included and not subject to selection.

``` r

set.seed(42173)
N <- 200
J <- 10

#### True parameter values
true.beta <- c(2.5, 1.3, 0.0, -0.9, 0.0, 0.0, 0.7, 0.0, 0.0, -1.1)
true.gamma <- c(  1,   0,   1,   0,   0,   1,   0,   0,   1)  # J-1 indicators
true.b.sd <- 1.5    # beta.sigma when gamma=1
true.sigma <- 0.8

#### Generate design matrix (intercept + J-1 standardised predictors)
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))

#### Generate response from the true model (gamma masks inactive coefficients)
mu.true <- X %*% (true.beta * c(1, true.gamma))
y <- rnorm(N, mu.true, true.sigma)

#### Assemble Data list
mon.names <- c("LP", "min.beta.sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), gamma=rep(0,J-1),
    b.sd=0, sigma=0))
pos.beta <- grep("beta", parm.names)
pos.gamma <- grep("gamma", parm.names)
pos.b.sd <- grep("b.sd", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    gamma <- rep(1,Data$J-1)
    b.sd <- rnorm(1)
    sigma <- runif(1)
    return(c(beta, gamma, b.sd, sigma))
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.gamma=pos.gamma,
    pos.b.sd=pos.b.sd, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Hyperparameters
    beta.sigma <- interval(parm[Data$pos.b.sd], 1e-100, Inf)
    parm[Data$pos.b.sd] <- beta.sigma
    ### Parameters
    beta <- parm[Data$pos.beta]
    gamma <- parm[Data$pos.gamma]
    beta.sigma <- rep(beta.sigma, Data$J-1)
    beta.sigma[gamma == 0] <- 0.1
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    ### Log-Hyperprior
    beta.sigma.prior <- sum(dhalfcauchy(beta.sigma, 25, log=TRUE))
    ### Log-Prior
    beta.prior <- sum(dnorm(beta, 0, c(sqrt(1000), beta.sigma, log=TRUE)))
    gamma.prior <- sum(dbern(gamma, 1/(Data$J-1), log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta*c(1,gamma)))
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    ### Log-Posterior
    LP <- LL + beta.prior + beta.sigma.prior + gamma.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, min(beta.sigma)),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), rep(1,J-1), rep(1,2))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation (Slice is recommended for discrete parameters)
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit using Slice (required for discrete gamma indicators)
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
summary(fit)

#### Parameter recovery: inclusion indicators
post.gamma <- round(fit$Summary2[pos.gamma, "Mean"], 2)
cat("True gamma:     ", true.gamma, "\n")
cat("Posterior mean: ", post.gamma, "\n")
cat("Selected (>0.5):", which(post.gamma >= 0.5), "\n")
cat("True active:    ", which(true.gamma == 1), "\n")

#### Parameter recovery: regression coefficients
cat("True beta:      ", true.beta, "\n")
cat("Posterior mean: ", round(fit$Summary2[pos.beta, "Mean"], 3), "\n")

#### Parameter recovery: beta.sigma and sigma
cat("b.sd  -- true:", true.b.sd,
    " post. mean:", round(fit$Summary2[pos.b.sd, "Mean"], 3), "\n")
cat("sigma -- true:", true.sigma,
    " post. mean:", round(fit$Summary2[pos.sigma, "Mean"], 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.gamma] <- true.gamma[seq_along(pos.gamma)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Phylogenetic and evolutionary models

### GTR+G substitution model (fixed tree) **\[NEW\]**

This example estimates the parameters of the general time-reversible
(GTR) nucleotide substitution model with discrete gamma rate
heterogeneity, conditioned on a fixed phylogenetic tree. The GTR model
[\[164\]](#ref164) parameterizes the instantaneous rate matrix \\Q\\
through six symmetric exchangeability parameters and four equilibrium
base frequencies; gamma rate heterogeneity across sites
[\[165\]](#ref165) adds one shape parameter. The likelihood is computed
via Felsenstein’s pruning algorithm [\[166\]](#ref166), using the
`phangorn` package [\[167\]](#ref167). Data are the 15 European
woodmouse cytochrome *b* sequences from the `ape` package
[\[168\]](#ref168), with a neighbor-joining tree fixed for inference.
This model is central to molecular phylogenetics and serves as the
foundation for more complex phylodynamic analyses involving molecular
clocks and coalescent demographic models.

#### Form

\\L(D \mid T, \theta) = \prod\_{k=1}^{K} L_k^{w_k}, \quad L_k =
\sum\_{r} \pi_r \prod\_{\text{branches}} P\_{ij}(b \cdot c)\\ \\P(t) =
\exp(Qt), \quad Q\_{ij} = \rho\_{ij} \pi_j \text{ for } i \neq j\\
\\\boldsymbol{\pi} \sim \text{Dirichlet}(1,1,1,1)\\ \\\rho_m \sim
\text{Exp}(1), \quad m = 1,\ldots,5, \quad \rho_6 = 1\\ \\\alpha \sim
\text{Exp}(1)\\

where \\K\\ is the number of unique site patterns, \\w_k\\ are pattern
weights, \\c\\ is the rate category from a discrete gamma distribution
with shape \\\alpha\\ and 4 categories, and \\\rho_6\\ (GT rate) is
fixed at 1 for identifiability.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the phylogenetic likelihood requires external tree and alignment
objects passed through
[`phangorn::pml()`](https://klausvigo.github.io/phangorn/reference/pml.html),
which falls outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

library(ape)
library(phangorn)

#### Load woodmouse sequences
data(woodmouse, package = "ape")
N <- ncol(woodmouse)

#### Build and root NJ tree
dm <- ape::dist.dna(woodmouse, model = "F84")
tree <- ape::root(ape::nj(dm), outgroup = "No305", resolve.root = TRUE)
pdat <- phangorn::as.phyDat(woodmouse)

#### ML reference values (from optim.pml)
fit_ml <- phangorn::optim.pml(
    phangorn::pml(tree, pdat), model = "GTR",
    optGamma = TRUE, rearrangement = "none",
    control = phangorn::pml.control(trace = 0))
cat("ML logLik:", fit_ml$logLik, "\n")
cat("ML bf:", round(fit_ml$bf, 4), "\n")
cat("ML Q:", round(fit_ml$Q, 4), "\n")
cat("ML shape:", round(fit_ml$shape, 4), "\n")

#### Helper: phylogenetic log-likelihood via Felsenstein pruning
phylo_loglik <- function(bf, Q, shape, tree, pdat, k = 4) {
    fit <- phangorn::pml(tree, pdat, bf = bf, Q = Q, k = k, shape = shape)
    site_ll <- rep(fit$siteLik, attr(pdat, "weight"))
    list(logLik = fit$logLik, site_ll = site_ll)
}

#### Assemble Data list
mon.names <- c("LP", "LL")
parm.names <- c("alpha_A", "alpha_C", "alpha_G",
                "q_AC", "q_AG", "q_AT", "q_CG", "q_CT",
                "shape")
PGF <- function(Data) {
    alphas <- rnorm(3, 0, 1)  # softmax log-ratios
    rates <- rexp(5, 1)
    shape <- rexp(1, 1)
    c(alphas, rates, shape)
}
Data <- list(N=N, PGF=PGF, mon.names=mon.names, parm.names=parm.names,
    tree=tree, pdat=pdat)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Base frequencies via softmax (unconstrained parameterization)
    alpha <- c(parm[1:3], 0)
    alpha <- alpha - max(alpha)
    bf    <- exp(alpha) / sum(exp(alpha))
    ### GTR rates
    q_AC  <- interval(parm[4], 1e-4, 100)
    q_AG  <- interval(parm[5], 1e-4, 100)
    q_AT  <- interval(parm[6], 1e-4, 100)
    q_CG  <- interval(parm[7], 1e-4, 100)
    q_CT  <- interval(parm[8], 1e-4, 100)
    parm[4:8] <- c(q_AC, q_AG, q_AT, q_CG, q_CT)
    Q <- c(q_AC, q_AG, q_AT, q_CG, q_CT, 1)
    shape <- interval(parm[9], 0.01, 100)
    parm[9] <- shape
    ### Log-Prior
    pi.prior <- sum(log(bf))  # Dirichlet(1,1,1,1) + Jacobian
    q.prior  <- sum(dexp(c(q_AC, q_AG, q_AT, q_CG, q_CT), 1, log = TRUE))
    s.prior  <- dexp(shape, 1, log = TRUE)
    ### Log-Likelihood
    res <- phylo_loglik(bf, Q, shape, Data$tree, Data$pdat)
    LL  <- res$logLik
    ### Log-Posterior
    LP <- LL + pi.prior + q.prior + s.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, LL),
        yhat=rep(LL / Data$N, Data$N), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(0, 0, 0, 1, 1, 1, 1, 1, 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit with AMWG (phangorn likelihood has no analytical gradients)
fit <- lucifer(Model, Data, Initial.Values,
    Iterations = 20000, Status = 5000, Thinning = 10,
    Algorithm = "AMWG",
    Specs = list(B = NULL, n = 0, Periodicity = 50))

#### Diagnostics
print(fit)
Consort(fit)

#### Parameter recovery: base frequencies (transform softmax alphas)
alpha_post <- c(fit$Summary2[1:3, "Mean"], 0)
alpha_post <- alpha_post - max(alpha_post)
bf_post <- exp(alpha_post) / sum(exp(alpha_post))
cat("ML bf:", round(fit_ml$bf, 4), "\n")
cat("Posterior mean bf:", round(bf_post, 4), "\n")

#### Parameter recovery: GTR rates
cat("ML Q:", round(fit_ml$Q, 4), "\n")
cat("Posterior mean Q:",
    round(c(fit$Summary2[4:8, "Mean"], 1), 4), "\n")

#### Parameter recovery: gamma shape
cat("ML shape:", round(fit_ml$shape, 4), "\n")
cat("Posterior mean shape:", round(fit$Summary2[9, "Mean"], 4), "\n")

#### Posterior predictive checks
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
```

------------------------------------------------------------------------

### Bayesian skyline coalescent **\[NEW\]**

The Bayesian skyline model [\[169\]](#ref169) reconstructs the
demographic history of a population from a genealogy by estimating
piecewise-constant effective population sizes across grouped coalescent
intervals. Under Kingman’s coalescent [\[170\]](#ref170), the waiting
time for \\k\\ lineages to coalesce is exponentially distributed with
rate \\\binom{k}{2}/N_e\\, yielding an analytical log-likelihood that
requires no external phylogenetic library. A random walk prior on the
log-population sizes imposes temporal smoothness. This model is a core
tool in phylodynamic inference, used by BEAST [\[171\]](#ref171) to
infer epidemic growth curves, population bottlenecks, and expansion
events from viral genealogies. This self-contained example uses
simulated coalescent data with a known demographic trajectory for
validation.

#### Form

\\\Delta t_j \mid N_j \sim
\text{Exp}\left(\frac{\binom{k_j}{2}}{N_j}\right), \quad
j=1,\ldots,n-1\\ \\\log N_e^{(1)} \sim \mathcal{N}(0, 25)\\ \\\log
N_e^{(m)} \mid \log N_e^{(m-1)} \sim \mathcal{N}(\log N_e^{(m-1)},
\tau^2), \quad m=2,\ldots,M\\ \\\tau \sim \mathcal{HC}(1)\\

where \\k_j\\ is the number of lineages in interval \\j\\, \\\Delta
t_j\\ is the interval duration, and \\M\\ groups of consecutive
intervals share a common \\N_e\\.

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because the coalescent likelihood is a custom function over
inter-coalescence waiting times and lineage counts, which cannot be
expressed in the standard observation-level likelihood format of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42)

#### True demographic trajectory (piecewise constant)
M <- 5
true.Ne <- c(10, 5, 2, 8, 15)  # population sizes across time
true.tau <- 0.5                  # smoothness

#### Simulate coalescent intervals
n_tips <- 30
n_int <- n_tips - 1
ints_per_group <- ceiling(n_int / M)
groups <- rep(seq_len(M), each = ints_per_group)[seq_len(n_int)]
lineages <- n_tips:2

intervals <- numeric(n_int)
for (j in seq_len(n_int)) {
    choose_k <- lineages[j] * (lineages[j] - 1) / 2
    intervals[j] <- rexp(1, rate = choose_k / true.Ne[groups[j]])
}

#### Coalescent log-likelihood function
coal_loglik <- function(Ne, intervals, lineages, groups) {
    choose_k <- lineages * (lineages - 1) / 2
    ll <- 0
    for (j in seq_along(intervals)) {
        Ne_j <- Ne[groups[j]]
        ll <- ll + log(choose_k[j]) - log(Ne_j) -
              choose_k[j] * intervals[j] / Ne_j
    }
    ll
}

#### Group midpoint times for plotting
cum_t <- c(0, cumsum(intervals))
group_times <- tapply(
    (cum_t[seq_len(n_int)] + cum_t[seq_len(n_int) + 1]) / 2,
    groups, mean)

#### Assemble Data list
mon.names <- c("LP")
parm.names <- c(paste0("log_Ne_", seq_len(M)), "tau")
PGF <- function(Data) c(rnorm(Data$M, 0, 2), runif(1, 0.1, 2))
Data <- list(M=M, N=n_int, PGF=PGF, mon.names=mon.names,
    parm.names=parm.names, intervals=intervals, lineages=lineages,
    groups=groups, group_times=as.numeric(group_times),
    true.Ne=true.Ne)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    log_Ne <- parm[seq_len(Data$M)]
    Ne <- exp(log_Ne)
    tau <- interval(parm[Data$M + 1], 0.01, 10)
    parm[Data$M + 1] <- tau
    ### Log-Prior: random walk on log-Ne
    first.prior <- dnorm(log_Ne[1], 0, 5, log = TRUE)
    rw.prior <- 0
    for (m in 2:Data$M) {
        rw.prior <- rw.prior +
            dnorm(log_Ne[m], log_Ne[m - 1], tau, log = TRUE)
    }
    tau.prior <- dhalfcauchy(tau, 1, log = TRUE)
    ### Log-Likelihood
    LL <- coal_loglik(Ne, Data$intervals, Data$lineages, Data$groups)
    ### Log-Posterior
    LP <- LL + first.prior + rw.prior + tau.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
        yhat=rep(LL / Data$N, Data$N), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0, M), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### Fit
fit <- lucifer(Model, Data, Initial.Values,
    Iterations = 30000, Status = 5000, Thinning = 10,
    Algorithm = "NUTS",
    Specs = list(A = 5000, delta = 0.8, epsilon = NULL, Lmax = 10))

#### Diagnostics
print(fit)
Consort(fit)

#### Parameter recovery: population sizes
post.Ne <- exp(fit$Posterior2[, seq_len(M)])
cat("True Ne:", true.Ne, "\n")
cat("Posterior median Ne:", round(apply(post.Ne, 2, median), 2), "\n")

#### Skyline plot
library(ggplot2)
sky_df <- data.frame(
    time = Data$group_times,
    Ne_median = apply(post.Ne, 2, median),
    Ne_lower  = apply(post.Ne, 2, quantile, 0.025),
    Ne_upper  = apply(post.Ne, 2, quantile, 0.975),
    Ne_true   = true.Ne
)
ggplot2::ggplot(sky_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(
        ggplot2::aes(ymin = Ne_lower, ymax = Ne_upper),
        fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = Ne_median),
        color = "steelblue", linewidth = 1.2) +
    ggplot2::geom_point(ggplot2::aes(y = Ne_true),
        color = "red", size = 3, shape = 17) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(x = "Time (coalescent units)",
        y = expression(N[e]),
        title = "Bayesian skyline reconstruction") +
    ggplot2::theme_minimal()

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[grep("Ne", parm.names)]  <- true.Ne[seq_along(grep("Ne", parm.names))]
ground_truth[grep("tau", parm.names)] <- true.tau[seq_along(grep("tau", parm.names))]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

  

## Model comparison and computational methods

### Approximate Bayesian Computation (ABC)

Approximate Bayesian Computation (ABC), also called likelihood-free
estimation, replaces the intractable likelihood with a distance-based
approximation between the observed data and data simulated from the
model. This example applies ABC to linear regression, where the
log-likelihood is replaced with the negative sum of absolute differences
between \\\textbf{y}\\ and \\\textbf{y}^{rep}\\. Although linear
regression has an analytic likelihood, it serves as a transparent
demonstration of the ABC mechanism within `lucifer`. ABC was pioneered
by Tavaré et al. [\[103\]](#ref103) and Pritchard et
al. [\[95\]](#ref95), with the regression adjustment of Beaumont et
al. [\[6\]](#ref6) improving efficiency. ABC is essential in population
genetics where the likelihood of coalescent models is intractable,
enabling inference on effective population sizes and migration rates.

#### Form

\\\textbf{y} = \mu + \epsilon\\ \\\mu = \textbf{X} \beta\\ \\\beta_j
\sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because it uses a custom distance-based likelihood approximation rather
than a standard distributional likelihood, which falls outside the
declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42205)
N <- 100
J <- 4

#### True parameter values
true.beta <- c(2.0, 1.5, -0.8, 0.4)

#### Design matrix with intercept
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#### Generate response
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- mu.true + rnorm(N, 0, 1.0)

#### Assemble Data list
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J)))
pos.beta <- grep("beta", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    return(beta)
    }
Data <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    ### Log-Likelihood Approximation
    mu <- as.vector(tcrossprod(Data$X, t(beta)))
    epsilon <- Data$y - mu
    sigma <- sd(epsilon)
    LL <- -sum(abs(epsilon))
    ### Log-Posterior Approximation
    LP <- LL + beta.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J))
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Parameter recovery
post.beta <- fit$Summary2[pos.beta, "Mean"]
cat("True beta:      ", true.beta, "\n")
cat("Posterior mean:  ", round(post.beta, 3), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta] <- true.beta[seq_along(pos.beta)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

### Conditional Predictive Ordinate

The conditional predictive ordinate (CPO) measures the predictive
adequacy of a model for each observation. By monitoring the record-level
inverse likelihood \\\mathrm{InvL}\_i\\ during MCMC, one obtains
\\\mathrm{CPO}\_i\\ as the inverse of the posterior mean of
\\\mathrm{InvL}\_i\\. Observations with inverse CPO (ICPO) above 40 are
potential outliers; above 70 are extreme. This example adds CPO
monitoring to a standard linear regression. The conditional predictive
ordinate was proposed by Geisser [\[45\]](#ref45) and developed by
Gelfand et al. [\[46\]](#ref46) as a leave-one-out cross-validation
diagnostic. CPO values identify influential observations in survival
analysis, flagging patients whose outcomes are poorly predicted by the
fitted model.

#### Form

\\\textbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\ \\\mu = \textbf{X}\beta\\
\\\beta_j \sim \mathcal{N}(0, 1000), \quad j=1,\dots,J\\ \\\sigma \sim
\mathcal{HC}(25)\\

#### model_spec() notation

This model is not compatible with
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
because CPO requires monitoring the record-level inverse likelihood
within the model function, a custom diagnostic calculation that falls
outside the declarative scope of
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md).

#### Ground truth and data

``` r

set.seed(42206)
N <- 50
J <- 4

#### True parameter values
true.beta  <- c(3.0, 1.0, -0.5, 0.8)
true.sigma <- 1.0

#### Design matrix with intercept
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
for (j in 2:J) X[,j] <- CenterScale(X[,j])

#### Generate response with two outliers
mu.true <- as.vector(tcrossprod(X, t(true.beta)))
y <- rnorm(N, mu.true, true.sigma)
y[6]  <- y[6] + 8   # possible outlier
y[8]  <- y[8] + 15  # extreme value

#### Assemble Data list
mon.names <- c("LP", as.parm.names(list(InvL=rep(0,N))))
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
    beta <- rnorm(Data$J)
    sigma <- runif(1)
    return(c(beta, sigma))
    }
Data <- list(J=J, N=N, PGF=PGF, X=X, mon.names=mon.names,
    parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
```

#### Model

``` r

Model <- function(parm, Data)
    {
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    ### Log-Prior
    beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
    sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
    ### Log-Likelihood
    mu <- tcrossprod(Data$X, t(beta))
    LL <- dnorm(Data$y, mu, sigma, log=TRUE)
    InvL <- 1 / exp(LL)
    LL <- sum(LL)
    ### Log-Posterior
    LP <- LL + beta.prior + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,InvL),
        yhat=rnorm(length(mu), mu, sigma), parm=parm)
    return(Modelout)
    }
```

#### Initial values

``` r

Initial.Values <- c(rep(0,J), 1)
```

#### Fitting and recovery

``` r

#### Algorithm recommendation
rx <- Prescribe(Model, Data, Initial.Values)
print(rx)

#### fit
### Run the recommended code from Prescribe (see print output above)
eval(parse(text = rx$Recommend$primary$code))
print(fit)
Consort(fit)

#### Compute CPO and identify outliers
ICPO <- colMeans(fit$Monitor[, -1])  # posterior mean of InvL
CPO  <- 1 / ICPO
cat("Records with ICPO > 40 (possible outliers):", which(ICPO > 40), "\n")
cat("Records with ICPO > 70 (extreme values):   ", which(ICPO > 70), "\n")

#### Posterior predictive checks
ground_truth <- setNames(rep(NA_real_, length(parm.names)),
                         parm.names)
ground_truth[pos.beta]  <- true.beta[seq_along(pos.beta)]
ground_truth[pos.sigma] <- true.sigma[seq_along(pos.sigma)]
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

------------------------------------------------------------------------

## References

**\[1\]** Albert, J.H. (1997). Bayesian testing and estimation of
association in a two-way contingency table. *Journal of the American
Statistical Association*, 92(438), 685-693. DOI:
[10.1080/01621459.1997.10474020](https://doi.org/10.1080/01621459.1997.10474020).

**\[2\]** Albert, J.H. and Chib, S. (1993). Bayesian analysis of binary
and polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669-679. DOI:
[10.1080/01621459.1993.10476321](https://doi.org/10.1080/01621459.1993.10476321).

**\[3\]** Banerjee, S., Gelfand, A.E., Finley, A.O. and Sang, H. (2008).
Gaussian predictive process models for large spatial data sets. *Journal
of the Royal Statistical Society: Series B*, 70(4), 825-848. DOI:
[10.1111/j.1467-9868.2008.00663.x](https://doi.org/10.1111/j.1467-9868.2008.00663.x).

**\[4\]** Banfield, J.D. and Raftery, A.E. (1993). Model-based Gaussian
and non-Gaussian clustering. *Biometrics*, 49(3), 803-821. DOI:
[10.2307/2532201](https://doi.org/10.2307/2532201).

**\[5\]** Baum, L.E., Petrie, T., Soules, G. and Weiss, N. (1970). A
maximization technique occurring in the statistical analysis of
probabilistic functions of Markov chains. *Annals of Mathematical
Statistics*, 41(1), 164-171. DOI:
[10.1214/aoms/1177697196](https://doi.org/10.1214/aoms/1177697196).

**\[6\]** Beaumont, M.A., Zhang, W. and Balding, D.J. (2002).
Approximate Bayesian computation in population genetics. *Genetics*,
162(4), 2025-2035. DOI:
[10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025).

**\[7\]** Ben-Akiva, M. (1973). *Structure of Passenger Travel Demand
Models*. PhD thesis, MIT.

**\[8\]** Berkson, J. (1944). Application of the logistic function to
bio-assay. *Journal of the American Statistical Association*, 39(227),
357-365. DOI:
[10.1080/01621459.1944.10500699](https://doi.org/10.1080/01621459.1944.10500699).

**\[9\]** Bernardo, J.M. and Smith, A.F.M. (1994). *Bayesian Theory*.
Wiley. ISBN:
[978-0471494645](https://www.wiley.com/en-us/Bayesian+Theory-p-9780471494645).

**\[10\]** Besag, J. (1974). Spatial interaction and the statistical
analysis of lattice systems. *Journal of the Royal Statistical Society:
Series B*, 36(2), 192-236. DOI:
[10.1111/j.2517-6161.1974.tb00999.x](https://doi.org/10.1111/j.2517-6161.1974.tb00999.x).

**\[11\]** Bhat, C.R. (2005). A multiple discrete-continuous extreme
value model: formulation and application to discretionary time-use
decisions. *Transportation Research Part B*, 39(8), 679-707. DOI:
[10.1016/j.trb.2004.08.003](https://doi.org/10.1016/j.trb.2004.08.003).

**\[12\]** Bhattacharya, A. and Dunson, D.B. (2011). Sparse Bayesian
infinite factor models. *Biometrika*, 98(2), 291-306. DOI:
[10.1093/biomet/asr013](https://doi.org/10.1093/biomet/asr013).

**\[13\]** Blei, D.M., Ng, A.Y. and Jordan, M.I. (2003). Latent
Dirichlet allocation. *Journal of Machine Learning Research*, 3,
993-1022. [JMLR](https://jmlr.org/papers/v3/blei03a.html).

**\[14\]** Bliss, C.I. (1934). The method of probits. *Science*,
79(2037), 38-39. DOI:
[10.1126/science.79.2037.38](https://doi.org/10.1126/science.79.2037.38).

**\[15\]** Bollerslev, T. (1986). Generalized autoregressive conditional
heteroskedasticity. *Journal of Econometrics*, 31(3), 307-327. DOI:
[10.1016/0304-4076(86)90063-1](https://doi.org/10.1016/0304-4076(86)90063-1).

**\[16\]** Box, G.E.P. and Jenkins, G.M. (1970). *Time Series Analysis:
Forecasting and Control*. Holden-Day. ISBN: 978-0816211043.

**\[17\]** Brunsdon, C., Fotheringham, A.S. and Charlton, M.E. (1996).
Geographically weighted regression: a method for exploring spatial
nonstationarity. *Geographical Analysis*, 28(4), 281-298. DOI:
[10.1111/j.1538-4632.1996.tb00936.x](https://doi.org/10.1111/j.1538-4632.1996.tb00936.x).

**\[18\]** Canova, F. and Ciccarelli, M. (2013). Panel vector
autoregressive models: a survey. In *Advances in Econometrics*, 32,
205-246. DOI:
[10.1108/S0731-9053(2013)0000031006](https://doi.org/10.1108/S0731-9053(2013)0000031006).

**\[19\]** Carlin, B.P., Gelfand, A.E. and Smith, A.F.M. (1992).
Hierarchical Bayesian analysis of changepoint problems. *Journal of the
Royal Statistical Society: Series C*, 41(2), 389-405. DOI:
[10.2307/2347570](https://doi.org/10.2307/2347570).

**\[20\]** Carvalho, C.M., Chang, J., Lucas, J.E., Nevins, J.R., Wang,
Q. and West, M. (2008). High-dimensional sparse factor modeling:
applications in gene expression genomics. *Journal of the American
Statistical Association*, 103(484), 1438-1456. DOI:
[10.1198/016214508000000869](https://doi.org/10.1198/016214508000000869).

**\[21\]** Carvalho, C.M., Polson, N.G. and Scott, J.G. (2010). The
horseshoe estimator for sparse signals. *Biometrika*, 97(2), 465-480.
DOI: [10.1093/biomet/asq017](https://doi.org/10.1093/biomet/asq017).

**\[22\]** Chen, R. and Tsay, R.S. (1993). Functional-coefficient
autoregressive models. *Journal of the American Statistical
Association*, 88(421), 298-308. DOI:
[10.1080/01621459.1993.10594322](https://doi.org/10.1080/01621459.1993.10594322).

**\[23\]** Chhikara, R.S. and Folks, J.L. (1989). *The Inverse Gaussian
Distribution: Theory, Methodology, and Applications*. Marcel Dekker.
ISBN: 978-0824779979.

**\[24\]** Chib, S. (1998). Estimation and comparison of multiple
change-point models. *Journal of Econometrics*, 86(2), 221-241. DOI:
[10.1016/S0304-4076(97)00115-2](https://doi.org/10.1016/S0304-4076(97)00115-2).

**\[25\]** Chib, S. and Greenberg, E. (1998). Analysis of multivariate
probit models. *Biometrika*, 85(2), 347-361. DOI:
[10.1093/biomet/85.2.347](https://doi.org/10.1093/biomet/85.2.347).

**\[26\]** Chib, S. and Winkelmann, R. (2001). Markov chain Monte Carlo
analysis of correlated count data. *Journal of Business and Economic
Statistics*, 19(4), 428-435. DOI:
[10.1198/07350010152596673](https://doi.org/10.1198/07350010152596673).

**\[27\]** Cox, D.R. (1972). Regression models and life-tables. *Journal
of the Royal Statistical Society: Series B*, 34(2), 187-220. DOI:
[10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x).

**\[28\]** Cressie, N. and Wikle, C.K. (2011). *Statistics for
Spatio-Temporal Data*. Wiley. ISBN:
[978-0471692744](https://www.wiley.com/en-us/Statistics+for+Spatio+Temporal+Data-p-9780471692744).

**\[29\]** Denison, D.G.T., Holmes, C.C., Mallick, B.K. and Smith,
A.F.M. (2002). *Bayesian Methods for Nonlinear Classification and
Regression*. Wiley. ISBN:
[978-0471490364](https://www.wiley.com/en-us/Bayesian+Methods+for+Nonlinear+Classification+and+Regression-p-9780471490364).

**\[30\]** Diebolt, J. and Robert, C.P. (1994). Estimation of finite
mixture distributions through Bayesian sampling. *Journal of the Royal
Statistical Society: Series B*, 56(2), 363-375. DOI:
[10.1111/j.2517-6161.1994.tb01985.x](https://doi.org/10.1111/j.2517-6161.1994.tb01985.x).

**\[31\]** Doan, T., Litterman, R. and Sims, C. (1984). Forecasting and
conditional projection using realistic prior distributions. *Econometric
Reviews*, 3(1), 1-100. DOI:
[10.1080/07474938408800053](https://doi.org/10.1080/07474938408800053).

**\[32\]** Draper, D. (1995). Assessment and propagation of model
uncertainty. *Journal of the Royal Statistical Society: Series B*,
57(1), 45-97. DOI:
[10.1111/j.2517-6161.1995.tb02015.x](https://doi.org/10.1111/j.2517-6161.1995.tb02015.x).

**\[33\]** Durbin, J. and Koopman, S.J. (2012). *Time Series Analysis by
State Space Methods*, 2nd ed. Oxford University Press. ISBN:
[978-0199641178](https://global.oup.com/academic/product/time-series-analysis-by-state-space-methods-9780199641178).

**\[34\]** Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89-121. DOI:
[10.1214/ss/1038425655](https://doi.org/10.1214/ss/1038425655).

**\[35\]** Engle, R.F. (1982). Autoregressive conditional
heteroscedasticity with estimates of the variance of United Kingdom
inflation. *Econometrica*, 50(4), 987-1007. DOI:
[10.2307/1912773](https://doi.org/10.2307/1912773).

**\[36\]** Engle, R.F. and Kroner, K.F. (1995). Multivariate
simultaneous generalized ARCH. *Econometric Theory*, 11(1), 122-150.
DOI:
[10.1017/S0266466600009063](https://doi.org/10.1017/S0266466600009063).

**\[37\]** Engle, R.F., Lilien, D.M. and Robins, R.P. (1987). Estimating
time varying risk premia in the term structure: the ARCH-M model.
*Econometrica*, 55(2), 391-407. DOI:
[10.2307/1913242](https://doi.org/10.2307/1913242).

**\[38\]** Escobar, M.D. and West, M. (1995). Bayesian density
estimation and inference using mixtures. *Journal of the American
Statistical Association*, 90(430), 577-588. DOI:
[10.1080/01621459.1995.10476550](https://doi.org/10.1080/01621459.1995.10476550).

**\[39\]** Ferguson, T.S. (1973). A Bayesian analysis of some
nonparametric problems. *Annals of Statistics*, 1(2), 209-230. DOI:
[10.1214/aos/1176342360](https://doi.org/10.1214/aos/1176342360).

**\[40\]** Ferrari, S.L.P. and Cribari-Neto, F. (2004). Beta regression
for modelling rates and proportions. *Journal of Applied Statistics*,
31(7), 799-815. DOI:
[10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501).

**\[41\]** Fisher, R.A. (1925). *Statistical Methods for Research
Workers*. Oliver and Boyd.

**\[42\]** Frühwirth-Schnatter, S. (2006). *Finite Mixture and Markov
Switching Models*. Springer. DOI:
[10.1007/978-0-387-35768-3](https://doi.org/10.1007/978-0-387-35768-3).

**\[43\]** Gardner, E.S. (2006). Exponential smoothing: the state of the
art, Part II. *International Journal of Forecasting*, 22(4), 637-666.
DOI:
[10.1016/j.ijforecast.2006.03.005](https://doi.org/10.1016/j.ijforecast.2006.03.005).

**\[44\]** Gauss, C.F. (1809). *Theoria motus corporum coelestium*.
Perthes et Besser.

**\[45\]** Geisser, S. (1993). *Predictive Inference: An Introduction*.
Chapman and Hall. ISBN:
[978-0412034718](https://doi.org/10.1007/978-1-4899-4467-2).

**\[46\]** Gelfand, A.E., Dey, D.K. and Chang, H. (1992). Model
determination using predictive distributions with implementation via
sampling-based methods. In *Bayesian Statistics 4*, pp. 147-167. Oxford
University Press.

**\[47\]** Gelman, A. and Hill, J. (2006). *Data Analysis Using
Regression and Multilevel/Hierarchical Models*. Cambridge University
Press. ISBN: [978-0521686891](https://doi.org/10.1017/CBO9780511790942).

**\[48\]** Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari,
A. and Rubin, D.B. (2013). *Bayesian Data Analysis*, 3rd ed. Chapman and
Hall/CRC. ISBN:
[978-1439840955](https://www.routledge.com/Bayesian-Data-Analysis/Gelman-Carlin-Stern-Dunson-Vehtari-Rubin/p/book/9781439840955).

**\[49\]** George, E.I. and McCulloch, R.E. (1993). Variable selection
via Gibbs sampling. *Journal of the American Statistical Association*,
88(423), 881-889. DOI:
[10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353).

**\[50\]** George, E.I., Sun, D. and Ni, S. (2008). Bayesian stochastic
search for VAR model restrictions. *Journal of Econometrics*, 142(1),
553-580. DOI:
[10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017).

**\[51\]** Glosten, L.R., Jagannathan, R. and Runkle, D.E. (1993). On
the relation between the expected value and the volatility of the
nominal excess return on stocks. *Journal of Finance*, 48(5), 1779-1801.
DOI:
[10.1111/j.1540-6261.1993.tb05128.x](https://doi.org/10.1111/j.1540-6261.1993.tb05128.x).

**\[52\]** Gneiting, T. (2002). Nonseparable, stationary covariance
functions for space-time data. *Journal of the American Statistical
Association*, 97(458), 590-600. DOI:
[10.1198/016214502760047113](https://doi.org/10.1198/016214502760047113).

**\[53\]** Green, P.J. (1995). Reversible jump Markov chain Monte Carlo
computation and Bayesian model determination. *Biometrika*, 82(4),
711-732. DOI:
[10.1093/biomet/82.4.711](https://doi.org/10.1093/biomet/82.4.711).

**\[54\]** Haavelmo, T. (1943). The statistical implications of a system
of simultaneous equations. *Econometrica*, 11(1), 1-12. DOI:
[10.2307/1905714](https://doi.org/10.2307/1905714).

**\[55\]** Harvey, A.C. (1989). *Forecasting, Structural Time Series
Models and the Kalman Filter*. Cambridge University Press. ISBN:
[978-0521405737](https://doi.org/10.1017/CBO9781107049994).

**\[56\]** Hilbe, J.M. (2011). *Negative Binomial Regression*, 2nd
ed. Cambridge University Press. ISBN:
[978-0521198158](https://doi.org/10.1017/CBO9780511973420).

**\[57\]** Hoerl, A.E. and Kennard, R.W. (1970). Ridge regression:
biased estimation for nonorthogonal problems. *Technometrics*, 12(1),
55-67. DOI:
[10.1080/00401706.1970.10488634](https://doi.org/10.1080/00401706.1970.10488634).

**\[58\]** Hyndman, R.J., Koehler, A.B., Snyder, R.D. and Grose, S.
(2002). A state space framework for automatic forecasting using
exponential smoothing methods. *International Journal of Forecasting*,
18(3), 439-454. DOI:
[10.1016/S0169-2070(01)00110-8](https://doi.org/10.1016/S0169-2070(01)00110-8).

**\[59\]** Ibrahim, J.G. and Chen, M.-H. (2000). Power prior
distributions for regression models. *Statistical Science*, 15(1),
46-60. DOI:
[10.1214/ss/1009212673](https://doi.org/10.1214/ss/1009212673).

**\[60\]** Ibrahim, J.G., Chen, M.-H. and Sinha, D. (2001). *Bayesian
Survival Analysis*. Springer. DOI:
[10.1007/978-1-4757-3447-8](https://doi.org/10.1007/978-1-4757-3447-8).

**\[61\]** Jöreskog, K.G. (1969). A general approach to confirmatory
maximum likelihood factor analysis. *Psychometrika*, 34(2), 183-202.
DOI: [10.1007/BF02289343](https://doi.org/10.1007/BF02289343).

**\[62\]** Kim, J., Allenby, G.M. and Rossi, P.E. (2002). Modeling
consumer demand for variety. *Marketing Science*, 21(3), 229-250. DOI:
[10.1287/mksc.21.3.229.143](https://doi.org/10.1287/mksc.21.3.229.143).

**\[63\]** Kim, S., Shephard, N. and Chib, S. (1998). Stochastic
volatility: likelihood inference and comparison with ARCH models.
*Review of Economic Studies*, 65(3), 361-393. DOI:
[10.1111/1467-937X.00050](https://doi.org/10.1111/1467-937X.00050).

**\[64\]** Klein, L.R. (1950). *Economic Fluctuations in the United
States, 1921-1941*. Wiley.

**\[65\]** Kotz, S., Kozubowski, T.J. and Podgórski, K. (2001). *The
Laplace Distribution and Generalizations*. Birkhäuser. DOI:
[10.1007/978-1-4612-0173-1](https://doi.org/10.1007/978-1-4612-0173-1).

**\[66\]** Koyck, L.M. (1954). *Distributed Lags and Investment
Analysis*. North-Holland.

**\[67\]** Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods
for Bayesian quantile regression. *Journal of Statistical Computation
and Simulation*, 81(11), 1565-1578. DOI:
[10.1080/00949655.2010.496117](https://doi.org/10.1080/00949655.2010.496117).

**\[68\]** Krige, D.G. (1951). A statistical approach to some basic mine
valuation problems on the Witwatersrand. *Journal of the Southern
African Institute of Mining and Metallurgy*, 52(6), 119-139.

**\[69\]** Lambert, D. (1992). Zero-inflated Poisson regression, with an
application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
DOI: [10.2307/1269547](https://doi.org/10.2307/1269547).

**\[70\]** Lange, K.L., Little, R.J.A. and Taylor, J.M.G. (1989). Robust
statistical modeling using the t distribution. *Journal of the American
Statistical Association*, 84(408), 881-896. DOI:
[10.1080/01621459.1989.10478852](https://doi.org/10.1080/01621459.1989.10478852).

**\[71\]** Leng, C., Tran, M.-N. and Nott, D. (2014). Bayesian adaptive
Lasso. *Annals of the Institute of Statistical Mathematics*, 66(2),
221-244. DOI:
[10.1007/s10463-013-0429-6](https://doi.org/10.1007/s10463-013-0429-6).

**\[72\]** Lindley, D.V. and Smith, A.F.M. (1972). Bayes estimates for
the linear model. *Journal of the Royal Statistical Society: Series B*,
34(1), 1-41. DOI:
[10.1111/j.2517-6161.1972.tb00885.x](https://doi.org/10.1111/j.2517-6161.1972.tb00885.x).

**\[73\]** Litterman, R.B. (1986). Forecasting with Bayesian vector
autoregressions: five years of experience. *Journal of Business and
Economic Statistics*, 4(1), 25-38. DOI:
[10.1080/07350015.1986.10509491](https://doi.org/10.1080/07350015.1986.10509491).

**\[74\]** Little, R.J.A. and Rubin, D.B. (2002). *Statistical Analysis
with Missing Data*, 2nd ed. Wiley. ISBN:
[978-0471183860](https://doi.org/10.1002/9781119013563).

**\[75\]** Liu, C. (2004). Robit regression: a simple robust alternative
to logistic and probit regression. In *Applied Bayesian Modeling and
Causal Inference from Incomplete-Data Perspectives*, pp. 227-238. Wiley.

**\[76\]** Lopes, H.F. and West, M. (2004). Bayesian model assessment in
factor analysis. *Statistica Sinica*, 14(1), 41-67.

**\[77\]** Lütkepohl, H. (2005). *New Introduction to Multiple Time
Series Analysis*. Springer. DOI:
[10.1007/978-3-540-27752-1](https://doi.org/10.1007/978-3-540-27752-1).

**\[78\]** Matheron, G. (1963). Principles of geostatistics. *Economic
Geology*, 58(8), 1246-1266. DOI:
[10.2113/gsecongeo.58.8.1246](https://doi.org/10.2113/gsecongeo.58.8.1246).

**\[79\]** McCullagh, P. (1980). Regression models for ordinal data.
*Journal of the Royal Statistical Society: Series B*, 42(2), 109-142.
DOI:
[10.1111/j.2517-6161.1980.tb01109.x](https://doi.org/10.1111/j.2517-6161.1980.tb01109.x).

**\[80\]** McCullagh, P. and Nelder, J.A. (1989). *Generalized Linear
Models*, 2nd ed. Chapman and Hall. DOI:
[10.1007/978-1-4899-3242-6](https://doi.org/10.1007/978-1-4899-3242-6).

**\[81\]** McCulloch, R. and Rossi, P.E. (1994). An exact likelihood
analysis of the multinomial probit model. *Journal of Econometrics*,
64(1-2), 207-240. DOI:
[10.1016/0304-4076(94)90064-7](https://doi.org/10.1016/0304-4076(94)90064-7).

**\[82\]** McFadden, D. (1974). Conditional logit analysis of
qualitative choice behavior. In *Frontiers in Econometrics*,
pp. 105-142. Academic Press.

**\[83\]** McFadden, D. (1978). Modeling the choice of residential
location. In *Spatial Interaction Theory and Planning Models*,
pp. 75-96. North-Holland.

**\[84\]** McFadden, D. and Train, K. (2000). Mixed MNL models for
discrete response. *Journal of Applied Econometrics*, 15(5), 447-470.
DOI:
[10.1002/1099-1255(200009/10)15:5\<447::AID-JAE570\>3.0.CO;2-1](https://doi.org/10.1002/1099-1255(200009/10)15:5%3C447::AID-JAE570%3E3.0.CO;2-1).

**\[85\]** McLachlan, G.J. and Peel, D. (2000). *Finite Mixture Models*.
Wiley. DOI: [10.1002/0471721182](https://doi.org/10.1002/0471721182).

**\[86\]** Muth, J.F. (1960). Optimal properties of exponentially
weighted forecasts. *Journal of the American Statistical Association*,
55(290), 299-306. DOI:
[10.1080/01621459.1960.10482064](https://doi.org/10.1080/01621459.1960.10482064).

**\[87\]** Muthén, B. (1984). A general structural equation model with
dichotomous, ordered categorical, and continuous latent variable
indicators. *Psychometrika*, 49(1), 115-132. DOI:
[10.1007/BF02294210](https://doi.org/10.1007/BF02294210).

**\[88\]** Muthén, B. and Asparouhov, T. (2012). Bayesian structural
equation modeling: a more flexible representation of substantive theory.
*Psychological Methods*, 17(3), 313-335. DOI:
[10.1037/a0026802](https://doi.org/10.1037/a0026802).

**\[89\]** Nelder, J.A. and Wedderburn, R.W.M. (1972). Generalized
linear models. *Journal of the Royal Statistical Society: Series A*,
135(3), 370-384. DOI:
[10.2307/2344614](https://doi.org/10.2307/2344614).

**\[90\]** O’Hara, R.B. and Sillanpää, M.J. (2009). A review of Bayesian
variable selection methods: what, how and which. *Bayesian Analysis*,
4(1), 85-117. DOI: [10.1214/09-BA403](https://doi.org/10.1214/09-BA403).

**\[91\]** Ord, K. (1975). Estimation methods for models of spatial
interaction. *Journal of the American Statistical Association*, 70(349),
120-126. DOI:
[10.1080/01621459.1975.10480272](https://doi.org/10.1080/01621459.1975.10480272).

**\[92\]** Park, T. and Casella, G. (2008). The Bayesian Lasso. *Journal
of the American Statistical Association*, 103(482), 681-686. DOI:
[10.1198/016214508000000337](https://doi.org/10.1198/016214508000000337).

**\[93\]** Pfeifer, P.E. and Deutsch, S.J. (1980). A three-stage
iterative procedure for space-time modeling. *Technometrics*, 22(1),
35-47. DOI:
[10.1080/00401706.1980.10486099](https://doi.org/10.1080/00401706.1980.10486099).

**\[94\]** Pregibon, D. (1980). Goodness of link tests for generalized
linear models. *Journal of the Royal Statistical Society: Series C*,
29(1), 15-24. DOI: [10.2307/2346405](https://doi.org/10.2307/2346405).

**\[95\]** Pritchard, J.K., Seielstad, M.T., Perez-Lezaun, A. and
Feldman, M.W. (1999). Population growth of human Y chromosomes: a study
of Y chromosome microsatellites. *Molecular Biology and Evolution*,
16(12), 1791-1798. DOI:
[10.1093/oxfordjournals.molbev.a026091](https://doi.org/10.1093/oxfordjournals.molbev.a026091).

**\[96\]** Rubin, D.B. and Schenker, N. (1986). Multiple imputation for
interval estimation from simple random samples with ignorable
nonresponse. *Journal of the American Statistical Association*, 81(394),
366-374. DOI:
[10.1080/01621459.1986.10478280](https://doi.org/10.1080/01621459.1986.10478280).

**\[97\]** Ruppert, D., Wand, M.P. and Carroll, R.J. (2003).
*Semiparametric Regression*. Cambridge University Press. ISBN:
[978-0521785167](https://doi.org/10.1017/CBO9780511755453).

**\[98\]** Scott, S.L. (2002). Bayesian methods for hidden Markov
models: recursive computing in the 21st century. *Journal of the
American Statistical Association*, 97(457), 337-351. DOI:
[10.1198/016214502753479464](https://doi.org/10.1198/016214502753479464).

**\[99\]** Skellam, J.G. (1948). A probability distribution derived from
the binomial distribution by regarding the probability of success as
variable between the sets of trials. *Journal of the Royal Statistical
Society: Series B*, 10(2), 257-261.

**\[100\]** Spearman, C. (1904). “General intelligence,” objectively
determined and measured. *American Journal of Psychology*, 15(2),
201-292. DOI: [10.2307/1412107](https://doi.org/10.2307/1412107).

**\[101\]** Stock, J.H. and Watson, M.W. (2002). Forecasting using
principal components from a large number of predictors. *Journal of the
American Statistical Association*, 97(460), 1167-1179. DOI:
[10.1198/016214502388618960](https://doi.org/10.1198/016214502388618960).

**\[102\]** Tanner, M.A. and Wong, W.H. (1987). The calculation of
posterior distributions by data augmentation. *Journal of the American
Statistical Association*, 82(398), 528-540. DOI:
[10.1080/01621459.1987.10478458](https://doi.org/10.1080/01621459.1987.10478458).

**\[103\]** Tavaré, S., Balding, D.J., Griffiths, R.C. and Donnelly, P.
(1997). Inferring coalescence times from DNA sequence data. *Genetics*,
145(2), 505-518. DOI:
[10.1093/genetics/145.2.505](https://doi.org/10.1093/genetics/145.2.505).

**\[104\]** Taylor, S.J. (1982). Financial returns modelled by the
product of two stochastic processes: a study of the daily sugar prices
1961-75. In *Time Series Analysis: Theory and Practice 1*, pp. 203-226.
North-Holland.

**\[105\]** Teräsvirta, T. (1994). Specification, estimation, and
evaluation of smooth transition autoregressive models. *Journal of the
American Statistical Association*, 89(425), 208-218. DOI:
[10.1080/01621459.1994.10476462](https://doi.org/10.1080/01621459.1994.10476462).

**\[106\]** Thurstone, L.L. (1947). *Multiple-Factor Analysis*.
University of Chicago Press.

**\[107\]** Tong, H. (1978). On a threshold model. In *Pattern
Recognition and Signal Processing*, pp. 575-586. Sijthoff and Noordhoff.

**\[108\]** Tong, H. and Lim, K.S. (1980). Threshold autoregression,
limit cycles and cyclical data. *Journal of the Royal Statistical
Society: Series B*, 42(3), 245-292. DOI:
[10.1111/j.2517-6161.1980.tb01126.x](https://doi.org/10.1111/j.2517-6161.1980.tb01126.x).

**\[109\]** Train, K.E. (2009). *Discrete Choice Methods with
Simulation*, 2nd ed. Cambridge University Press. ISBN:
[978-0521766555](https://doi.org/10.1017/CBO9780511805271).

**\[110\]** Tweedie, M.C.K. (1957). Statistical properties of inverse
Gaussian distributions. *Annals of Mathematical Statistics*, 28(2),
362-377. DOI:
[10.1214/aoms/1177706964](https://doi.org/10.1214/aoms/1177706964).

**\[111\]** West, M. (2003). Bayesian factor regression models in the
“large p, small n” paradigm. In *Bayesian Statistics 7*, pp. 723-732.
Oxford University Press.

**\[112\]** West, M. and Harrison, J. (1997). *Bayesian Forecasting and
Dynamic Models*, 2nd ed. Springer. DOI:
[10.1007/b98971](https://doi.org/10.1007/b98971).

**\[113\]** Whittle, P. (1954). On stationary processes in the plane.
*Biometrika*, 41(3-4), 434-449. DOI:
[10.1093/biomet/41.3-4.434](https://doi.org/10.1093/biomet/41.3-4.434).

**\[114\]** Wilks, S.S. (1932). Certain generalizations in the analysis
of variance. *Biometrika*, 24(3-4), 471-494. DOI:
[10.1093/biomet/24.3-4.471](https://doi.org/10.1093/biomet/24.3-4.471).

**\[115\]** Williams, D.A. (1975). The analysis of binary responses from
toxicological experiments involving reproduction and teratogenicity.
*Biometrics*, 31(4), 949-952. DOI:
[10.2307/2529820](https://doi.org/10.2307/2529820).

**\[116\]** Yu, K. and Moyeed, R.A. (2001). Bayesian quantile
regression. *Statistics and Probability Letters*, 54(4), 437-447. DOI:
[10.1016/S0167-7152(01)00124-9](https://doi.org/10.1016/S0167-7152(01)00124-9).

**\[117\]** Yule, G.U. (1927). On a method of investigating
periodicities in disturbed series, with special reference to Wolfer’s
sunspot numbers. *Philosophical Transactions of the Royal Society of
London A*, 226, 267-298. DOI:
[10.1098/rsta.1927.0007](https://doi.org/10.1098/rsta.1927.0007).

**\[118\]** Zakoian, J.-M. (1994). Threshold heteroskedastic models.
*Journal of Economic Dynamics and Control*, 18(5), 931-955. DOI:
[10.1016/0165-1889(94)90039-6](https://doi.org/10.1016/0165-1889(94)90039-6).

**\[119\]** Zeger, S.L. (1988). A regression model for time series of
counts. *Biometrika*, 75(4), 621-629. DOI:
[10.1093/biomet/75.4.621](https://doi.org/10.1093/biomet/75.4.621).

**\[120\]** Zellner, A. (1962). An efficient method of estimating
seemingly unrelated regressions and tests for aggregation bias. *Journal
of the American Statistical Association*, 57(298), 348-368. DOI:
[10.1080/01621459.1962.10480664](https://doi.org/10.1080/01621459.1962.10480664).

**\[121\]** Zellner, A. (1971). *An Introduction to Bayesian Inference
in Econometrics*. Wiley. ISBN: 978-0471169376.

**\[122\]** Zellner, A. (1986). On assessing prior distributions and
Bayesian regression analysis with g-prior distributions. In *Bayesian
Inference and Decision Techniques*, pp. 233-243. North-Holland.

**\[123\]** Zou, H. (2006). The adaptive lasso and its oracle
properties. *Journal of the American Statistical Association*, 101(476),
1418-1429. DOI:
[10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735).

**\[124\]** Einstein, A. (1905). Uber die von der molekularkinetischen
Theorie der Warme geforderte Bewegung von in ruhenden Flussigkeiten
suspendierten Teilchen. *Annalen der Physik*, 17(8), 549-560. DOI:
[10.1002/andp.19053220806](https://doi.org/10.1002/andp.19053220806).

**\[125\]** Pearson, K. (1905). The problem of the random walk.
*Nature*, 72(1865), 294. DOI:
[10.1038/072294b0](https://doi.org/10.1038/072294b0).

**\[126\]** Fuller, W.A. (1987). *Measurement Error Models*. Wiley.
ISBN:
[978-0471861874](https://www.wiley.com/en-us/Measurement+Error+Models-p-9780471861874).

**\[127\]** Carroll, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu,
C.M. (2006). *Measurement Error in Nonlinear Models: A Modern
Perspective*, 2nd ed. Chapman and Hall/CRC. DOI:
[10.1201/9781420010138](https://doi.org/10.1201/9781420010138).

**\[128\]** Lawless, J.F. (2003). *Statistical Models and Methods for
Lifetime Data*, 2nd ed. Wiley. DOI:
[10.1002/9781118033005](https://doi.org/10.1002/9781118033005).

**\[129\]** Schuster, A. (1898). On the investigation of hidden
periodicities with application to a suggested investigation of the
period of 26 days. *Terrestrial Magnetism*, 3(1), 13-41. DOI:
[10.1029/TM003i001p00013](https://doi.org/10.1029/TM003i001p00013).

**\[130\]** Jaynes, E.T. (1987). Bayesian spectrum and chirp analysis.
In *Maximum Entropy and Bayesian Spectral Analysis and Estimation
Problems*, pp. 1-37. Springer. DOI:
[10.1007/978-94-009-3961-5_1](https://doi.org/10.1007/978-94-009-3961-5_1).

**\[131\]** Gompertz, B. (1825). On the nature of the function
expressive of the law of human mortality, and on a new mode of
determining the value of life contingencies. *Philosophical Transactions
of the Royal Society of London*, 115, 513-583.

**\[132\]** Laird, A.K. (1964). Dynamics of tumour growth. *British
Journal of Cancer*, 18(3), 490-502. DOI:
[10.1038/bjc.1964.55](https://doi.org/10.1038/bjc.1964.55).

**\[133\]** Bernal, J.L., Cummins, S. and Gasparrini, A. (2017).
Interrupted time series regression for the evaluation of public health
interventions: a tutorial. *International Journal of Epidemiology*,
46(1), 348-355. DOI:
[10.1093/ije/dyw098](https://doi.org/10.1093/ije/dyw098).

**\[134\]** Wagner, A.K., Soumerai, S.B., Zhang, F. and Ross-Degnan, D.
(2002). Segmented regression analysis of interrupted time series studies
in medication use research. *Journal of Clinical Pharmacy and
Therapeutics*, 27(4), 299-309. DOI:
[10.1046/j.1365-2710.2002.00430.x](https://doi.org/10.1046/j.1365-2710.2002.00430.x).

**\[135\]** Raudenbush, S.W. and Bryk, A.S. (2002). *Hierarchical Linear
Models: Applications and Data Analysis Methods*, 2nd ed. Sage. ISBN:
[978-0761919049](https://us.sagepub.com/en-us/nam/hierarchical-linear-models/book9230).

**\[136\]** Barry, D. and Hartigan, J.A. (1993). A Bayesian analysis for
change point problems. *Journal of the American Statistical
Association*, 88(421), 309-319. DOI:
[10.1080/01621459.1993.10594323](https://doi.org/10.1080/01621459.1993.10594323).

**\[137\]** Mullahy, J. (1986). Specification and testing of some
modified count data models. *Journal of Econometrics*, 33(3), 341-365.
DOI:
[10.1016/0304-4076(86)90002-3](https://doi.org/10.1016/0304-4076(86)90002-3).

**\[138\]** Zuur, A.F., Ieno, E.N., Walker, N.J., Saveliev, A.A. and
Smith, G.M. (2009). *Mixed Effects Models and Extensions in Ecology with
R*. Springer. DOI:
[10.1007/978-0-387-87458-6](https://doi.org/10.1007/978-0-387-87458-6).

**\[139\]** Hijazi, R.H. and Jernigan, R.W. (2009). Modelling
compositional data using Dirichlet regression models. *Journal of
Applied Probability and Statistics*, 4(1), 77-91.

**\[140\]** Maier, M.J. (2014). DirichletReg: Dirichlet regression for
compositional data in R. Research Report Series, Institute for
Statistics and Mathematics, WU Vienna University of Economics and
Business.

**\[141\]** Tobin, J. (1958). Estimation of relationships for limited
dependent variables. *Econometrica*, 26(1), 24-36. DOI:
[10.2307/1907382](https://doi.org/10.2307/1907382).

**\[142\]** Chib, S. (1992). Bayes inference in the Tobit censored
regression model. *Journal of Econometrics*, 51(1-2), 79-99. DOI:
[10.1016/0304-4076(92)90030-U](https://doi.org/10.1016/0304-4076(92)90030-U).

**\[143\]** Williams, D.A. (1982). Extra-binomial variation in logistic
linear models. *Journal of the Royal Statistical Society: Series C*,
31(2), 144-148. DOI: [10.2307/2347977](https://doi.org/10.2307/2347977).

**\[144\]** Mosimann, J.E. (1962). On the compound multinomial
distribution, the multivariate beta-distribution, and correlations among
proportions. *Biometrika*, 49(1-2), 65-82. DOI:
[10.1093/biomet/49.1-2.65](https://doi.org/10.1093/biomet/49.1-2.65).

**\[145\]** Jorgensen, B. (1997). *The Theory of Dispersion Models*.
Chapman and Hall. ISBN:
[978-0412997112](https://www.routledge.com/The-Theory-of-Dispersion-Models/Jorgensen/p/book/9780412997112).

**\[146\]** Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of
Tweedie exponential dispersion model densities. *Statistics and
Computing*, 15(4), 267-280. DOI:
[10.1007/s11222-005-4070-y](https://doi.org/10.1007/s11222-005-4070-y).

**\[147\]** Uhlenbeck, G.E. and Ornstein, L.S. (1930). On the theory of
the Brownian motion. *Physical Review*, 36(5), 823-841. DOI:
[10.1103/PhysRev.36.823](https://doi.org/10.1103/PhysRev.36.823).

**\[148\]** Vasicek, O. (1977). An equilibrium characterization of the
term structure. *Journal of Financial Economics*, 5(2), 177-188. DOI:
[10.1016/0304-405X(77)90016-2](https://doi.org/10.1016/0304-405X(77)90016-2).

**\[149\]** Rocha, A.V. and Cribari-Neto, F. (2009). Beta autoregressive
moving average models. *TEST*, 18(3), 529-545. DOI:
[10.1007/s11749-008-0112-z](https://doi.org/10.1007/s11749-008-0112-z).

**\[150\]** Pankratz, A. (1991). *Forecasting with Dynamic Regression
Models*. Wiley. ISBN:
[978-0471615286](https://www.wiley.com/en-us/Forecasting+with+Dynamic+Regression+Models-p-9780471615286).

**\[151\]** Aitken, A.C. (1935). On least squares and linear combination
of observations. *Proceedings of the Royal Society of Edinburgh*, 55,
42-48. DOI:
[10.1017/S0370164600014346](https://doi.org/10.1017/S0370164600014346).

**\[152\]** Kedem, B. and Fokianos, K. (2002). *Regression Models for
Time Series Analysis*. Wiley. DOI:
[10.1002/0471266981](https://doi.org/10.1002/0471266981).

**\[153\]** Birnbaum, A. (1968). Some latent trait models and their use
in inferring an examinee’s ability. In Lord, F.M. and Novick, M.R.
(eds), *Statistical Theories of Mental Test Scores*, pp. 395-479.
Addison-Wesley.

**\[154\]** Albert, J.H. (1992). Bayesian estimation of normal ogive
item response curves using Gibbs sampling. *Journal of Educational
Statistics*, 17(3), 251-269. DOI:
[10.3102/10769986017003251](https://doi.org/10.3102/10769986017003251).

**\[155\]** Shmueli, G., Minka, T.P., Kadane, J.B., Borle, S. and
Boatwright, P. (2005). A useful distribution for fitting discrete data:
revival of the Conway-Maxwell-Poisson distribution. *Journal of the
Royal Statistical Society: Series C*, 54(1), 127-142. DOI:
[10.1111/j.1467-9876.2005.00474.x](https://doi.org/10.1111/j.1467-9876.2005.00474.x).

**\[156\]** Dawid, A.P. and Skene, A.M. (1979). Maximum likelihood
estimation of observer error-rates using the EM algorithm. *Journal of
the Royal Statistical Society: Series C*, 28(1), 20-28. DOI:
[10.2307/2346806](https://doi.org/10.2307/2346806).

**\[157\]** Lang, S. and Brezger, A. (2004). Bayesian P-splines.
*Journal of Computational and Graphical Statistics*, 13(1), 183-212.
DOI: [10.1198/1061860043010](https://doi.org/10.1198/1061860043010).

**\[158\]** Wood, S.N. (2006). Low-rank scale-invariant tensor product
smooths for generalized additive modelling. *Biometrics*, 62(4),
1025-1036. DOI:
[10.1111/j.1541-0420.2006.00574.x](https://doi.org/10.1111/j.1541-0420.2006.00574.x).

**\[159\]** Gauch, H.G. (1992). *Statistical Analysis of Regional Yield
Trials: AMMI Analysis of Factorial Designs*. Elsevier. ISBN:
[978-0444892409](https://www.elsevier.com/books/statistical-analysis-of-regional-yield-trials/gauch/978-0-444-89240-9).

**\[160\]** Aitchison, J. (1986). *The Statistical Analysis of
Compositional Data*. Chapman and Hall. ISBN:
[978-0412280603](https://doi.org/10.1007/978-94-009-4109-0).

**\[161\]** Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and
Barcelo-Vidal, C. (2003). Isometric logratio transformations for
compositional data analysis. *Mathematical Geology*, 35(3), 279-300.
DOI: [10.1023/A:1023818214614](https://doi.org/10.1023/A:1023818214614).

**\[162\]** LeSage, J.P. and Pace, R.K. (2009). *Introduction to Spatial
Econometrics*. CRC Press. ISBN:
[978-1420064247](https://www.routledge.com/Introduction-to-Spatial-Econometrics/LeSage-Pace/p/book/9781420064247).

**\[163\]** Moller, J., Syversveen, A.R. and Waagepetersen, R.P. (1998).
Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25(3),
451-482. DOI:
[10.1111/1467-9469.00115](https://doi.org/10.1111/1467-9469.00115).

**\[164\]** Tavare, S. (1986). Some probabilistic and statistical
problems in the analysis of DNA sequences. *Lectures on Mathematics in
the Life Sciences*, 17, 57-86.

**\[165\]** Yang, Z. (1994). Maximum likelihood phylogenetic estimation
from DNA sequences with variable rates over sites: approximate methods.
*Journal of Molecular Evolution*, 39(3), 306-314. DOI:
[10.1007/BF00160154](https://doi.org/10.1007/BF00160154).

**\[166\]** Felsenstein, J. (1981). Evolutionary trees from DNA
sequences: a maximum likelihood approach. *Journal of Molecular
Evolution*, 17(6), 368-376. DOI:
[10.1007/BF01734359](https://doi.org/10.1007/BF01734359).

**\[167\]** Schliep, K.P. (2011). phangorn: phylogenetic analysis in R.
*Bioinformatics*, 27(4), 592-593. DOI:
[10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706).

**\[168\]** Paradis, E. and Schliep, K. (2019). ape 5.0: an environment
for modern phylogenetics and evolutionary analyses in R.
*Bioinformatics*, 35(3), 526-528. DOI:
[10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633).

**\[169\]** Drummond, A.J., Rambaut, A., Shapiro, B. and Pybus, O.G.
(2005). Bayesian coalescent inference of past population dynamics from
molecular sequences. *Molecular Biology and Evolution*, 22(5),
1185-1192. DOI:
[10.1093/molbev/msi103](https://doi.org/10.1093/molbev/msi103).

**\[170\]** Kingman, J.F.C. (1982). The coalescent. *Stochastic
Processes and their Applications*, 13(3), 235-248. DOI:
[10.1016/0304-4149(82)90011-4](https://doi.org/10.1016/0304-4149(82)90011-4).

**\[171\]** Suchard, M.A., Lemey, P., Baele, G., Ayres, D.L., Drummond,
A.J. and Rambaut, A. (2018). Bayesian phylogenetic and phylodynamic data
integration using BEAST 1.10. *Virus Evolution*, 4(1), vey016. DOI:
[10.1093/ve/vey016](https://doi.org/10.1093/ve/vey016).

**\[172\]** Samia, N.I., Stramer, O., Saitoh, T. and Stenseth, N.C.
(2024). Climate-driven context-dependent structure of population cycles.
*Royal Society Open Science*, 11(8), 240047. DOI:
[10.1098/rsos.240047](https://doi.org/10.1098/rsos.240047).

**\[173\]** Hjellvik, V. and Tjostheim, D. (1999). Modelling panels of
intercorrelated autoregressive time series. *Biometrika*, 86(3),
573-590. DOI:
[10.1093/biomet/86.3.573](https://doi.org/10.1093/biomet/86.3.573).
