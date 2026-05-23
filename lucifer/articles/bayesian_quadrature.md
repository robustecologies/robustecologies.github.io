# Bayesian quadrature in lucifer

## Overview

Bayesian quadrature (BQ) treats numerical integration as a statistical
inference problem. Rather than returning a point estimate of an
integral, BQ models the integrand with a Gaussian process surrogate and
computes a posterior distribution *over the integral value itself*. This
paradigm, introduced by O’Hagan (1991) [\[1\]](#ref1) and developed
substantially by Rasmussen and Ghahramani (2003) [\[2\]](#ref2),
provides three advantages over classical quadrature:
uncertainty-calibrated marginal likelihood estimates for model
comparison, superior sample efficiency when model evaluations are
expensive, and principled active learning for selecting where to
evaluate the model next.

The
[`BayesianQuadrature()`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)
function in lucifer implements four modern BQ methods that share the
standard `Model(parm, Data)` interface with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
and
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).
All four methods model the unnormalized posterior density \\f(\theta) =
\exp(\mathrm{LP}(\theta))\\ with a Gaussian process equipped with a
squared exponential kernel, then exploit closed-form kernel mean
embeddings against a Gaussian reference measure to compute the integral
posterior analytically.

## Mathematical foundations

### The integration problem

The central computational task in Bayesian inference is evaluating the
marginal likelihood (or normalizing constant):

\\Z = \int \exp\bigl(\mathrm{LP}(\theta)\bigr) \\ d\theta = \int
L(\theta) \pi(\theta) \\ d\theta\\

where \\L(\theta)\\ is the likelihood and \\\pi(\theta)\\ the prior.
Classical quadrature rules like Gauss-Hermite (as implemented in
[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md))
compute \\Z\\ by deterministic weighted sums \\\hat{Z} = \sum_i w_i
f(\theta_i)\\ with no uncertainty quantification on the estimate. MCMC
avoids computing \\Z\\ altogether, but many tasks, notably Bayes factor
computation and model selection, require it.

### Gaussian process surrogate

BQ places a Gaussian process prior over the integrand \\f\\:

\\f \sim \mathcal{GP}\bigl(0, k(\theta, \theta')\bigr)\\

with squared exponential (SE) kernel featuring automatic relevance
determination (ARD):

\\k(\theta, \theta') = \sigma_f^2 \exp\left(-\frac{1}{2} \sum\_{d=1}^D
\frac{(\theta_d - \theta'\_d)^2}{\ell_d^2}\right)\\

After observing \\f\\ at points \\\Theta = \\\theta_1, \dots,
\theta_n\\\\ with values \\\mathbf{y} = (f(\theta_1), \dots,
f(\theta_n))^\top\\, the GP posterior is:

\\f \mid \mathbf{y} \sim \mathcal{GP}\bigl(\mu_n(\cdot), k_n(\cdot,
\cdot)\bigr)\\

where \\\mu_n(\theta) = \mathbf{k}(\theta)^\top K^{-1} \mathbf{y}\\ and
\\k_n(\theta, \theta') = k(\theta, \theta') - \mathbf{k}(\theta)^\top
K^{-1} \mathbf{k}(\theta')\\, with \\K\_{ij} = k(\theta_i, \theta_j)\\
and \\\mathbf{k}(\theta)\_i = k(\theta_i, \theta)\\.

### Kernel mean embedding

The integral \\Z = \int f(\theta) \pi(\theta) d\theta\\ is a linear
functional of \\f\\, so under the GP model it inherits a Gaussian
posterior:

\\Z \mid \mathbf{y} \sim \mathcal{N}\bigl(\hat{Z}, \hat{V}\bigr)\\

where \\\hat{Z} = \boldsymbol{\mu}\_k^\top K^{-1} \mathbf{y}\\ and
\\\hat{V} = z_0 - \boldsymbol{\mu}\_k^\top K^{-1} \boldsymbol{\mu}\_k\\.
The kernel mean embedding vector \\\boldsymbol{\mu}\_k\\ and prior
integral \\z_0\\ admit closed-form expressions when \\\pi(\theta) =
\mathcal{N}(\mu_0, \Sigma_0)\\ and the kernel is SE:

\\\[\boldsymbol{\mu}\_k\]\_i = \sigma_f^2 \|\Lambda^{-1}\Sigma_0 +
I\|^{-1/2} \exp\left(-\frac{1}{2}(\theta_i - \mu_0)^\top (\Lambda +
\Sigma_0)^{-1} (\theta_i - \mu_0)\right)\\

\\z_0 = \sigma_f^2 \|\Lambda^{-1} (2\Sigma_0) + I\|^{-1/2}\\

where \\\Lambda = \mathrm{diag}(\ell_1^2, \dots, \ell_D^2)\\. These are
computed in C++ with OpenMP parallelization.

## Algorithms

### BQ (vanilla Bayesian quadrature)

The baseline method of Rasmussen and Ghahramani (2003) [\[2\]](#ref2).
At each iteration, BQ generates candidate evaluation points from the
reference measure, selects the candidate that maximally reduces the
integral posterior variance \\\hat{V}\\, evaluates the model there, and
updates the GP. Hyperparameters \\(\sigma_f, \boldsymbol{\ell})\\ are
re-optimized periodically via marginal likelihood maximization with
multi-start L-BFGS-B.

``` r

library(lucifer)

# Normal-normal conjugate model
Model <- function(parm, Data) {
  mu <- parm[1]
  LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
  LP <- LL + dnorm(mu, 0, 1000, log = TRUE)
  list(LP = LP, Dev = -2 * LL,
       Monitor = LP, yhat = rep(mu, length(Data$y)),
       parm = parm)
}

set.seed(42)
y <- rnorm(50, mean = 3)
Data <- list(y = y, mon.names = "LP", parm.names = "mu")

fit_bq <- BayesianQuadrature(Model, parm = 3, Data = Data,
  Covar = matrix(0.1), Iterations = 50, Algorithm = "BQ")

print(fit_bq)
summary(fit_bq)
```

### WSABI-L (warped sequential active Bayesian integration)

Gunter et al. (2014) [\[3\]](#ref3) observed that the likelihood
integrand \\f(\theta) = \exp(\mathrm{LP}(\theta))\\ is strictly
non-negative, but standard GP models can predict negative values.
WSABI-L applies a square-root warping \\g(\theta) = \sqrt{2(f(\theta) -
\alpha)}\\ where \\\alpha \< \min f\\, fits a GP to \\g\\, and
linearizes the inverse warping \\f = \alpha + g^2/2\\ for closed-form
integral moment propagation. This enforces non-negativity and is
recommended for marginal likelihood estimation.

``` r

fit_wsabi <- BayesianQuadrature(Model, parm = 3, Data = Data,
  Covar = matrix(0.1), Iterations = 50, Algorithm = "WSABI-L")

# Compare LML estimates
cat("BQ  LML:", fit_bq$LML.BQ, "\n")
cat("WSABI LML:", fit_wsabi$LML.BQ, "\n")
```

### FWBQ (Frank-Wolfe Bayesian quadrature)

Briol et al. (2015) [\[4\]](#ref4) reformulated BQ node selection as
Frank-Wolfe optimization in a reproducing kernel Hilbert space. At each
step, the algorithm selects the point \\\theta\_{t+1}\\ maximizing the
Frank-Wolfe gradient \\g_t(\theta) = \mu_k(\theta) - \sum_i w_i
k(\theta_i, \theta)\\ and recomputes quadrature weights as \\\mathbf{w}
= K^{-1} \boldsymbol{\mu}\_k\\. FWBQ is less sensitive to GP
hyperparameters than vanilla BQ and achieves provable exponential
convergence rates under regularity conditions on the integrand. It can
start with as few as one evaluation point.

``` r

fit_fwbq <- BayesianQuadrature(Model, parm = 3, Data = Data,
  Covar = matrix(0.1), Iterations = 50, Algorithm = "FWBQ")

print(fit_fwbq)
```

### BatchBQ (batch Bayesian quadrature)

When model evaluations can be parallelized, BatchBQ selects batches of
\\B\\ points via greedy sequential conditioning: the first point
maximizes the standard acquisition function, then the GP is
conditionally updated (without model evaluation) and the next point is
selected, and so on. After the batch is assembled, all \\B\\ points are
evaluated in parallel using the `CPUs` and `Type` arguments.

``` r

fit_batch <- BayesianQuadrature(Model, parm = 3, Data = Data,
  Covar = matrix(0.1), Iterations = 20, Algorithm = "BatchBQ",
  Specs = list(batch_size = 4), CPUs = 2)
```

## Example: logistic regression

This example demonstrates BQ on a two-parameter logistic regression
model where analytical solutions are unavailable.

``` r

# Logistic regression model
LogisticModel <- function(parm, Data) {
  beta0 <- parm[1]
  beta1 <- parm[2]
  eta <- beta0 + beta1 * Data$x
  p <- 1 / (1 + exp(-eta))
  p <- pmax(pmin(p, 1 - 1e-10), 1e-10)  # numerical safety
  LL <- sum(Data$y * log(p) + (1 - Data$y) * log(1 - p))
  LP <- LL + dnorm(beta0, 0, 10, log = TRUE) + dnorm(beta1, 0, 10, log = TRUE)
  list(LP = LP, Dev = -2 * LL,
       Monitor = c(LP, beta0, beta1),
       yhat = p, parm = parm)
}

set.seed(42)
n <- 100
x <- rnorm(n)
p_true <- 1 / (1 + exp(-(1 + 2 * x)))
y_log <- rbinom(n, 1, p_true)

Data_log <- list(y = y_log, x = x,
  mon.names = c("LP", "beta0", "beta1"),
  parm.names = c("beta0", "beta1"))

# Fit with BQ (start near MLE for efficiency)
mle <- coef(glm(y_log ~ x, family = binomial))
fit_log <- BayesianQuadrature(LogisticModel, parm = as.numeric(mle),
  Data = Data_log, Iterations = 60, Algorithm = "BQ",
  sir = TRUE, Samples = 1000)

print(fit_log)
summary(fit_log)
plot(fit_log, Data = Data_log)
```

## Example: multivariate normal

A three-dimensional multivariate normal model with known covariance,
providing an analytical benchmark for the BQ integral estimate.

``` r

MVNModel <- function(parm, Data) {
  mu <- parm
  diff <- sweep(Data$Y, 2, mu)
  LL <- -0.5 * sum(diff %*% Data$Sigma_inv * diff) -
    (nrow(Data$Y) / 2) * determinant(Data$Sigma, logarithm = TRUE)$modulus
  LP <- as.numeric(LL)
  list(LP = LP, Dev = -2 * LL,
       Monitor = LP, yhat = rep(0, nrow(Data$Y)),
       parm = parm)
}

set.seed(42)
k <- 3
n <- 100
Sigma_true <- diag(k) * 0.5 + 0.5
Y <- MASS::mvrnorm(n, mu = c(1, 2, 3), Sigma = Sigma_true)

Data_mvn <- list(Y = Y, Sigma = Sigma_true,
  Sigma_inv = solve(Sigma_true),
  mon.names = "LP",
  parm.names = paste0("mu", 1:k))

# BQ with initial values near posterior mode
fit_mvn <- BayesianQuadrature(MVNModel, parm = colMeans(Y),
  Data = Data_mvn, Covar = Sigma_true / n * 4,
  Iterations = 80, Algorithm = "BQ", sir = TRUE)

summary(fit_mvn)
```

## Practical guidance

Bayesian quadrature is most effective for problems where (1) the
parameter dimension is moderate (\\d \leq 20\\), (2) the model
evaluation is expensive enough that minimizing the number of evaluations
matters, and (3) an uncertainty-calibrated integral estimate is needed,
as in Bayes factor computation. For higher-dimensional problems,
[`IterativeQuadrature()`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)
with algorithm `"CAGH"` or MCMC via
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
will be more efficient.

The `parm` argument should ideally be set near the posterior mode; for
example, from a Laplace approximation or maximum likelihood estimate. BQ
actively learns the integrand surface around this starting point, but it
cannot efficiently explore very distant regions. The `Covar` matrix
defines the Gaussian reference measure and should cover the bulk of the
posterior mass; setting it to a scaled version of the posterior
covariance (e.g., from
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md))
works well.

WSABI-L is the recommended default for marginal likelihood estimation
because it enforces non-negativity of the likelihood integrand. FWBQ is
preferable when GP hyperparameter sensitivity is a concern, since its
convergence guarantees hold under weaker assumptions. BatchBQ is useful
when multiple cores are available and individual model evaluations are
expensive.

## References

**\[1\]** O’Hagan, A. (1991). Bayes-Hermite quadrature. *Journal of the
Royal Statistical Society: Series B*, 53(1), 145-169.

**\[2\]** Rasmussen, C.E. and Ghahramani, Z. (2003). Bayesian Monte
Carlo. *Advances in Neural Information Processing Systems*, 15.

**\[3\]** Gunter, T., Osborne, M.A., Garnett, R., Hennig, P., and
Roberts, S.J. (2014). Sampling for inference in probabilistic models
with fast Bayesian quadrature. *Advances in Neural Information
Processing Systems*, 27.

**\[4\]** Briol, F.-X., Oates, C.J., Girolami, M., and Osborne,
M.A. (2015). Frank-Wolfe Bayesian quadrature: probabilistic integration
with theoretical guarantees. *Advances in Neural Information Processing
Systems*, 28.

**\[5\]** Briol, F.-X., Oates, C.J., Girolami, M., Osborne, M.A., and
Sejdinovic, D. (2019). Probabilistic integration: a role in statistical
computation? *Statistical Science*, 34(1), 1-22.
[doi:10.1214/18-STS660](https://doi.org/10.1214/18-STS660)

**\[6\]** Rasmussen, C.E. and Williams, C.K.I. (2006). *Gaussian
Processes for Machine Learning*. MIT Press. ISBN 0-262-18253-X.
