# Population Monte Carlo

The `PMC` function updates a model with Population Monte Carlo. Given a
model specification, data, and initial values, `PMC` maximizes the
logarithm of the unnormalized joint posterior density and provides
samples of the marginal posterior distributions, deviance, and other
monitored variables.

## Usage

``` r
PMC(
  Model,
  Data,
  Initial.Values,
  Covar = NULL,
  Iterations = 10,
  Thinning = 1,
  alpha = NULL,
  M = 1,
  N = 1000,
  nu = 9,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  A model specification function. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  A list of data. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Initial.Values:

  Either a vector of initial values, one for each of \\K\\ parameters,
  or in the case of a mixture of \\M\\ components, a \\M \times K\\
  matrix of initial values. If all initial values are zero in this
  vector, or in the first row of a matrix, then
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
  is used to optimize initial values, in which case all mixture
  components receive the same initial values and covariance matrix from
  the object of class `laplace`. Parameters must be continuous.

- Covar:

  A \\K \times K\\ covariance matrix for \\K\\ parameters, or for
  multiple mixture components, a \\K \times K \times M\\ array of \\M\\
  covariance matrices, where \\M\\ is the number of mixture components.
  Defaults to `NULL`, in which case a scaled identity matrix (with the
  same scale as in
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md))
  is applied to all mixture components.

- Iterations:

  The number of iterations during which PMC will update the model.
  Updating the model for only one iteration is the same as applying
  non-adaptive importance sampling.

- Thinning:

  The number by which the posterior is thinned. To have 1,000 posterior
  samples with `M=3` mixture components and `N=10000` samples each,
  `Thinning=30`. For more information, see the
  [`Thin`](https://robustecologies.github.io/lucifer/reference/Thin.md)
  function.

- alpha:

  A vector of length \\M\\, the number of mixture components. \\\alpha\\
  is the probability of each mixture component. The default value is
  `NULL`, which assigns an equal probability to each component.

- M:

  The number \\M\\ of multivariate t distribution mixture components.

- N:

  The number \\N\\ of samples per mixture component. The required number
  of samples increases with the number \\K\\ of parameters. These
  samples are also called walkers or particles.

- nu:

  The degrees of freedom parameter \\\nu\\ for the multivariate t
  distribution for each mixture component. If a multivariate normal
  distribution is preferred, then set \\\nu \> 1e4\\.

- CPUs:

  The number of central processing units (CPUs) for parallel processing.
  For example, when a user has a quad-core computer, `CPUs=4`.

- Type:

  Defaults to `"PSOCK"` and uses the Simple Network of Workstations
  (SNOW) for parallelization. Alternatively, `Type="MPI"` may be
  specified to use Message Passing Interface (MPI) for parallelization.

## Value

An object of class `pmc` with the following components:

- alpha:

  A \\M \times T\\ matrix of the probabilities of mixture components.

- Call:

  The matched call of `PMC`.

- Covar:

  The \\K \times K \times T \times M\\ proposal covariance array.

- Deviance:

  A vector of the deviance of the model.

- DIC:

  A vector of three values: Dbar, pD, and DIC.

- ESSN:

  A vector of normalized effective sample sizes per iteration.

- Initial.Values:

  The vector or matrix of initial values.

- Iterations:

  The number of iterations.

- LML:

  An approximation of the logarithm of the marginal likelihood.

- M:

  The number of mixture components.

- Minutes:

  The number of minutes that `PMC` was running.

- Model:

  The model specification function.

- N:

  The number of un-thinned samples per mixture component.

- nu:

  The degrees of freedom parameter.

- Mu:

  A \\T \times K \times M\\ array of means for the importance sampling
  distribution.

- Monitor:

  A matrix of thinned samples of monitored variables.

- Parameters:

  The number \\K\\ of parameters.

- Perplexity:

  A vector of normalized perplexity per iteration.

- Posterior1:

  An \\N \times K \times T \times M\\ array of un-thinned posterior
  samples.

- Posterior2:

  A \\S \times K\\ matrix of thinned posterior samples.

- Summary:

  A matrix summarizing the marginal posterior distributions.

- Thinned.Samples:

  The number of thinned samples.

- Thinning:

  The amount of thinning requested.

- W:

  A \\N \times T\\ matrix of normalized importance weights.

## Details

The `PMC` function uses the adaptive importance sampling algorithm of
Wraith et al. (2009), also called Mixture PMC or M-PMC (Cappe et al.,
2008). Iterative adaptive importance sampling was introduced in the
1980s. Modern PMC was introduced (Cappe et al., 2004), and extended to
multivariate Gaussian or t-distributed mixtures (Cappe et al., 2008).
This version uses a multivariate t distribution for each mixture
component, and also allows a multivariate normal distribution when the
degrees of freedom, \\\nu \> 1e4\\. At each iteration, a mixture
distribution is sampled with importance sampling, and the samples (or
populations) are adapted to improve the importance sampling. Adaptation
is a variant of EM (Expectation-Maximization). The sample is
self-normalized, and is an example of self-normalized importance
sampling (SNIS), or self-importance sampling. The vector \\\alpha\\
contains the probability of each mixture component. These, as well as
multivariate t distribution mixture parameters (except \\\nu\\), are
adapted each iteration.

Advantages of PMC over MCMC include: it is difficult to assess
convergence of MCMC chains, and this is not necessary in PMC (Wraith et
al., 2009); MCMC chains have autocorrelation that effectively reduces
posterior samples, whereas PMC produces independent samples that are not
reduced with autocorrelation; PMC has been reported to produce samples
with less variance than MCMC; it is difficult to parallelize MCMC, while
PMC can parallelize the independent Monte Carlo samples during each
iteration; and the multivariate mixture in PMC can represent a
multimodal posterior.

Disadvantages of PMC, compared to MCMC, include: in PMC, the required
number of samples at each iteration increases quickly with respect to an
increase in parameters, so MCMC is more suitable for models with large
numbers of parameters; PMC is more sensitive to initial values than
MCMC, especially as the number of parameters increases; and PMC is more
sensitive to the initial covariance matrix (or matrices for mixture
components) than adaptive MCMC.

Since PMC requires better initial information than iterative quadrature,
Laplace Approximation, MCMC, and Variational Bayes, it is not
recommended to begin updating a model that has little prior information
with PMC, especially when the model has more than a few parameters.
Instead, iterative quadrature, Laplace Approximation, MCMC, or
Variational Bayes should be used. However, once convergence is found or
assumed, it is recommended to attempt to update the model with PMC,
given the latest parameters and covariance matrix from iterative
quadrature, Laplace Approximation, MCMC, or Variational Bayes.

Convergence is assessed by observing two outputs: normalized effective
sample size (`ESSN`) and normalized perplexity (`Perplexity`). PMC is
considered to have converged when these diagnostics stabilize (Wraith et
al., 2009), or when the normalized perplexity becomes sufficiently close
to 1 (Cappe et al., 2008).

## References

Cappe, O., Douc, R., Guillin, A., Marin, J.M., and Robert, C. (2008).
"Adaptive Importance Sampling in General Mixture Classes". *Statistics
and Computing*, 18, p. 587–600.

Cappe, O., Guillin, A., Marin, J.M., and Robert, C. (2004). "Population
Monte Carlo". *Journal of Computational and Graphical Statistics*, 13,
p. 907–929.

Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian Data
Analysis, Texts in Statistical Science, 2nd ed.". Chapman and Hall,
London.

Wraith, D., Kilbinger, M., Benabed, K., Cappe, O., Cardoso, J.F., Fort,
G., Prunet, S., and Robert, C.P. (2009). "Estimation of Cosmological
Parameters Using Adaptive Importance Sampling". *Physical Review D*,
80(2), p. 023507.

## See also

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`PMC.RAM`](https://robustecologies.github.io/lucifer/reference/PMC.RAM.md),
[`Thin`](https://robustecologies.github.io/lucifer/reference/Thin.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
     beta <- rnorm(Data$J)
     sigma <- runif(1)
     return(c(beta, sigma))
     }
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

Model <- function(parm, Data)
     {
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     LP <- LL + beta.prior + sigma.prior
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

set.seed(666)
Initial.Values <- GIV(Model, MyData, PGF=TRUE)
Fit <- PMC(Model, MyData, Initial.Values, Covar=NULL, Iterations=5,
     Thinning=1, alpha=NULL, M=1, N=100, CPUs=1)
Fit
print(Fit)
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
plot(Fit, BurnIn=0, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData, CPUs=1)
summary(Pred, Discrep="Chi-Square")
plot(Pred, Style="Fitted")
} # }
```
