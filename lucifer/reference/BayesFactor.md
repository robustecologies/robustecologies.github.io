# Bayes factor

This function calculates Bayes factors for two or more fitted objects of
class `demonoid`, `iterquad`, `laplace`, `pmc`, or `vb` that were
estimated respectively with the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
functions, and indicates the strength of evidence in favor of the
hypothesis (that each model, \\\mathcal{M}\_i\\, is better than another
model, \\\mathcal{M}\_j\\).

## Usage

``` r
BayesFactor(x)
```

## Arguments

- x:

  This is a list of two or more fitted objects of class `demonoid`,
  `iterquad`, `laplace`, `pmc`, or `vb`. The components are named in
  order beginning with model 1, `M1`, and \\k\\ models are usually
  represented as \\\mathcal{M}\_1,\dots,\mathcal{M}\_k\\.

## Value

`BayesFactor` returns an object of class `bayesfactor` that is a list
with the following components:

- B:

  This is a matrix of Bayes factors.

- Hypothesis:

  This is the hypothesis, and is stated as "row \> column", indicating
  that the model associated with the row of an element in matrix `B` is
  greater than the model associated with the column of that element.

- Strength.of.Evidence:

  This is the strength of evidence in favor of the hypothesis.

- Posterior.Probability:

  This is a vector of the posterior probability of each model, given
  flat priors.

## Details

Introduced by Harold Jeffreys, a "Bayes factor" is a Bayesian
alternative to frequentist hypothesis testing that is most often used
for the comparison of multiple models by hypothesis testing, usually to
determine which model better fits the data (Jeffreys, 1961). Bayes
factors are notoriously difficult to compute, and the Bayes factor is
only defined when the marginal density of \\\textbf{y}\\ under each
model is proper (see
[`is.proper`](https://robustecologies.github.io/lucifer/reference/is.proper.md)).
However, the Bayes factor is easy to approximate with the
Laplace-Metropolis estimator (Lewis and Raftery, 1997) and other methods
of approximating the logarithm of the marginal likelihood (for more
information, see
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md)).

Hypothesis testing with Bayes factors is more robust than frequentist
hypothesis testing, since the Bayesian form avoids model selection bias,
evaluates evidence in favor of the null hypothesis, includes model
uncertainty, and allows non-nested models to be compared (though of
course the model must have the same dependent variable). Also,
frequentist significance tests become biased in favor of rejecting the
null hypothesis with sufficiently large sample size.

The Bayes factor for comparing two models may be approximated as the
ratio of the marginal likelihood of the data in model 1 and model 2.
Formally, the Bayes factor in this case is

\$\$B =
\frac{p(\textbf{y}\|\mathcal{M}\_1)}{p(\textbf{y}\|\mathcal{M}\_2)} =
\frac{\int
p(\textbf{y}\|\Theta_1,\mathcal{M}\_1)p(\Theta_1\|\mathcal{M}\_1)d\Theta_1}{\int
p(\textbf{y}\|\Theta_2,\mathcal{M}\_2)p(\Theta_2\|\mathcal{M}\_2)d\Theta_2}\$\$

where \\p(\textbf{y}\|\mathcal{M}\_1)\\ is the marginal likelihood of
the data in model 1.

The
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
functions each return the
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md), the
approximate logarithm of the marginal likelihood of the data, in each
fitted object of class `iterquad`, `laplace`, `demonoid`, `pmc`, or
`vb`. The `BayesFactor` function calculates matrix `B`, a matrix of
Bayes factors, where each element of matrix `B` is a comparison of two
models. Each Bayes factor is calculated as the exponentiated difference
of `LML` of model 1 (\\\mathcal{M}\_1\\) and `LML` of model 2
(\\\mathcal{M}\_2\\), and the hypothesis for each element of matrix `B`
is that the model associated with the row is greater than the model
associated with the column. For example, element `B[3,2]` is the Bayes
factor that model 3 is greater than model 2. The "Strength of Evidence"
aids in the interpretation (Jeffreys, 1961).

Each Bayes factor, `B`, is the posterior odds in favor of the hypothesis
divided by the prior odds in favor of the hypothesis, where the
hypothesis is usually \\\mathcal{M}\_1 \> \mathcal{M}\_2\\. For example,
when `B[3,2]=2`, the data favor \\\mathcal{M}\_3\\ over
\\\mathcal{M}\_2\\ with 2:1 odds.

It is also popular to consider the natural logarithm of the Bayes
factor. The scale of the logged Bayes factor is the same above and below
one, which is more appropriate for visual comparisons. For example, when
comparing two Bayes factors at 0.5 and 2, the logarithm of these Bayes
factors is -0.69 and 0.69.

Gelman finds Bayes factors generally to be irrelevant, because they
compute the relative probabilities of the models conditional on one of
them being true. Gelman prefers approaches that measure the distance of
the data to each of the approximate models (Gelman et al., 2004, p.
180), such as with posterior predictive checks (see the
[`predict.iterquad`](https://robustecologies.github.io/lucifer/reference/predict.iterquad.md)
function regarding iterative quadrature,
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md)
function in the context of Laplace Approximation,
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md)
function in the context of MCMC,
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md)
function in the context of PMC, or
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md)
function in the context of Variational Bayes). Kass et al. (1995)
asserts this can be done without assuming one model is the true model.

## References

Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian Data
Analysis, Texts in Statistical Science, 2nd ed.". Chapman and Hall,
London.

Jeffreys, H. (1961). "Theory of Probability, Third Edition". Oxford
University Press: Oxford, England.

Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". *Journal of the
American Statistical Association*, 90(430), p. 773–795.

Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
Posterior Simulation with the Laplace-Metropolis Estimator". *Journal of
the American Statistical Association*, 92, p. 648–655.

## See also

[`is.bayesfactor`](https://robustecologies.github.io/lucifer/reference/is.class.md),
[`is.proper`](https://robustecologies.github.io/lucifer/reference/is.proper.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`predict.iterquad`](https://robustecologies.github.io/lucifer/reference/predict.iterquad.md),
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# The following example fits a model as Fit1, then adds a predictor, and
# fits another model, Fit2. The two models are compared with Bayes
# factors.

library(lucifer)

##############################  Demon Data  ###############################
data(demonsnacks)
J <- 2
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10]+1)))
X[,2] <- CenterScale(X[,2])

#########################  Data List Preparation  #########################
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

##########################  Model Specification  ##########################
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

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

########################  Laplace Approximation  ##########################
Fit1 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=10000)
Fit1

##############################  Demon Data  ###############################
data(demonsnacks)
J <- 3
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(7,8)]))
X[,2] <- CenterScale(X[,2])
X[,3] <- CenterScale(X[,3])

#########################  Data List Preparation  #########################
mon.names <- c("sigma","mu[1]")
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) return(c(rnormv(Data$J,0,10), rhalfcauchy(1,5)))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

########################  Laplace Approximation  ##########################
Fit2 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
     Iterations=10000)
Fit2

#############################  Bayes Factor  ##############################
Model.list <- list(M1=Fit1, M2=Fit2)
BayesFactor(Model.list)
} # }
```
