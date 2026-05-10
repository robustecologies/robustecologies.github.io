# The Bayesian bootstrap

This function performs the Bayesian bootstrap of Rubin (1981), returning
either bootstrapped weights or statistics.

## Usage

``` r
BayesianBootstrap(X, n = 1000, Method = "weights", Status = NULL)
```

## Arguments

- X:

  This is a vector or matrix of data. When a matrix is supplied,
  sampling is based on the first column.

- n:

  This is the number of bootstrapped replications.

- Method:

  When `Method="weights"` (which is the default), a matrix of row
  weights is returned. Otherwise, a function is accepted. The function
  specifies the statistic to be bootstrapped. The first argument of the
  function should be a matrix of data, and the second argument should be
  a vector of weights.

- Status:

  This determines the periodicity of status messages. When `Status=100`,
  for example, a status message is displayed every 100 replications.
  Otherwise, `Status` defaults to `NULL`, and status messages are not
  displayed.

## Value

When `Method="weights"`, this function returns a \\N \times n\\ matrix
of weights, where the number of rows \\N\\ is equal to the number of
rows in `X`.

For statistics, a matrix or array is returned, depending on the number
of dimensions. The replicates are indexed by row in a matrix or in the
first dimension of the array.

## Details

The term "bootstrap" comes from the German novel *Adventures of Baron
Munchausen* by Rudolph Raspe, in which the hero saves himself from
drowning by pulling on his own bootstraps. The idea of the statistical
bootstrap is to evaluate properties of an estimator through the
empirical, rather than theoretical, CDF.

Rubin (1981) introduced the Bayesian bootstrap. In contrast to the
frequentist bootstrap which simulates the sampling distribution of a
statistic estimating a parameter, the Bayesian bootstrap simulates the
posterior distribution.

The data, \\\textbf{X}\\, are assumed to be independent and identically
distributed (IID), and to be a representative sample of the larger
(bootstrapped) population. Given that the data has \\N\\ rows in one
bootstrap replication, the row weights are sampled from a Dirichlet
distribution with all \\N\\ concentration parameters equal to \\1\\ (a
uniform distribution over an open standard \\N-1\\ simplex). The
distributions of a parameter inferred from considering many samples of
weights are interpretable as posterior distributions on that parameter.

The Bayesian bootstrap is useful for estimating marginal posterior
covariance and standard deviations for the posterior modes of
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
especially when the model dimension (the number of parameters) is large
enough that estimating the
[`Hessian`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
matrix of second partial derivatives is too computationally demanding.

Just as with the frequentist bootstrap, inappropriate use of the
Bayesian bootstrap can lead to inappropriate inferences. The Bayesian
bootstrap violates the likelihood principle, because the evaluation of a
statistic of interest depends on data sets other than the observed data
set.

The `BayesianBootstrap` function has many uses, including creating test
statistics on the population data given the observed data (supported
here), imputation (with this variation:
[`ABB`](https://robustecologies.github.io/lucifer/reference/ABB.md)),
validation, and more.

## References

Rubin, D.B. (1981). "The Bayesian Bootstrap". *The Annals of
Statistics*, 9(1), p. 130–134.

## See also

[`ABB`](https://robustecologies.github.io/lucifer/reference/ABB.md),
[`Hessian`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)

#Example 1: Samples
x <- 1:2
BB <- BayesianBootstrap(X=x, n=100, Method="weights"); BB

#Example 2: Mean, Univariate
x <- 1:2
BB <- BayesianBootstrap(X=x, n=100, Method=weighted.mean); BB

#Example 3: Mean, Multivariate
data(demonsnacks)
BB <- BayesianBootstrap(X=demonsnacks, n=100,
     Method=function(x,w) apply(x, 2, weighted.mean, w=w)); BB

#Example 4: Correlation
dye <- c(1.15, 1.70, 1.42, 1.38, 2.80, 4.70, 4.80, 1.41, 3.90)
efp <- c(1.38, 1.72, 1.59, 1.47, 1.66, 3.45, 3.87, 1.31, 3.75)
X <- matrix(c(dye,efp), length(dye), 2)
colnames(X) <- c("dye","efp")
BB <- BayesianBootstrap(X=X, n=100,
     Method=function(x,w) cov.wt(x, w, cor=TRUE)$cor); BB
} # }
```
