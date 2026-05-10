# Pre-fit algorithm advisor

Profiles a Bayesian model to characterize its computational properties
and recommends the most suitable inference strategy from the full suite
of algorithms available in lucifer: MCMC, variational Bayes, Laplace
approximation, iterative quadrature, PMC, SMC, ABC, and SBI.

## Usage

``` r
Prescribe(
  Model,
  Data,
  Initial.Values,
  n.profile = 200,
  methods = c("MCMC", "VB", "Laplace", "IQ", "PMC", "SMC", "ABC", "SBI"),
  verbose = TRUE
)
```

## Arguments

- Model:

  A model specification function compatible with
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  A named list of data. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  for the required structure.

- Initial.Values:

  A numeric vector of initial parameter values.

- n.profile:

  Integer. Number of random perturbation evaluations used to assess
  posterior surface roughness and multimodality. Default 200.

- methods:

  Character vector specifying which inference categories to consider.
  Default includes all available categories.

- verbose:

  Logical. If `TRUE`, prints progress messages during profiling. Default
  `TRUE`.

## Value

An object of class `prescription` containing:

- Profile:

  Named list of model characteristics (dimensionality, evaluation speed,
  gradient availability, constraints, conditioning, multimodality score,
  likelihood-free status, discrete parameters, sample size).

- Scores:

  Data frame with columns Method, Category, Subcategory, Score, and Rank
  for all evaluated methods.

- Recommend:

  List with `primary` (top recommendation including method name,
  category, score, justification, and ready-to-paste R code) and
  `alternatives` (top 5 alternatives).

- Warnings:

  Character vector of diagnostic warnings.

- Call:

  The matched call.

- Minutes:

  Wall-clock time for the profiling phase.

## Details

The function proceeds in two phases. The profiling phase evaluates the
model approximately `n.profile + 2d` times (where \\d\\ is the number of
parameters) to characterize dimensionality, evaluation speed, gradient
availability, constraint structure, posterior conditioning, and
multimodality. No fitting is performed, so the computational cost is
typically negligible compared to actual inference.

The scoring phase assigns each candidate method a composite score via
multiplicative factors. A zero factor disqualifies a method entirely.
Scores reflect the expected quality-adjusted efficiency of each method
given the problem's profile. The method with the highest score is
recommended as the primary choice, with alternatives ranked by
decreasing score.

## References

Betancourt, M. (2017). A Conceptual Introduction to Hamiltonian Monte
Carlo. *arXiv preprint arXiv:1701.02434*.

Blei, D. M., Kucukelbir, A. and McAuliffe, J. D. (2017). Variational
Inference: A Review for Statisticians. *Journal of the American
Statistical Association*, 112(518), 859-877.
[doi:10.1080/01621459.2017.1285773](https://doi.org/10.1080/01621459.2017.1285773)

## See also

[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md),
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md),
[`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define a simple normal-normal model
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rep(mu, Data$N)
  return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
    yhat = yhat, parm = parm))
}
Data <- list(N = 50, y = rnorm(50, 5, 2),
  mon.names = "LP", parm.names = c("mu", "log.sigma"))
IV <- c(0, 0)

# Get recommendation
rx <- Prescribe(Model, Data, IV)
print(rx)
summary(rx)
plot(rx)
} # }
```
