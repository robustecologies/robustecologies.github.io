# Bayesian mode inference

Estimates the number and location of modes in a distribution using
Sparse Finite Mixture (SFM) models with full posterior uncertainty
quantification. The method fits a mixture model via Gibbs MCMC with a
sparsity-inducing Dirichlet prior that automatically prunes unnecessary
components, then detects modes from the fitted mixture at each MCMC draw
to produce a posterior distribution over both the number and location of
modes.

## Usage

``` r
BayesMode(
  x,
  family = "auto",
  K = 10L,
  iter = 5000L,
  burnin = 2000L,
  thin = 1L,
  alpha.prior = c(1, 200),
  n.threads = 1L,
  verbose = TRUE,
  parm = NULL,
  ...
)
```

## Arguments

- x:

  A numeric vector of observations, or a posterior object of class
  `demonoid`, `laplace`, `iterquad`, `vb`, or `pmc`. When a posterior
  object is provided, mode inference is applied independently to each
  parameter's marginal posterior.

- family:

  Character string specifying the mixture component distribution.
  `"auto"` selects based on data properties (integer non-negative
  \\\to\\ Poisson, positive continuous \\\to\\ Gamma, general continuous
  \\\to\\ Normal). Options: `"normal"`, `"skew_normal"`, `"student_t"`,
  `"gamma"`, `"lognormal"`, `"beta"`, `"poisson"`, `"shifted_poisson"`,
  `"negbin"`, `"geometric"`, `"binomial"`. Default `"auto"`.

- K:

  Maximum number of mixture components. The SFM framework deliberately
  overfits, then prunes via the sparsity prior on \\\alpha\\. Default
  10.

- iter:

  Total number of MCMC iterations. Default 5000.

- burnin:

  Number of burn-in iterations to discard. Default 2000.

- thin:

  Thinning interval. Default 1.

- alpha.prior:

  Numeric vector of length 2 specifying the Gamma hyperprior
  \\(\alpha_a, \alpha_b)\\ on the Dirichlet concentration parameter.
  Default `c(1, 200)`, giving \\E(\alpha) = 1/200\\ which strongly
  favors parsimony.

- n.threads:

  Number of OpenMP threads for the C++ engine. Default 1.

- verbose:

  Logical; print progress. Default `TRUE`.

- parm:

  Character vector of parameter names to analyze when `x` is a posterior
  object. Default `NULL` (all parameters).

- ...:

  Additional arguments passed to the family-specific prior
  specification.

## Value

An object of class `bayes_mode` (when `x` is a numeric vector) or
`bayes_mode_multi` (when `x` is a posterior object). See `print`,
`summary`, and `plot` methods for output.

## Details

The Sparse Finite Mixture approach of Malsiner-Walli et al. (2016) is
extended following Cross et al. (2024) to provide Bayesian inference on
the number and location of modes. The algorithm proceeds in two stages.
First, a mixture of \\K\\ components is fitted via Gibbs MCMC with a
hierarchical Dirichlet prior on the weights: \\\pi \sim
\text{Dirichlet}(\alpha, \ldots, \alpha)\\ where \\\alpha \sim
\text{Gamma}(a\_\alpha, b\_\alpha)\\. Small \\E(\alpha)\\ encourages
sparsity, effectively removing superfluous components. Second, at each
post-burn-in MCMC draw, modes of the fitted mixture density are detected
using family-appropriate algorithms: fixed-point iteration for Normal
mixtures (Carreira-Perpinan, 2000), discrete enumeration for count
distributions, and grid search with golden-section refinement for
general continuous distributions.

The result is a full posterior distribution over the number of modes,
enabling formal probabilistic statements such as \\P(\text{bimodal}) =
0.85\\, and posterior distributions over mode locations with credible
intervals.

## References

Cross, J.L., Hoogerheide, L., Ulker, S. and van Dijk, H.K. (2024).
Bayesian mode inference for discrete distributions in economics and
finance. *Economics Letters*, 235, 111579.
[doi:10.1016/j.econlet.2024.111579](https://doi.org/10.1016/j.econlet.2024.111579)

Malsiner-Walli, G., Fruhwirth-Schnatter, S. and Grun, B. (2016).
Model-based clustering based on sparse finite Gaussian mixtures.
*Statistics and Computing*, 26, 303-324.
[doi:10.1007/s11222-014-9500-2](https://doi.org/10.1007/s11222-014-9500-2)

Carreira-Perpinan, M.A. (2000). Mode-finding for mixtures of Gaussian
distributions. *IEEE Transactions on Pattern Analysis and Machine
Intelligence*, 22(11), 1318-1323.
[doi:10.1109/34.888716](https://doi.org/10.1109/34.888716)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Bimodal continuous data
set.seed(42)
x <- c(rnorm(300, 0, 1), rnorm(300, 6, 1))
bm <- BayesMode(x, family = "normal", K = 10, iter = 5000)
print(bm)
summary(bm)
plot(bm)
plot(bm, type = "n_modes")
plot(bm, type = "mixture")

# Discrete count data
x <- c(rpois(200, 3), rpois(200, 10))
bm <- BayesMode(x, family = "poisson")
print(bm)
plot(bm)

# Posterior mode inference from MCMC output
# result <- lucifer(Model, Data, ...)
# bm <- BayesMode(result)
# print(bm)
} # }
```
