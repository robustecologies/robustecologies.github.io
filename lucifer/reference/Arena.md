# Universal cross-method benchmarking

Compares heterogeneous Bayesian inference methods on the same problem,
computing unified efficiency, accuracy, and reliability metrics across
MCMC, variational Bayes, Laplace approximation, iterative quadrature,
PMC, SMC, ABC, and SBI fit objects. Returns an S3 object of class
`arena` with print, summary, and plot methods.

## Usage

``` r
Arena(
  x = NULL,
  Model = NULL,
  Data = NULL,
  Initial.Values = NULL,
  methods = NULL,
  reference = NULL,
  n.samples = 1000,
  verbose = TRUE
)
```

## Arguments

- x:

  A named list of fit objects (Mode A). Each element must be an object
  of class `demonoid`, `vb`, `laplace`, `iterquad`, `pmc`, `smc`, `abc`,
  or `sbi`. Names are used as method labels in output.

- Model:

  Model specification function for Mode B auto-run. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  Data list for Mode B. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Initial.Values:

  Initial parameter vector for Mode B.

- methods:

  For Mode B, a named list of method specifications. Each element should
  be a list with at minimum `type` (one of "MCMC", "VB", "Laplace",
  "IQ", "PMC", "SMC", "ABC", "SBI") and the relevant fitting arguments.

- reference:

  A reference for accuracy comparison: a fit object, a list with `mean`
  and `sd` vectors (analytical truth), or `NULL` for consensus reference
  computed from all methods.

- n.samples:

  Integer. Posterior samples to retain per method. Default 1000.

- verbose:

  Logical. If `TRUE`, print progress. Default `TRUE`.

## Value

An object of class `arena` containing:

- Fits:

  The original fit objects.

- Metrics:

  Data frame with per-method metrics including ESS, ESS/second,
  convergence, DIC, LML, KLD and Wasserstein to reference, and composite
  ranks.

- Posteriors:

  Per-method posterior summaries and samples.

- Pairwise:

  Pairwise KLD, Wasserstein, and mean-difference matrices between all
  methods.

- Pareto:

  Data frame identifying the Pareto frontier in the time-accuracy plane.

- Reference:

  The reference distribution used (type, mean, sd, samples).

- Parameters:

  Parameter names.

- Call:

  The matched call.

- Minutes:

  Total wall-clock time in minutes.

## Details

Arena operates in two modes. In Mode A (pre-computed), the user passes a
named list of existing fit objects via `x`. Each object is processed
through a class-specific extraction adapter that normalizes posterior
samples, ESS, convergence status, wall-clock time, and information
criteria into a standardized format, enabling apples-to- apples
comparison across fundamentally different inference paradigms.

In Mode B (auto-run), the user supplies `Model`, `Data`,
`Initial.Values`, and a named list of method specifications via
`methods`. Arena runs each method sequentially, catches failures
gracefully (reporting which methods failed), and then proceeds as in
Mode A on the successful fits. This mode is useful for systematic
benchmarking where all methods should start from identical conditions.

The metrics framework evaluates three orthogonal dimensions of inference
quality:

**Efficiency.** Measured as ESS per second of wall-clock time. For
stochastic methods, the minimum ESS across all parameters is used,
reflecting the bottleneck dimension. For deterministic methods (Laplace,
IQ), ESS is estimated from the SIR resampling step or from the effective
number of independent quadrature points.

**Accuracy.** Quantified by the marginal Kullback-Leibler divergence and
Wasserstein-1 distance between each method's posterior and a reference
distribution. When an analytical reference with known mean and standard
deviation is supplied, accuracy reflects absolute fidelity. When no
reference is provided, Arena constructs a consensus reference by
ESS-weighted averaging of posterior moments across all methods,
incorporating between-method variance: \$\$\mu\_{\mathrm{ref}} = \sum_m
w_m \hat{\mu}\_m, \quad \sigma^2\_{\mathrm{ref}} = \sum_m w_m
(\hat{\sigma}\_m^2 + (\hat{\mu}\_m - \mu\_{\mathrm{ref}})^2)\$\$ where
\\w_m = \mathrm{ESS}\_m / \sum_j \mathrm{ESS}\_j\\.

**Reliability.** Assessed via convergence diagnostics: \\\hat{R}\\ for
MCMC, ELBO convergence for VB, optimization convergence for Laplace/IQ,
and particle ESS for SMC. A method is marked as converged only if all
diagnostics pass their respective thresholds.

Pairwise KLD and Wasserstein matrices between all method posteriors are
also computed, enabling cluster analysis of which methods agree. The
Pareto frontier identifies methods that are not dominated in the
time-accuracy tradeoff, and a composite rank aggregates efficiency,
accuracy, and reliability into an overall recommendation.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). "Rank-normalization, folding, and localization: an improved Rhat
for assessing convergence of MCMC." *Bayesian Analysis*, 16(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

Rice, J.R. (1976). "The algorithm selection problem." *Advances in
Computers*, 15, 65–118.
[doi:10.1016/S0065-2458(08)60520-3](https://doi.org/10.1016/S0065-2458%2808%2960520-3)

Wolpert, D.H. and Macready, W.G. (1997). "No free lunch theorems for
optimization." *IEEE Transactions on Evolutionary Computation*, 1(1),
67–82. [doi:10.1109/4235.585893](https://doi.org/10.1109/4235.585893)

## See also

[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`RobustBayes`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
[`plot.arena`](https://robustecologies.github.io/lucifer/reference/plot.arena.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data for a simple regression
set.seed(42)
N <- 30
x <- rnorm(N)
y <- 0.5 + 1.2 * x + rnorm(N, 0, 0.8)

Data <- list(
  mon.names = "LP",
  parm.names = c("beta0", "beta1", "log.sigma"),
  y = y, x = x, N = N
)

Model <- function(parm, Data) {
  beta0 <- parm[1]; beta1 <- parm[2]; sigma <- exp(parm[3])
  mu <- beta0 + beta1 * Data$x
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + sum(dnorm(parm[1:2], 0, 10, log = TRUE)) +
             dhalfcauchy(sigma, 5, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}

IV <- c(0, 0, 0)

# Mode A: compare pre-computed fits
fit_nuts <- lucifer(Model, Data, IV, "NUTS", 5000, Status = 1000)
fit_amwg <- lucifer(Model, Data, IV, "AMWG", 5000, Status = 1000)
fit_vb   <- VariationalBayes(Model, Data, IV, Method = "Pathfinder")
fit_la   <- LaplaceApproximation(Model, Data, IV, sir = TRUE)

ar <- Arena(list(NUTS = fit_nuts, AMWG = fit_amwg,
                 VB = fit_vb, Laplace = fit_la))
print(ar)
summary(ar)
plot(ar)
plot(ar, type = "accuracy")
plot(ar, type = "pareto")
plot(ar, type = "heatmap")

# Mode B: auto-run with method specifications
ar2 <- Arena(Model = Model, Data = Data, Initial.Values = IV,
             methods = list(
               NUTS = list(type = "MCMC", Algorithm = "NUTS",
                           Iterations = 5000),
               AMWG = list(type = "MCMC", Algorithm = "AMWG",
                           Iterations = 5000),
               VB   = list(type = "VB", Method = "Pathfinder")))
plot(ar2, type = "forest")
} # }
```
