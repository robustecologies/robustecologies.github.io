# Maximum likelihood estimation using Bayesian MCMC via data cloning

Data cloning (Lele, Dennis & Lutscher 2007) uses Bayesian MCMC to
compute frequentist maximum likelihood estimates and standard errors for
complex hierarchical models. The observed data likelihood is raised to
the Kth power by multiplying the log-likelihood by K, causing the
posterior to concentrate on the MLE as K grows. The posterior mean
converges to the MLE and K times the posterior variance converges to the
inverse Fisher information matrix.

## Usage

``` r
DataCloning(
  Model,
  Data,
  Initial.Values,
  Covar = NULL,
  Iterations = 10000,
  Status = 100,
  Thinning = 10,
  Algorithm = "NUTS",
  Specs = NULL,
  Debug = list(DB.chol = FALSE, DB.eigen = FALSE, DB.MCSE = FALSE, DB.Model = TRUE),
  Chains = 3L,
  CPUs = 3L,
  Type = NULL,
  Grad = NULL,
  n_clones = c(1, 2, 4, 8, 16),
  n_priors = 1L,
  prior_generator = NULL,
  update_prior = FALSE,
  conf_level = 0.95,
  store_posteriors = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- Model:

  A function taking arguments `parm` and `Data`, returning a list with
  components `LP`, `Dev`, `Monitor`, `yhat`, and `parm`. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  for details.

- Data:

  A named list containing at minimum `mon.names` and `parm.names`
  character vectors.

- Initial.Values:

  A numeric vector of starting parameter values.

- Covar:

  Proposal covariance matrix or vector. Passed to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Iterations:

  Integer. MCMC iterations per clone level. Default 10000.

- Status:

  Integer. Printing frequency. Default 100.

- Thinning:

  Integer. Thinning interval. Default 10.

- Algorithm:

  Character. MCMC algorithm name or abbreviation. Default `"NUTS"`.

- Specs:

  List of algorithm specifications. Default `NULL`.

- Debug:

  List of debug flags. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Chains:

  Integer. Number of parallel chains. Default 3.

- CPUs:

  Integer. Number of CPUs. Default 3.

- Type:

  Character. Parallelization type (`"FORK"` or `"PSOCK"`). Default
  `NULL` (auto-detect).

- Grad:

  Optional gradient function. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- n_clones:

  Integer vector of clone sizes, sorted ascending. Default
  `c(1, 2, 4, 8, 16)`.

- n_priors:

  Integer \>= 1. Number of distinct priors for estimability testing.
  Default 1.

- prior_generator:

  Function `function(prior_index, Model, Data)` returning a modified
  Model function with a different prior. When `NULL` and `n_priors > 1`,
  automatic location-shifted priors are generated from the K=1 posterior
  of the original model.

- update_prior:

  Logical. If `TRUE`, use posterior mean and covariance from clone K as
  initial values and proposal covariance for clone K+1 (warm starting).
  Default `FALSE`.

- conf_level:

  Numeric. Confidence level for Wald intervals. Default 0.95.

- store_posteriors:

  Logical. If `TRUE`, retain full posterior matrices for each K. Default
  `TRUE`.

- verbose:

  Logical. Console progress output. Default `TRUE`.

- ...:

  Additional arguments passed to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Value

An object of class `data_cloning` containing:

- MLE:

  Numeric vector of maximum likelihood estimates.

- SE:

  Numeric vector of asymptotic standard errors.

- Asymptotic.Cov:

  Matrix. K_max times the posterior covariance from the largest clone
  run, estimating the inverse Fisher information.

- Fisher.Info:

  Matrix. Estimated Fisher information (inverse of `Asymptotic.Cov`).

- Wald.CI:

  Matrix with columns `lower` and `upper`.

- per_clone:

  List of per-K summary lists.

- eigenvalue_diagnostic:

  Eigenvalue convergence diagnostics.

- convergence:

  Scaled variance convergence assessment.

- mse_r2:

  MSE and R-squared per K.

- estimability:

  ANOVA estimability results (NULL if n_priors == 1).

- final_fit:

  The `demonoid` object from the largest K run.

## Details

**Theoretical foundation.** Given observed data \\y\\ and a parametric
model with likelihood \\L(\theta \mid y)\\, define the K-cloned
posterior as \$\$\pi_K(\theta \mid y) \propto \[L(\theta \mid y)\]^K \\
\pi(\theta)\$\$ where \\\pi(\theta)\\ is the prior. Walker (1969)
established the asymptotic behaviour of posterior distributions under
repeated observations; data cloning exploits the same mechanism by
artificially amplifying the information content of a fixed dataset. As
\\K \to \infty\\, three key properties emerge: (i) the posterior mean
converges to the MLE, \\\bar{\theta}\_K \to
\hat{\theta}\_{\mathrm{MLE}}\\; (ii) the scaled posterior covariance
converges to the inverse Fisher information, \\K \\
\mathrm{Var}\_K(\theta) \to I(\hat{\theta})^{-1}\\; and (iii) the
posterior becomes asymptotically normal, \\\theta \mid y^{(K)}
\xrightarrow{d} N(\hat{\theta}\_{\mathrm{MLE}},\\ K^{-1}
I(\hat{\theta})^{-1})\\. This connection to Robert's (1993) prior
feedback framework means that data cloning can be understood as a
deterministic annealing scheme in which the likelihood temperature is
systematically increased.

**Implementation via LP/Dev decomposition.** The lucifer Model interface
guarantees that the log-posterior decomposes as \\\mathrm{LP} =
\ell(\theta) + \log \pi(\theta)\\ and the deviance satisfies
\\\mathrm{Dev} = -2\\\ell(\theta)\\, where \\\ell(\theta)\\ is the
log-likelihood. The cloned model therefore computes \\\mathrm{LP}\_K =
K\\\ell + \log \pi\\ and \\\mathrm{Dev}\_K = K \cdot \mathrm{Dev}\\
without any actual data duplication. A thin wrapper around the user's
Model function implements this transformation; all 55+ MCMC algorithms,
multi-chain parallel infrastructure, and modern diagnostics (split
R-hat, bulk/tail ESS, MCSE) work transparently through the wrapper.

**Convergence diagnostics.** Three complementary diagnostics are
computed automatically. The eigenvalue diagnostic (Lele et al. 2010)
computes the standardized largest eigenvalue ratio \\\lambda_S(K) =
\lambda\_{\max}(\Sigma_K) / \lambda\_{\max}(\Sigma_1)\\, which should
decrease at rate \\1/K\\ for identifiable models. The scaled variance
check monitors \\K \cdot \mathrm{Var}\_K(\theta_j)\\ for each parameter;
stabilization (constant value across K) indicates that the asymptotic
variance has been reached. The MSE/R-squared diagnostic (Lele et al.
2010) computes the mean Mahalanobis distance of posterior draws from the
posterior mean; under convergence this should approximate \\\chi^2(p)\\
where \\p\\ is the number of parameters.

**Estimability and identifiability.** When `n_priors > 1`, the Campbell
and Lele (2014) two-stage ANOVA test is applied to assess parameter
estimability. A parameter is *estimable* if its posterior mean changes
with K (concentrates on MLE) but does not depend on the choice of prior.
Non-estimable parameters indicate structural non-identifiability in the
model. The test proceeds by fitting the data cloning sequence under
multiple distinct priors (supplied via `prior_generator` or generated
automatically); a two-way ANOVA then tests for the effects of clone size
and prior on each posterior mean. When no `prior_generator` is supplied,
the function runs K=1 with the original model first, then creates
informative alternative priors centered at \\\bar{\theta}\_1 \pm 3 \\
\mathrm{SD}\_1\\ with standard deviation equal to the K=1 posterior SD.
This location-shift approach ensures that the alternative priors
genuinely pull the posterior away from the MLE at low K, creating a
detectable prior effect that disappears at high K for estimable
parameters. The ANOVA classification is cross-validated against the
eigenvalue and scaled-variance diagnostics to guard against false
negatives when the original prior is already sufficiently diffuse.

**Warm starting.** Setting `update_prior = TRUE` passes the posterior
mean and covariance from clone K as the initial values and proposal
covariance for clone K+1. This substantially improves sampling
efficiency at high K, where the posterior becomes very narrow and a
sampler starting from the prior region would waste most iterations in
adaptation.

**Numerical considerations.** At large K the cloned log-likelihood \\K
\cdot \ell(\theta)\\ can overflow to `-Inf`; the function validates that
the cloned LP is finite at the initial values and skips any K with a
warning if not. The Fisher information matrix is estimated as
\\I(\hat{\theta}) = \[K \cdot \Sigma_K\]^{-1}\\; when the condition
number of the asymptotic covariance exceeds
\\1/\sqrt{\epsilon\_{\mathrm{mach}}}\\, a pseudoinverse via SVD is used
instead. Note that the DIC from cloned runs is meaningless since
\\\mathrm{Dev}\_K = K \cdot \mathrm{Dev}\\; it should not be interpreted
as a model comparison criterion.

## References

Lele, S.R., Dennis, B. & Lutscher, F. (2007). Data cloning: easy maximum
likelihood estimation for complex ecological models using Bayesian
Markov chain Monte Carlo methods. *Ecology Letters*, 10, 551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

Lele, S.R., Nadeem, K. & Schmuland, B. (2010). Estimability and
likelihood inference for generalized linear mixed models using data
cloning. *J. Amer. Statist. Assoc.*, 105(492), 1617-1625.
[doi:10.1198/jasa.2010.tm09757](https://doi.org/10.1198/jasa.2010.tm09757)

Campbell, D. & Lele, S. (2014). An ANOVA test for parameter estimability
using data cloning with application to statistical inference for dynamic
systems. *Comput. Statist. Data Anal.*, 70, 257-267.
[doi:10.1016/j.csda.2013.09.013](https://doi.org/10.1016/j.csda.2013.09.013)

Ponciano, J.M., Taper, M.L., Dennis, B. & Lele, S.R. (2009).
Hierarchical models in ecology: confidence intervals, hypothesis
testing, and model selection using data cloning. *Ecology*, 90, 356-362.
[doi:10.1890/07-1960.1](https://doi.org/10.1890/07-1960.1)

Walker, A.M. (1969). On the asymptotic behaviour of posterior
distributions. *J. Roy. Statist. Soc. B*, 31, 80-88.
[doi:10.1111/j.2517-6161.1969.tb00767.x](https://doi.org/10.1111/j.2517-6161.1969.tb00767.x)

Robert, C.P. (1993). Prior feedback: a Bayesian approach to maximum
likelihood estimation. *Comput. Statist.*, 8, 279-294.

Ponciano, J.M., Burleigh, J.G., Braun, E.L. & Taper, M.L. (2012).
Assessing parameter identifiability in phylogenetic models using data
cloning. *Syst. Biol.*, 61(6), 955-972.
[doi:10.1093/sysbio/sys055](https://doi.org/10.1093/sysbio/sys055)

## See also

[`print.data_cloning`](https://robustecologies.github.io/lucifer/reference/print.data_cloning.md)
for concise output,
[`summary.data_cloning`](https://robustecologies.github.io/lucifer/reference/summary.data_cloning.md)
for extended diagnostics,
[`plot.data_cloning`](https://robustecologies.github.io/lucifer/reference/plot.data_cloning.md)
for visualization (10 plot types). Related inference functions:
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`SMC`](https://robustecologies.github.io/lucifer/reference/SMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Normal model: MLE should converge to sample mean
set.seed(666)
y <- rnorm(100, mean = 3, sd = 1)
Data <- list(
    mon.names = "LP",
    parm.names = "mu",
    y = y, n = length(y)
)
Model <- function(parm, Data) {
    mu <- parm[1]
    LL <- sum(dnorm(Data$y, mu, 1, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rep(mu, Data$n), parm = parm)
}

# Run data cloning
dc <- DataCloning(Model, Data, Initial.Values = 0,
                  n_clones = c(1, 2, 4, 8, 20),
                  Iterations = 5000, Algorithm = "NUTS",
                  n_priors = 3)

# Inspect results
print(dc)
summary(dc)

# Visualize diagnostics (includes convergence, eigenvalue,
# scaled variance, and density in a 3x2 dashboard)
plot(dc, type = "diagnostics")

# Individual plot types
plot(dc, type = "convergence")
plot(dc, type = "eigenvalue")

# Plot likelihood profile
plot(dc, Model = Model, Data = Data, type = "profile", n_grid = 100)
} # }
```
