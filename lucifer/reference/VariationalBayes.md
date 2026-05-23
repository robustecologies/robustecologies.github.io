# Variational Bayes

The `VariationalBayes` function is a numerical approximation method for
deterministically estimating the marginal posterior distributions
(target distributions) in a Bayesian model with approximated
distributions by minimizing the Kullback-Leibler Divergence
([`KLD`](https://robustecologies.github.io/lucifer/reference/KLD.md))
between the target and its approximation.

## Usage

``` r
VariationalBayes(
  Model,
  parm,
  Data,
  Covar = NULL,
  Interval = 1e-06,
  Iterations = 1000,
  Method = "Salimans2",
  Samples = 1000,
  sir = TRUE,
  Stop.Tolerance = 1e-05,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  The model from a user-defined function. The user-defined function is
  where the model is specified. `VariationalBayes` passes two arguments
  to the model function, `parms` and `Data`. For more information, see
  the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function and "lucifer tutorial" vignette.

- parm:

  A vector of initial values equal in length to the number of
  parameters. `VariationalBayes` will attempt to optimize these initial
  values for the parameters, where the optimized values are the
  posterior means, for later use with the
  [`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  or [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
  function. The
  [`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md)
  function may be used to randomly generate initial values. Parameters
  must be continuous.

- Data:

  A list of data. The list of data must include `mon.names` which
  contains monitored variable names, and `parm.names` which contains
  parameter names. `VariationalBayes` must be able to determine the
  sample size of the data, and will look for a scalar sample size
  variable `n` or `N`. If not found, it will look for variable `y` or
  `Y`, and attempt to take its number of rows as sample size.
  `VariationalBayes` needs to determine sample size due to the
  asymptotic nature of this method. Sample size should be at least
  \\\sqrt{J}\\ with \\J\\ exchangeable parameters.

- Covar:

  Defaults to `NULL`, but may otherwise accept a \\K \times K\\
  covariance matrix (where \\K\\ is the number of dimensions or
  parameters) of the parameters. When the model is updated for the first
  time and prior variance or covariance is unknown, then `Covar=NULL`
  should be used. Once `VariationalBayes` has finished updating, it may
  be desired to continue updating where it left off, in which case the
  covariance matrix from the last run can be input into the next run.
  For `Method="SVGD"`, `Covar` may be an integer specifying the number
  of particles (default 20). For `Method="Pathfinder"`, `Covar` may be
  an integer specifying the number of independent paths (default 4).

- Interval:

  An interval for estimating approximate gradients. The logarithm of the
  unnormalized joint posterior density of the Bayesian model is
  evaluated at the current parameter value, and again at the current
  parameter value plus this interval.

- Iterations:

  An integer that determines the number of iterations that
  `VariationalBayes` will attempt to maximize the logarithm of the
  unnormalized joint posterior density. Defaults to 1000.
  `VariationalBayes` will stop before this number of iterations if the
  tolerance is less than or equal to the `Stop.Tolerance` criterion. The
  required amount of computer memory increases with `Iterations`. If
  computer memory is exceeded, then all will be lost.

- Method:

  A character string specifying the variational inference algorithm.
  Accepts `"Salimans1"` (gradient-only fixed-form VB), `"Salimans2"`
  (gradient + Hessian fixed-form VB, the default), `"ADVI.mf"`
  (automatic differentiation VI, mean-field), `"ADVI.fr"` (automatic
  differentiation VI, full-rank), `"BBVI"` (black-box VI with score
  function estimator), `"SVGD"` (Stein variational gradient descent),
  `"NGD"` (natural gradient descent), or `"Pathfinder"` (quasi-Newton VI
  along L-BFGS trajectory). See Details.

- Samples:

  The number of posterior samples to be taken with sampling importance
  resampling via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function, which occurs only when `sir=TRUE`. Note that the number of
  samples should increase with the number and intercorrelations of the
  parameters.

- sir:

  Logical indicating whether or not Sampling Importance Resampling (SIR)
  is conducted via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function to draw independent posterior samples. Defaults to `TRUE`.
  Even when `TRUE`, posterior samples are drawn only when
  `VariationalBayes` has converged. Posterior samples are required for
  many other functions, including `plot.vb` and `predict.vb`. The only
  time that it is advantageous for `sir=FALSE` is when
  `VariationalBayes` is used to help the initial values for
  [`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  or
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
  and it is unnecessary for time to be spent on sampling. Less time can
  be spent on sampling by increasing `CPUs`, which parallelizes the
  sampling.

- Stop.Tolerance:

  Any positive number, defaults to 1.0E-5. Tolerance is calculated each
  iteration, and the criteria varies by algorithm. The algorithm is
  considered to have converged to the user-specified `Stop.Tolerance`
  when the tolerance is less than or equal to the value of
  `Stop.Tolerance`, and the algorithm terminates at the end of the
  current iteration. Often, multiple criteria are used, in which case
  the maximum of all criteria becomes the tolerance.

- CPUs:

  An integer that specifies the number of central processing units
  (CPUs) of the multicore computer or computer cluster. Defaults to
  `CPUs=1`, in which parallel processing does not occur. Parallelization
  occurs only for sampling with
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  when `sir=TRUE`.

- Type:

  The type of parallel processing to perform, accepting either
  `Type="PSOCK"` or `Type="MPI"`.

## Value

An object of class `vb` that is a list with the following components:

- Call:

  The matched call of `VariationalBayes`.

- Converged:

  Logical indicator of whether or not `VariationalBayes` converged
  within the specified `Iterations` according to the supplied
  `Stop.Tolerance` criterion.

- Covar:

  The estimated covariance matrix. The `Covar` matrix may be scaled and
  input into the `Covar` argument of the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  or [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
  function for further estimation, or the diagonal of this matrix may be
  used to represent the posterior variance of the parameters, provided
  the algorithm converged and matrix inversion was successful. To scale
  this matrix for use with Lucifer or PMC, multiply it by \\2.38^2/d\\,
  where \\d\\ is the number of initial values.

- Deviance:

  A vector of the iterative history of the deviance in the
  `VariationalBayes` function, as it sought convergence.

- History:

  An array of the iterative history of the parameters in the
  `VariationalBayes` function, as it sought convergence. The first
  matrix is for means and the second matrix is for variances.

- Initial.Values:

  The vector of initial values that was originally given to
  `VariationalBayes` in the `parm` argument.

- LML:

  An approximation of the logarithm of the marginal likelihood of the
  data (see the
  [`LML`](https://robustecologies.github.io/lucifer/reference/LML.md)
  function for more information).

- LP.Final:

  The final scalar value for the logarithm of the unnormalized joint
  posterior density.

- LP.Initial:

  The initial scalar value for the logarithm of the unnormalized joint
  posterior density.

- Method:

  The variational inference algorithm used.

- Minutes:

  The number of minutes that `VariationalBayes` was running.

- Monitor:

  When `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix.

- Posterior:

  When `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix.

- Step.Size.Final:

  The final, scalar step size value at the end of the `VariationalBayes`
  algorithm.

- Step.Size.Initial:

  The initial, scalar step size.

- Summary1:

  A summary matrix of point-estimated posterior means and variances with
  columns Mean, SD, LB, and UB (95% probability interval).

- Summary2:

  A summary matrix of posterior samples drawn with sampling importance
  resampling
  ([`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md))
  when `sir=TRUE`, with columns Mean, SD, LB, and UB (95% probability
  interval).

- Tolerance.Final:

  The last tolerance of the `VariationalBayes` algorithm.

- Tolerance.Stop:

  The `Stop.Tolerance` criterion.

## Details

Variational Bayes (VB) is a family of numerical approximation algorithms
for deterministic posterior inference. The core idea is to posit a
tractable family of distributions \\q(\theta; \lambda)\\ parameterized
by variational parameters \\\lambda\\, then minimize the
Kullback-Leibler divergence \\\mathrm{KL}(q \\ p)\\ between the
approximation and the true posterior. Equivalently, VB maximizes the
evidence lower bound (ELBO): \$\$\mathcal{L}(\lambda) =
\mathbb{E}\_q\[\log p(\theta, y)\] - \mathbb{E}\_q\[\log q(\theta;
\lambda)\]\$\$ where the first term is the expected joint log-density
and the second is the entropy of the approximation.

The `VariationalBayes` function provides eight general-purpose
algorithms spanning the modern variational inference landscape:

**Salimans1** implements the first algorithm of Salimans and Knowles
(2013), which uses only gradient information (no Hessian). The precision
matrix is updated via the outer product of gradients, \\P \leftarrow
(1-w)P + w \\ g g^\top\\, making it \\O(J+1)\\ Model evaluations per
iteration rather than the \\O(J^2/2)\\ required by Salimans2. This is
preferred when \\J\\ is large or when the posterior is not twice
differentiable. The same stochastic averaging and convergence scheme as
Salimans2 is used.

**Salimans2** implements the second algorithm of Salimans and Knowles
(2013), using both gradient and Hessian in a stochastic linear
regression framework for fixed-form Gaussian VB. This is the most
efficient method when \\J\\ is small and the posterior is smooth. The
step size is \\w = 1/\sqrt{T}\\ where \\T\\ is the total number of
iterations. Convergence is checked in the second half of iterations
after Polyak averaging stabilizes.

**ADVI.mf** implements the mean-field variant of Automatic
Differentiation Variational Inference (Kucukelbir et al., 2017). The
variational family is \\q(\theta) = \prod_j N(\theta_j; \mu_j,
\sigma_j^2)\\ with \\\sigma_j = \mathrm{softplus}(\omega_j) = \log(1 +
e^{\omega_j})\\ for unconstrained optimization. Gradients are estimated
via the reparameterization trick: draw \\\epsilon \sim N(0,I)\\, set
\\\theta = \mu + \sigma \odot \epsilon\\, and differentiate through the
sampling. Optimization uses Adam (Kingma and Ba, 2015) with default
learning rate 0.01.

**ADVI.fr** implements the full-rank variant of ADVI, replacing the
diagonal covariance with \\q(\theta) = N(\mu, LL^\top)\\ where \\L\\ is
a lower-triangular Cholesky factor. The diagonal of \\L\\ is
parameterized through softplus for positivity. This captures posterior
correlations at the cost of \\J(J+1)/2\\ variational parameters. A
warning is issued when \\J \> 50\\ due to the quadratic scaling.

**BBVI** implements Black-Box Variational Inference (Ranganath et al.,
2014), which uses the score function (REINFORCE) estimator rather than
the reparameterization trick. The gradient of the ELBO is estimated as
\\\nabla\_\lambda \mathcal{L} \approx (1/S) \sum_s f_s \nabla\_\lambda
\log q(\theta_s; \lambda)\\ where \\f_s = \log p(\theta_s, y) - \log
q(\theta_s; \lambda)\\. Control variates (per-component baselines)
reduce variance. The default uses \\S = 10\\ Monte Carlo samples per
iteration.

**SVGD** implements Stein Variational Gradient Descent (Liu and Wang,
2016), a particle-based method that maintains \\K\\ particles and
updates them via the functional gradient of KL divergence in a
reproducing kernel Hilbert space. The update for particle \\k\\ is
\\\theta_k \leftarrow \theta_k + \epsilon \\ \hat{\phi}^\*(\theta_k)\\
where \\\hat{\phi}^\*(\theta) = (1/K) \sum_j \[k(\theta_j, \theta)
\nabla\_{\theta_j} \log p(\theta_j, y) + \nabla\_{\theta_j} k(\theta_j,
\theta)\]\\ with RBF kernel and median heuristic bandwidth. The default
uses \\K = 20\\ particles. This method can capture multimodal posteriors
but requires \\K \times (J+1)\\ Model evaluations per iteration.

**NGD** implements Natural Gradient Variational Inference (Amari, 1998;
Khan and Lin, 2017) on the mean-field Gaussian family. Natural gradients
account for the Fisher information geometry of the variational family,
yielding faster convergence than ordinary gradients. The natural
parameters \\(\eta_1 = \Sigma^{-1}\mu,\\ \eta_2 =
-\tfrac{1}{2}\Sigma^{-1})\\ are updated with a Robbins-Monro schedule
\\\rho_t = (t + \tau)^{-\alpha}\\ where \\\tau = 10\\ and \\\alpha =
0.75\\. This is structurally related to Salimans2 but with an
information-geometric interpretation and different step-size scheduling.

**Pathfinder** implements the quasi-Newton variational inference method
of Zhang et al. (2022). Rather than iteratively optimizing variational
parameters, Pathfinder runs an L-BFGS optimization of the log-posterior,
stores the full trajectory \\\\\theta_k, \nabla_k\\\\, and at each
iterate reconstructs the L-BFGS inverse Hessian approximation \\H_k\\ to
define a Gaussian approximation \\q_k(\theta) = N(\theta_k, H_k)\\. The
ELBO is estimated at each iterate and the best approximation \\k^\* =
\arg\max_k \text{ELBO}\_k\\ is selected. In multi-path mode (the
default, with \\M = 4\\ paths), multiple independent L-BFGS runs are
launched with jittered initial values, the resulting draws are pooled
into a mixture \\q(\theta) = (1/M) \sum_m N(\theta_m, H_m)\\, importance
weights \\\log p(\theta) - \log q(\theta)\\ are computed, and Pareto
smoothed importance sampling (PSIS) is applied to produce the final
posterior draws. The Pareto \\\hat{k}\\ diagnostic quantifies the
quality of the approximation: values below 0.5 indicate reliable
estimates, 0.5–0.7 are acceptable, 0.7–1.0 are questionable, and values
above 1.0 signal that the approximation is unreliable. Pathfinder
typically requires 1–2 orders of magnitude fewer gradient evaluations
than ADVI while producing comparable sample quality. The number of paths
is controlled by `Covar` (integer, default 4); set `Covar = 1` for
single-path mode.

**Method selection guidance.** For small \\J\\ (fewer than ~20
parameters) with smooth posteriors, `Salimans2` or `NGD` are
recommended. For moderate \\J\\ (20–100), `ADVI.mf` or `Salimans1` avoid
the quadratic cost of the Hessian. `ADVI.fr` is appropriate when
posterior correlations are important and \\J \le 50\\. `BBVI` is the
most generic (requiring only sampling and density evaluation from
\\q\\), but has higher variance. `SVGD` is uniquely suited for
multimodal posteriors but is the most expensive per iteration.
`Pathfinder` is the recommended starting point for most problems: it is
fast, captures correlations via the L-BFGS inverse Hessian, and the
Pareto k diagnostic indicates when the approximation is unreliable.

## References

Amari, S. (1998). "Natural Gradient Works Efficiently in Learning".
*Neural Computation*, 10(2), p. 251–276.

Attias, H. (1999). "Inferring Parameters and Structure of Latent
Variable Models by Variational Bayes". In *Uncertainty in Artificial
Intelligence*.

Attias, H. (2000). "A Variational Bayesian Framework for Graphical
Models". In *Neural Information Processing Systems*.

Ghahramani, Z. and Beal, M. (2001). "Propagation Algorithms for
Variational Bayesian Learning". In *Neural Information Processing
Systems*, p. 507–513.

Jaakkola, T. (1997). "Variational Methods for Inference and Estimation
in Graphical Models". PhD thesis, Massachusetts Institute of Technology.

Khan, M.E. and Lin, W. (2017). "Conjugate-Computation Variational
Inference: Converting Variational Inference in Non-Conjugate Models to
Inferences in Conjugate Models". In *Proceedings of the 20th
International Conference on Artificial Intelligence and Statistics
(AISTATS)*.

Kingma, D.P. and Ba, J. (2015). "Adam: A Method for Stochastic
Optimization". In *Proceedings of the 3rd International Conference on
Learning Representations (ICLR)*.

Kucukelbir, A., Tran, D., Ranganath, R., Gelman, A., and Blei, D.M.
(2017). "Automatic Differentiation Variational Inference". *Journal of
Machine Learning Research*, 18(14), p. 1–45.

Liu, Q. and Wang, D. (2016). "Stein Variational Gradient Descent: A
General Purpose Bayesian Inference Algorithm". In *Advances in Neural
Information Processing Systems 29 (NeurIPS)*.

Neal, R. and Hinton, G. (1999). "A View of the EM Algorithm that
Justifies Incremental, Sparse, and Other Variants". In Learning in
Graphical Models, p. 355–368. MIT Press, 1999.

Ranganath, R., Gerrish, S., and Blei, D.M. (2014). "Black Box
Variational Inference". In *Proceedings of the 17th International
Conference on Artificial Intelligence and Statistics (AISTATS)*.

Salimans, T. and Knowles, D.A. (2013). "Fixed-Form Variational Posterior
Approximation through Stochastic Linear Regression". *Bayesian
Analysis*, 8(4), p. 837–882.

Saul, L. and Jordan, M. (1996). "Exploiting Tractable Substructures in
Intractable Networks". *Neural Information Processing Systems*.

Saul, L., Jaakkola, T., and Jordan, M. (1996). "Mean Field Theory for
Sigmoid Belief Networks". *Journal of Artificial Intelligence Research*,
4, p. 61–76.

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
"Pareto Smoothed Importance Sampling". *Journal of Machine Learning
Research*, 25(72), p. 1–58.

Wiegerinck, W. (2000). "Variational Approximations Between Mean Field
Theory and the Junction Tree Algorithm". In *Uncertainty in Artificial
Intelligence*.

Xing, E., Jordan, M., and Russell, S. (2003). "A Generalized Mean Field
Algorithm for Variational Inference in Exponential Families". In
*Uncertainty in Artificial Intelligence*.

Zhang, L., Carpenter, B., Gelman, A., and Vehtari, A. (2022).
"Pathfinder: Parallel Quasi-Newton Variational Inference". *Journal of
Machine Learning Research*, 23(306), p. 1–49.

## See also

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)

# Setup: regression on demonsnacks data
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

mon.names <- "mu[1]"
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

Model <- function(parm, Data) {
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
     LP <- LL + beta.prior + sigma.prior
     list(LP=LP, Dev=-2*LL, Monitor=mu[1],
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     }

IV <- rep(0, J+1)

# Fixed-form VB (Salimans2, default)
Fit1 <- VariationalBayes(Model, IV, Data=MyData, Iterations=1000,
     Method="Salimans2", Stop.Tolerance=1e-3)
print(Fit1)

# ADVI mean-field
Fit2 <- VariationalBayes(Model, IV, Data=MyData, Iterations=1000,
     Method="ADVI.mf", Stop.Tolerance=1e-3)
print(Fit2)

# SVGD with 20 particles
Fit3 <- VariationalBayes(Model, IV, Data=MyData, Iterations=500,
     Method="SVGD", Covar=20, Stop.Tolerance=1e-3)
print(Fit3)

# Pathfinder with 4 paths (default)
Fit4 <- VariationalBayes(Model, IV, Data=MyData, Iterations=1000,
     Method="Pathfinder", Covar=4, Stop.Tolerance=1e-5)
print(Fit4)
summary(Fit4)
} # }
```
