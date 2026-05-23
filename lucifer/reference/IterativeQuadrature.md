# Iterative quadrature

The `IterativeQuadrature` function iteratively approximates the first
two moments of marginal posterior distributions of a Bayesian model with
deterministic integration.

## Usage

``` r
IterativeQuadrature(
  Model,
  parm,
  Data,
  Covar = NULL,
  Iterations = 100,
  Algorithm = "CAGH",
  Specs = NULL,
  Samples = 1000,
  sir = TRUE,
  Stop.Tolerance = c(1e-05, 1e-15),
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  this required argument receives the model from a user-defined
  function. The user-defined function is where the model is specified.
  `IterativeQuadrature` passes two arguments to the model function,
  `parms` and `Data`. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function and "lucifer tutorial" vignette.

- parm:

  this argument requires a vector of initial values equal in length to
  the number of parameters. `IterativeQuadrature` will attempt to
  approximate these initial values for the parameters as means (or
  posterior modes) of normal integrals. The
  [`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md)
  function may be used to randomly generate initial values. Parameters
  must be continuous.

- Data:

  this required argument accepts a list of data. The list of data must
  include `mon.names` which contains monitored variable names, and
  `parm.names` which contains parameter names.

- Covar:

  this argument accepts a \\J \times J\\ covariance matrix for \\J\\
  initial values. When a covariance matrix is not supplied, a scaled
  identity matrix is used.

- Iterations:

  this argument accepts an integer that determines the number of
  iterations that `IterativeQuadrature` will attempt to approximate the
  posterior with normal integrals. `Iterations` defaults to 100.
  `IterativeQuadrature` will stop before this number of iterations if
  the tolerance is less than or equal to the `Stop.Tolerance` criterion.
  The required amount of computer memory increases with `Iterations`. If
  computer memory is exceeded, then all will be lost.

- Algorithm:

  this optional argument accepts a quoted string that specifies the
  iterative quadrature algorithm. The default method is
  `Algorithm="CAGH"`. Options include `"AGHSG"` for Adaptive
  Gauss-Hermite Sparse Grid, and `"CAGH"` for Componentwise Adaptive
  Gaussian-Hermite.

- Specs:

  this argument accepts a list of specifications for an algorithm.

- Samples:

  this argument indicates the number of posterior samples to be taken
  with sampling importance resampling via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function, which occurs only when `sir=TRUE`. Note that the number of
  samples should increase with the number and intercorrelations of the
  parameters.

- sir:

  this logical argument indicates whether or not Sampling Importance
  Resampling (SIR) is conducted via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function to draw independent posterior samples. This argument defaults
  to `TRUE`. Even when `TRUE`, posterior samples are drawn only when
  `IterativeQuadrature` has converged. Posterior samples are required
  for many other functions, including `plot.iterquad` and
  `predict.iterquad`. Less time can be spent on sampling by increasing
  `CPUs`, if available, which parallelizes the sampling.

- Stop.Tolerance:

  this argument accepts a vector of two positive numbers, and defaults
  to `c(1e-5, 1e-15)`. Tolerance is calculated each iteration, and the
  criteria varies by algorithm. The algorithm is considered to have
  converged to the user-specified `Stop.Tolerance` when the tolerance is
  less than or equal to the value of `Stop.Tolerance`, and the algorithm
  terminates at the end of the current iteration. Unless stated
  otherwise, the first element is the stop tolerance for the change in
  \\\mu\\, the second element is the stop tolerance for the change in
  mean integration error, and the first tolerance must be met before the
  second tolerance is considered.

- CPUs:

  this argument accepts an integer that specifies the number of central
  processing units (CPUs) of the multicore computer or computer cluster.
  This argument defaults to `CPUs=1`, in which parallel processing does
  not occur. When multiple CPUs are specified, model function
  evaluations are parallelized across the nodes, and sampling with
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md) is
  parallelized when `sir=TRUE`.

- Type:

  this argument specifies the type of parallel processing to perform,
  accepting either `Type="PSOCK"` or `Type="MPI"`.

## Value

`IterativeQuadrature` returns an object of class `iterquad` that is a
list with the following components:

- Algorithm:

  the name of the iterative quadrature algorithm.

- Call:

  the matched call of `IterativeQuadrature`.

- Converged:

  a logical indicator of whether or not `IterativeQuadrature` converged
  within the specified `Iterations` according to the supplied
  `Stop.Tolerance` criterion. Convergence does not indicate that the
  global maximum has been found, but only that the tolerance was less
  than or equal to the `Stop.Tolerance` criteria.

- Covar:

  the estimated covariance matrix. The `Covar` matrix may be scaled and
  input into the `Covar` argument of the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  or [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
  function for further estimation. To scale this matrix for use with
  Laplace's Demon or PMC, multiply it by \\2.38^2/d\\, where \\d\\ is
  the number of initial values.

- Deviance:

  a vector of the iterative history of the deviance in the
  `IterativeQuadrature` function, as it sought convergence.

- History:

  a matrix of the iterative history of the parameters in the
  `IterativeQuadrature` function, as it sought convergence.

- Initial.Values:

  the vector of initial values that was originally given to
  `IterativeQuadrature` in the `parm` argument.

- LML:

  an approximation of the logarithm of the marginal likelihood of the
  data (see the
  [`LML`](https://robustecologies.github.io/lucifer/reference/LML.md)
  function for more information). When the model has converged and
  `sir=TRUE`, the NSIS method is used. When the model has converged and
  `sir=FALSE`, the LME method is used. This is the logarithmic form of
  equation 4 in Lewis and Raftery (1997). As a rough estimate of Kass
  and Raftery (1995), the LME-based LML is worrisome when the sample
  size of the data is less than five times the number of parameters, and
  `LML` should be adequate in most problems when the sample size of the
  data exceeds twenty times the number of parameters (p. 778). The LME
  is inappropriate with hierarchical models. However `LML` is estimated,
  it is useful for comparing multiple models with the `BayesFactor`
  function.

- LP.Final:

  the final scalar value for the logarithm of the unnormalized joint
  posterior density.

- LP.Initial:

  the initial scalar value for the logarithm of the unnormalized joint
  posterior density.

- LPw:

  the latest matrix of the logarithm of the unnormalized joint posterior
  density. It is weighted and normalized so that each column sums to
  one.

- M:

  the final \\N \times J\\ matrix of quadrature weights that have been
  corrected for non-standard normal distributions, where \\N\\ is the
  number of nodes and \\J\\ is the number of parameters.

- Minutes:

  the number of minutes that `IterativeQuadrature` was running,
  including the initial checks as well as drawing posterior samples and
  creating summaries.

- Monitor:

  when `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix. The
  rows of the matrix are the samples, and the columns are the monitored
  variables.

- N:

  the final number of nodes.

- Posterior:

  when `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix. The
  rows of the matrix are the samples, and the columns are the
  parameters.

- Summary1:

  a summary matrix that summarizes the point-estimated posterior means.
  Uncertainty around the posterior means is estimated from the
  covariance matrix. Rows are parameters. The following columns are
  included: Mean, SD (Standard Deviation), LB (Lower Bound), and UB
  (Upper Bound). The bounds constitute a 95% probability interval.

- Summary2:

  a summary matrix that summarizes the posterior samples drawn with
  sampling importance resampling
  ([`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md))
  when `sir=TRUE`, given the point-estimated posterior modes and the
  covariance matrix. Rows are parameters. The following columns are
  included: Mean, SD (Standard Deviation), LB (Lower Bound), and UB
  (Upper Bound). The bounds constitute a 95% probability interval.

- Tolerance.Final:

  the last `Tolerance` of the `IterativeQuadrature` algorithm.

- Tolerance.Stop:

  the `Stop.Tolerance` criteria.

- Z:

  the final \\N \times J\\ matrix of the conditional mean, where \\N\\
  is the number of nodes and \\J\\ is the number of parameters.

## Details

Quadrature is a historical term in mathematics that means determining
area. Mathematicians of ancient Greece, according to the Pythagorean
doctrine, understood determination of area of a figure as the process of
geometrically constructing a square having the same area (squaring),
thus the name quadrature for this process.

In medieval Europe, quadrature meant the calculation of area by any
method. With the invention of integral calculus, quadrature has been
applied to the computation of a univariate definite integral. Numerical
integration is a broad family of algorithms for calculating the
numerical value of a definite integral. Numerical quadrature is a
synonym for quadrature applied to one-dimensional integrals.
Multivariate quadrature, also called cubature, is the application of
quadrature to multidimensional integrals.

A quadrature rule is an approximation of the definite integral of a
function, usually stated as a weighted sum of function values at
specified points within the domain of integration. The specified points
are referred to as abscissae, abscissas, integration points, or nodes,
and have associated weights. The calculation of the nodes and weights of
the quadrature rule differs by the type of quadrature. There are
numerous types of quadrature algorithms. Bayesian forms of quadrature
usually use Gauss-Hermite quadrature (Naylor and Smith, 1982), and
placing a Gaussian Process on the function is a common extension
(O'Hagan, 1991; Rasmussen and Ghahramani, 2003) that is called "Bayesian
Quadrature". Often, these and other forms of quadrature are also
referred to as model-based integration.

Gauss-Hermite quadrature uses Hermite polynomials to calculate the rule.
However, there are two versions of Hermite polynomials, which result in
different kernels in different fields. In physics, the kernel is
`exp(-x^2)`, while in probability the kernel is `exp(-x^2/2)`. The
weights are a normal density. If the parameters of the normal
distribution, \\\mu\\ and \\\sigma^2\\, are estimated from data, then it
is referred to as adaptive Gauss-Hermite quadrature, and the parameters
are the conditional mean and conditional variance. Outside of
Gauss-Hermite quadrature, adaptive quadrature implies that a difficult
range in the integrand is subdivided with more points until it is
well-approximated. Gauss-Hermite quadrature performs well when the
integrand is smooth, and assumes normality or multivariate normality.
Adaptive Gauss-Hermite quadrature has been demonstrated to outperform
Gauss-Hermite quadrature in speed and accuracy.

A goal in quadrature is to minimize integration error, which is the
error between the evaluations and the weights of the rule. Therefore, a
goal in Bayesian Gauss-Hermite quadrature is to minimize integration
error while approximating a marginal posterior distribution that is
assumed to be smooth and normally-distributed. This minimization often
occurs by increasing the number of nodes until a change in mean
integration error is below a tolerance, rather than minimizing
integration error itself, since the target may be only approximately
normally distributed, or minimizing the sum of integration error, which
would change with the number of nodes.

To approximate integrals in multiple dimensions, one approach applies
\\N\\ nodes of a univariate quadrature rule to multiple dimensions
(using the
[`GaussHermiteCubeRule`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
function for example) via the product rule, which results in many more
multivariate nodes. This requires the number of function evaluations to
grow exponentially as dimension increases. Multidimensional quadrature
is usually limited to less than ten dimensions, both due to the number
of nodes required, and because the accuracy of multidimensional
quadrature algorithms decreases as the dimension increases. Three
methods may overcome this curse of dimensionality in varying degrees:
componentwise quadrature, sparse grids, and Monte Carlo.

Componentwise quadrature is the iterative application of univariate
quadrature to each parameter. It is applicable with high-dimensional
models, but sacrifices the ability to calculate the conditional
covariance matrix, and calculates only the variance of each parameter.

Sparse grids were originally developed by Smolyak for multidimensional
quadrature. A sparse grid is based on a one-dimensional quadrature rule.
Only a subset of the nodes from the product rule is included, and the
weights are appropriately rescaled. Although a sparse grid is more
efficient because it reduces the number of nodes to achieve the same
accuracy, the user must contend with increasing the accuracy of the
grid, and it remains inapplicable to high-dimensional integrals.

Monte Carlo is a large family of sampling-based algorithms. O'Hagan
(1987) asserts that Monte Carlo is frequentist, inefficient, regards
irrelevant information, and disregards relevant information. Quadrature,
he maintains (O'Hagan, 1992), is the most Bayesian approach, and also
the most efficient. In high dimensions, he concedes, a popular subset of
Monte Carlo algorithms is currently the best for cheap model function
evaluations. These algorithms are called Markov chain Monte Carlo
(MCMC). High-dimensional models with expensive model evaluation
functions, however, are not well-suited to MCMC. A large number of MCMC
algorithms is available in the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function.

Following are some reasons to consider iterative quadrature rather than
MCMC. Once an MCMC sampler finds equilibrium, it must then draw enough
samples to represent all targets. Iterative quadrature does not need to
continue drawing samples. Multivariate quadrature is consistently
reported as more efficient than MCMC when its assumptions hold, though
multivariate quadrature is limited to small dimensions. High-dimensional
models therefore default to MCMC, between the two. Componentwise
quadrature algorithms like CAGH, however, may also be more efficient
with clock-time than MCMC in high dimensions, especially against
componentwise MCMC algorithms. Another reason to consider iterative
quadrature is that assessing convergence in MCMC is a difficult topic,
but not for iterative quadrature. A user of iterative quadrature does
not have to contend with effective sample size and autocorrelation,
assessing stationarity, acceptance rates, diminishing adaptation, etc.
Stochastic sampling in MCMC is less efficient when samples occur in
close proximity (such as when highly autocorrelated), whereas in
quadrature the nodes are spread out by design.

In general, the conditional means and conditional variances progress
smoothly to the target in multidimensional quadrature. For componentwise
quadrature, movement to the target is not smooth, and often resembles a
Markov chain or optimization algorithm.

Iterative quadrature is often applied after
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
to obtain a more reliable estimate of parameter variance or covariance
than the negative inverse of the
[`Hessian`](https://robustecologies.github.io/lucifer/reference/Matrices.md)
matrix of second derivatives, which is suitable only when the contours
of the logarithm of the unnormalized joint posterior density are
approximately ellipsoidal (Naylor and Smith, 1982, p. 224).

When `Algorithm="AGH"`, the Naylor and Smith (1982) algorithm is used.
The AGH algorithm uses multivariate quadrature with the physicist's (not
the probabilist's) kernel.

There are four algorithm specifications: `N` is the number of univariate
nodes, `Nmax` is the maximum number of univariate nodes, `Packages`
accepts any package required for the model function when parallelized,
and `Dyn.libs` accepts dynamic libraries for parallelization, if
required. The number of univariate nodes begins at \\N\\ and increases
by one each iteration. The number of multivariate nodes grows quickly
with \\N\\. Naylor and Smith (1982) recommend beginning with as few
nodes as \\N=3\\. Any of the following events will cause \\N\\ to
increase by 1 when \\N\\ is less than `Nmax`:

- All LP weights are zero (and non-finite weights are set to zero)

- \\\mu\\ does not result in an increase in LP

- All elements in \\\Sigma\\ are not finite

- The square root of the sum of the squared changes in \\\mu\\ is less
  than or equal to the `Stop.Tolerance`

Tolerance includes two metrics: change in mean integration error and
change in parameters. Including the change in parameters for tolerance
was not mentioned in Naylor and Smith (1982).

Naylor and Smith (1982) consider a transformation due to correlation.
This is not included here.

The AGH algorithm does not currently handle constrained parameters, such
as with the
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md)
function. If a parameter is constrained and changes during a model
evaluation, this changes the node and the multivariate weight. This is
currently not corrected.

An advantage of AGH over componentwise adaptive quadrature is that AGH
estimates covariance, where a componentwise algorithm ignores it. A
disadvantage of AGH over a componentwise algorithm is that the number of
nodes increases so quickly with dimension, that AGH is limited to
small-dimensional models.

When `Algorithm="AGHSG"`, the Naylor and Smith (1982) algorithm is
applied to a sparse grid, rather than a traditional multivariate
quadrature rule. This is identical to the AGH algorithm above, except
that a sparse grid replaces the multivariate quadrature rule.

The sparse grid reduces the number of nodes. The cost of reducing the
number of nodes is that the user must consider the accuracy, \\K\\.

There are four algorithm specifications: `K` is the accuracy (as a
positive integer), `Kmax` is the maximum accuracy, `Packages` accepts
any package required for the model function when parallelized, and
`Dyn.libs` accepts dynamic libraries for parallelization, if required.
These arguments represent accuracy rather than the number of univariate
nodes, but otherwise are similar to the AGH algorithm.

When `Algorithm="CAGH"`, a componentwise version of the adaptive
Gauss-Hermite quadrature of Naylor and Smith (1982) is used. Each
iteration, each marginal posterior distribution is approximated
sequentially, in a random order, with univariate quadrature. The
conditional mean and conditional variance are also approximated each
iteration, making it an adaptive algorithm.

There are four algorithm specifications: `N` is the number of nodes,
`Nmax` is the maximum number of nodes, `Packages` accepts any package
required for the model function when parallelized, and `Dyn.libs`
accepts dynamic libraries for parallelization, if required. The number
of nodes begins at \\N\\. All parameters have the same number of nodes.
Any of the following events will cause \\N\\ to increase by 1 when \\N\\
is less than `Nmax`, and these conditions refer to all parameters (not
individually):

- Any LP weights are not finite

- All LP weights are zero

- \\\mu\\ does not result in an increase in LP

- The square root of the sum of the squared changes in \\\mu\\ is less
  than or equal to the `Stop.Tolerance`

It is recommended to begin with `N=3` and set `Nmax` between 10 and 100.
As long as CAGH does not experience problematic weights, and as long as
CAGH is improving LP with \\\mu\\, the number of nodes does not
increase. When CAGH becomes either universally problematic or
universally stable, then \\N\\ slowly increases until the sum of both
the mean integration error and the sum of the squared changes in \\\mu\\
is less than the `Stop.Tolerance` for two consecutive iterations.

If the highest LP occurs at the lowest or highest node, then the value
at that node becomes the conditional mean, rather than calculating it
from all weighted samples; this facilitates movement when the current
integral is poorly centered toward a well-centered integral. If all
weights are zero, then a random proposal is generated with a small
variance.

Tolerance includes two metrics: change in mean integration error and
change in parameters, as the square root of the sum of the squared
differences.

When a parameter constraint is encountered, the node and weight of the
quadrature rule is recalculated.

An advantage of CAGH over multidimensional adaptive quadrature is that
CAGH may be applied in large dimensions. Disadvantages of CAGH are that
only variance, not covariance, is estimated, and ignoring covariance may
be problematic.

## References

Naylor, J.C. and Smith, A.F.M. (1982). "Applications of a Method for the
Efficient Computation of Posterior Distributions". *Applied Statistics*,
31(3), p. 214–225.

O'Hagan, A. (1987). "Monte Carlo is Fundamentally Unsound". *The
Statistician*, 36, p. 247–249.

O'Hagan, A. (1991). "Bayes-Hermite Quadrature". *Journal of Statistical
Planning and Inference*, 29, p. 245–260.

O'Hagan, A. (1992). "Some Bayesian Numerical Analysis". In Bernardo,
J.M., Berger, J.O., David, A.P., and Smith, A.F.M., editors, *Bayesian
Statistics*, 4, p. 356–363, Oxford University Press.

Rasmussen, C.E. and Ghahramani, Z. (2003). "Bayesian Monte Carlo". In
Becker, S. and Obermayer, K., editors, *Advances in Neural Information
Processing Systems*, 15, MIT Press, Cambridge, MA.

## See also

[`GaussHermiteCubeRule`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
[`GaussHermiteQuadRule`](https://robustecologies.github.io/lucifer/reference/Math.md),
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md),
[`Hermite`](https://robustecologies.github.io/lucifer/reference/Math.md),
[`Hessian`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md),
[`SparseGrid`](https://robustecologies.github.io/lucifer/reference/Matrices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10] + 1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])

mon.names <- "mu[1]"
parm.names <- as.parm.names(list(beta = rep(0, J), sigma = 0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) c(rnorm(Data$J), runif(1))
MyData <- list(J = J, PGF = PGF, X = X, mon.names = mon.names,
     parm.names = parm.names, pos.beta = pos.beta,
     pos.sigma = pos.sigma, y = y)

Model <- function(parm, Data) {
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
     parm[Data$pos.sigma] <- sigma
     beta.prior <- sum(dnormv(beta, 0, 1000, log = TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
     LP <- LL + beta.prior + sigma.prior
     list(LP = LP, Dev = -2*LL, Monitor = mu[1],
          yhat = rnorm(length(mu), mu, sigma), parm = parm)
}

Initial.Values <- rep(0, J + 1)

# Adaptive Gauss-Hermite
Fit <- IterativeQuadrature(Model, Initial.Values, MyData, Covar = NULL,
     Iterations = 100, Algorithm = "AGH",
     Specs = list(N = 5, Nmax = 7, Packages = NULL, Dyn.libs = NULL),
     CPUs = 1)

# Adaptive Gauss-Hermite Sparse Grid
Fit <- IterativeQuadrature(Model, Initial.Values, MyData, Covar = NULL,
     Iterations = 100, Algorithm = "AGHSG",
     Specs = list(K = 5, Kmax = 7, Packages = NULL, Dyn.libs = NULL),
     CPUs = 1)

# Componentwise Adaptive Gauss-Hermite
Fit <- IterativeQuadrature(Model, Initial.Values, MyData, Covar = NULL,
     Iterations = 100, Algorithm = "CAGH",
     Specs = list(N = 3, Nmax = 10, Packages = NULL, Dyn.libs = NULL),
     CPUs = 1)

print(Fit)
summary(Fit)
plot(Fit, MyData, PDF = FALSE)

Pred <- predict(Fit, Model, MyData, CPUs = 1)
summary(Pred, Discrep = "Chi-Square")
plot(Pred, Style = "Fitted")
} # }
```
