# Laplace approximation

The `LaplaceApproximation` function deterministically maximizes the
logarithm of the unnormalized joint posterior density with one of
several optimization algorithms. The goal of Laplace Approximation is to
estimate the posterior mode and variance of each parameter. This
function is useful for optimizing initial values and estimating a
covariance matrix to be input into the
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
function, or sometimes for model estimation in its own right.

## Usage

``` r
LaplaceApproximation(
  Model,
  parm,
  Data,
  Interval = 1e-06,
  Iterations = 100,
  Method = "SPG",
  Samples = 1000,
  CovEst = "Hessian",
  sir = TRUE,
  Stop.Tolerance = 1e-05,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- Model:

  this required argument receives the model from a user-defined
  function. The user-defined function is where the model is specified.
  `LaplaceApproximation` passes two arguments to the model function,
  `parms` and `Data`. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function and "lucifer tutorial" vignette.

- parm:

  this argument requires a vector of initial values equal in length to
  the number of parameters. `LaplaceApproximation` will attempt to
  optimize these initial values for the parameters, where the optimized
  values are the posterior modes, for later use with the
  [`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
  or the
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
  function. The
  [`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md)
  function may be used to randomly generate initial values. Parameters
  must be continuous.

- Data:

  this required argument accepts a list of data. The list of data must
  include `mon.names` which contains monitored variable names, and
  `parm.names` which contains parameter names. `LaplaceApproximation`
  must be able to determine the sample size of the data, and will look
  for a scalar sample size variable `n` or `N`. If not found, it will
  look for variable `y` or `Y`, and attempt to take its number of rows
  as sample size. `LaplaceApproximation` needs to determine sample size
  due to the asymptotic nature of this method. Sample size should be at
  least \\\sqrt{J}\\ with \\J\\ exchangeable parameters.

- Interval:

  this argument receives an interval for estimating approximate
  gradients. The logarithm of the unnormalized joint posterior density
  of the Bayesian model is evaluated at the current parameter value, and
  again at the current parameter value plus this interval.

- Iterations:

  this argument accepts an integer that determines the number of
  iterations that `LaplaceApproximation` will attempt to maximize the
  logarithm of the unnormalized joint posterior density. `Iterations`
  defaults to 100. `LaplaceApproximation` will stop before this number
  of iterations if the tolerance is less than or equal to the
  `Stop.Tolerance` criterion. The required amount of computer memory
  increases with `Iterations`. If computer memory is exceeded, then all
  will be lost.

- Method:

  this optional argument accepts a quoted string that specifies the
  method used for Laplace Approximation. The default method is
  `Method="SPG"`. Options include `"AGA"` for adaptive gradient ascent,
  `"BFGS"` for the Broyden-Fletcher-Goldfarb-Shanno algorithm, `"BHHH"`
  for the algorithm of Berndt et al., `"CG"` for conjugate gradient,
  `"DFP"` for the Davidon-Fletcher-Powell algorithm, `"HAR"` for
  adaptive hit-and-run, `"HJ"` for Hooke-Jeeves, `"LBFGS"` for
  limited-memory BFGS, `"LM"` for Levenberg-Marquardt, `"NM"` for
  Nelder-Mead, `"NR"` for Newton-Raphson, `"PSO"` for Particle Swarm
  Optimization, `"Rprop"` for resilient backpropagation, `"SGD"` for
  Stochastic Gradient Descent, `"SOMA"` for the Self-Organizing
  Migration Algorithm, `"SPG"` for Spectral Projected Gradient, `"SR1"`
  for Symmetric Rank-One, and `"TR"` for Trust Region.

- Samples:

  this argument indicates the number of posterior samples to be taken
  with sampling importance resampling via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function, which occurs only when `sir=TRUE`. Note that the number of
  samples should increase with the number and intercorrelations of the
  parameters.

- CovEst:

  this argument accepts a quoted string that indicates how the
  covariance matrix is estimated after the model finishes. This
  covariance matrix is used to obtain the standard deviation of each
  parameter, and may also be used for posterior sampling via Sampling
  Importance Resampling (SIR) (see the `sir` argument below), if
  converged. By default, the covariance matrix is approximated as the
  negative inverse of the `"Hessian"` matrix of second derivatives,
  estimated with Richardson extrapolation. Alternatives include
  `CovEst="Identity"`, `CovEst="OPG"`, or `CovEst="Sandwich"`. When
  `CovEst="Identity"`, the covariance matrix is not estimated, and is
  merely assigned an identity matrix. When `LaplaceApproximation` is
  performed internally by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  an identity matrix is returned and scaled. When `CovEst="OPG"`, the
  covariance matrix is approximated with the inverse of the sum of the
  outer products of the gradient, which requires `X`, and either `y` or
  `Y` in the list of data. For OPG, a partial derivative is taken for
  each row in `X`, and each element in `y` or row in `Y`. Therefore,
  this requires \\N + NJ\\ model evaluations for a data set with \\N\\
  records and \\J\\ variables. The OPG method is an asymptotic
  approximation of the Hessian, and usually requires fewer calculations
  with a small data set, or more with large data sets. Both methods
  require a matrix inversion, which becomes costly as dimension grows.
  The Richardson-based Hessian method is more accurate, but requires
  more calculation in large dimensions. An alternative approach to
  obtaining covariance is to use the
  [`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md)
  on the data, or sample the posterior with iterative quadrature
  ([`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md)),
  MCMC
  ([`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)),
  or
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

- sir:

  this logical argument indicates whether or not Sampling Importance
  Resampling (SIR) is conducted via the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function to draw independent posterior samples. This argument defaults
  to `TRUE`. Even when `TRUE`, posterior samples are drawn only when
  `LaplaceApproximation` has converged. Posterior samples are required
  for many other functions, including `plot.laplace` and
  `predict.laplace`. The only time that it is advantageous for
  `sir=FALSE` is when `LaplaceApproximation` is used to help the initial
  values for
  [`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
  or
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
  and it is unnecessary for time to be spent on sampling. Less time can
  be spent on sampling by increasing `CPUs`, which parallelizes the
  sampling.

- Stop.Tolerance:

  this argument accepts any positive number and defaults to 1.0E-5.
  Tolerance is calculated each iteration, and the criteria varies by
  algorithm. The algorithm is considered to have converged to the
  user-specified `Stop.Tolerance` when the tolerance is less than or
  equal to the value of `Stop.Tolerance`, and the algorithm terminates
  at the end of the current iteration. Often, multiple criteria are
  used, in which case the maximum of all criteria becomes the tolerance.
  For example, when partial derivatives are taken, it is commonly
  required that the Euclidean norm of the partial derivatives is a
  criterion, and another common criterion is the Euclidean norm of the
  differences between the current and previous parameter values. Several
  algorithms have other, specific tolerances.

- CPUs:

  this argument accepts an integer that specifies the number of central
  processing units (CPUs) of the multicore computer or computer cluster.
  This argument defaults to `CPUs=1`, in which parallel processing does
  not occur. Parallelization occurs only for sampling with
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  when `sir=TRUE`.

- Type:

  this argument specifies the type of parallel processing to perform,
  accepting either `Type="PSOCK"` or `Type="MPI"`.

## Value

`LaplaceApproximation` returns an object of class `laplace` that is a
list with the following components:

- Call:

  the matched call of `LaplaceApproximation`.

- Converged:

  a logical indicator of whether or not `LaplaceApproximation` converged
  within the specified `Iterations` according to the supplied
  `Stop.Tolerance` criterion. Convergence does not indicate that the
  global maximum has been found, but only that the tolerance was less
  than or equal to the `Stop.Tolerance` criterion.

- Covar:

  the covariance matrix estimated according to the `CovEst` argument.
  The `Covar` matrix may be scaled and input into the `Covar` argument
  of the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  or [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
  function for further estimation, or the diagonal of this matrix may be
  used to represent the posterior variance of the parameters, provided
  the algorithm converged and matrix inversion was successful. To scale
  this matrix for use with Lucifer or PMC, multiply it by \\2.38^2/d\\,
  where \\d\\ is the number of initial values.

- Deviance:

  a vector of the iterative history of the deviance in the
  `LaplaceApproximation` function, as it sought convergence.

- History:

  a matrix of the iterative history of the parameters in the
  `LaplaceApproximation` function, as it sought convergence.

- Initial.Values:

  the vector of initial values that was originally given to
  `LaplaceApproximation` in the `parm` argument.

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

- Minutes:

  the number of minutes that `LaplaceApproximation` was running,
  including the initial checks as well as drawing posterior samples and
  creating summaries.

- Monitor:

  when `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix. The
  rows of the matrix are the samples, and the columns are the monitored
  variables.

- Posterior:

  when `sir=TRUE`, a number of independent posterior samples equal to
  `Samples` is taken, and the draws are stored here as a matrix. The
  rows of the matrix are the samples, and the columns are the
  parameters.

- Step.Size.Final:

  the final, scalar `Step.Size` value at the end of the
  `LaplaceApproximation` algorithm.

- Step.Size.Initial:

  the initial, scalar `Step.Size`.

- Summary1:

  a summary matrix that summarizes the point-estimated posterior modes.
  Uncertainty around the posterior modes is estimated from the
  covariance matrix. Rows are parameters. The following columns are
  included: Mode, SD (Standard Deviation), LB (Lower Bound), and UB
  (Upper Bound). The bounds constitute a 95% probability interval.

- Summary2:

  a summary matrix that summarizes the posterior samples drawn with
  sampling importance resampling
  ([`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md))
  when `sir=TRUE`, given the point-estimated posterior modes and the
  covariance matrix. Rows are parameters. The following columns are
  included: Mode, SD (Standard Deviation), LB (Lower Bound), and UB
  (Upper Bound). The bounds constitute a 95% probability interval.

- Tolerance.Final:

  the last `Tolerance` of the `LaplaceApproximation` algorithm.

- Tolerance.Stop:

  the `Stop.Tolerance` criterion.

## Details

The Laplace Approximation or Laplace Method is a family of asymptotic
techniques used to approximate integrals. Laplace's method accurately
approximates unimodal posterior moments and marginal posterior
distributions in many cases. Since it is not applicable in all cases, it
is recommended here that Laplace Approximation is used cautiously in its
own right, or preferably, it is used before MCMC.

After introducing the Laplace Approximation (Laplace, 1774, p. 366–367),
a proof was published later (Laplace, 1814) as part of a mathematical
system of inductive reasoning based on probability. Laplace used this
method to approximate posterior moments.

Since its introduction, the Laplace Approximation has been applied
successfully in many disciplines. In the 1980s, the Laplace
Approximation experienced renewed interest, especially in statistics,
and some improvements in its implementation were introduced (Tierney et
al., 1986; Tierney et al., 1989). Only since the 1980s has the Laplace
Approximation been seriously considered by statisticians in practical
applications.

There are many variations of Laplace Approximation, with an effort
toward replacing Markov chain Monte Carlo (MCMC) algorithms as the
dominant form of numerical approximation in Bayesian inference. The
run-time of Laplace Approximation is a little longer than Maximum
Likelihood Estimation (MLE), usually shorter than variational Bayes, and
much shorter than MCMC (Azevedo and Shachter, 1994).

The speed of Laplace Approximation depends on the optimization algorithm
selected, and typically involves many evaluations of the objective
function per iteration (where an MCMC algorithm with a multivariate
proposal usually evaluates once per iteration), making many MCMC
algorithms faster per iteration. The attractiveness of Laplace
Approximation is that it typically improves the objective function
better than iterative quadrature, MCMC, and PMC when the parameters are
in low-probability regions. Laplace Approximation is also typically
faster than MCMC and PMC because it is seeking point-estimates, rather
than attempting to represent the target distribution with enough
simulation draws. Laplace Approximation extends MLE, but shares similar
limitations, such as its asymptotic nature with respect to sample size
and that marginal posterior distributions are Gaussian. Bernardo and
Smith (2000) note that Laplace Approximation is an attractive family of
numerical approximation algorithms, and will continue to develop.

`LaplaceApproximation` seeks a global maximum of the logarithm of the
unnormalized joint posterior density. The approach differs by `Method`.
The
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function uses the `LaplaceApproximation` algorithm to optimize initial
values and save time for the user.

Most optimization algorithms assume that the logarithm of the
unnormalized joint posterior density is defined and differentiable. Some
methods calculate an approximate gradient for each initial value as the
difference in the logarithm of the unnormalized joint posterior density
due to a slight increase in the parameter.

When `Method="AGA"`, the direction and distance for each parameter is
proposed based on an approximate truncated gradient and an adaptive step
size. The step size parameter, which is often plural and called rate
parameters in other literature, is adapted each iteration with the
univariate version of the Robbins-Monro stochastic approximation in
Garthwaite (2010). The step size shrinks when a proposal is rejected and
expands when a proposal is accepted.

Gradient ascent is criticized for sometimes being relatively slow when
close to the maximum, and its asymptotic rate of convergence is inferior
to other methods. However, compared to other popular optimization
algorithms such as Newton-Raphson, an advantage of the gradient ascent
is that it works in infinite dimensions, requiring only sufficient
computer memory. Although Newton-Raphson converges in fewer iterations,
calculating the inverse of the negative Hessian matrix of
second-derivatives is more computationally expensive and subject to
singularities. Therefore, gradient ascent takes longer to converge, but
is more generalizable.

When `Method="BFGS"`, the BFGS algorithm is used, which was proposed by
Broyden (1970), Fletcher (1970), Goldfarb (1970), and Shanno (1970),
independently. BFGS may be the most efficient and popular quasi-Newton
optimization algorithm. As a quasi-Newton algorithm, the Hessian matrix
is approximated using rank-one updates specified by (approximate)
gradient evaluations. Since BFGS is very popular, there are many
variations of it. This is a version by Nash that has been adapted from
the Rvmmin package, and is used in the `optim` function of base R. The
approximate Hessian is not guaranteed to converge to the Hessian. When
BFGS is used, the approximate Hessian is not used to calculate the final
covariance matrix.

When `Method="BHHH"`, the algorithm of Berndt et al. (1974) is used,
which is commonly pronounced B-triple H. The BHHH algorithm is a
quasi-Newton method that includes a step-size parameter, partial
derivatives, and an approximation of a covariance matrix that is
calculated as the inverse of the sum of the outer product of the
gradient (OPG), calculated from each record. The OPG method becomes more
costly with data sets with more records. Since partial derivatives must
be calculated per record of data, the list of data has special
requirements with this method, and must include design matrix `X`, and
dependent variable `y` or `Y`. Records must be row-wise. An advantage of
BHHH over NR (see below) is that the covariance matrix is necessarily
positive definite, and guaranteed to provide an increase in LP each
iteration (given a small enough step-size), even in convex areas. The
covariance matrix is better approximated with larger data sample sizes,
and when closer to the maximum of LP. Disadvantages of BHHH include that
it can give small increases in LP, especially when far from the maximum
or when LP is highly non-quadratic.

When `Method="CG"`, a nonlinear conjugate gradient algorithm is used. CG
uses partial derivatives, but does not use the Hessian matrix or any
approximation of it. CG usually requires more iterations to reach
convergence than other algorithms that use the Hessian or an
approximation. However, since the Hessian becomes computationally
expensive as the dimension of the model grows, CG is applicable to large
dimensional models when `CovEst="Hessian"` is avoided. CG was originally
developed by Hestenes and Stiefel (1952), though this version is adapted
from the `Rcgminu` function in package Rcgmin.

When `Method="DFP"`, the Davidon-Fletcher-Powell algorithm is used. DFP
was the first popular, multidimensional, quasi-Newton optimization
algorithm. The DFP update of an approximate Hessian matrix maintains
symmetry and positive-definiteness. The approximate Hessian is not
guaranteed to converge to the Hessian. When DFP is used, the approximate
Hessian is not used to calculate the final covariance matrix. Although
DFP is very effective, it was superseded by the BFGS algorithm.

When `Method="HAR"`, a hit-and-run algorithm with a multivariate
proposal and adaptive length is used. The length parameter is adapted
each iteration with the univariate version of the Robbins-Monro
stochastic approximation in Garthwaite (2010). The length shrinks when a
proposal is rejected and expands when a proposal is accepted. This is
the same algorithm as the HARM or Hit-And-Run Metropolis MCMC algorithm
with adaptive length, except that a Metropolis step is not used.

When `Method="HJ"`, the Hooke-Jeeves (1961) algorithm is used. This was
adapted from the `HJK` algorithm in package dfoptim. Hooke-Jeeves is a
derivative-free, direct search method. Each iteration involves two
steps: an exploratory move and a pattern move. The exploratory move
explores local behavior, and the pattern move takes advantage of pattern
direction. It is sometimes described as a hill-climbing algorithm. If
the solution improves, it accepts the move, and otherwise rejects it.
Step size decreases with each iteration. The decreasing step size can
trap it in local maxima, where it gets stuck and converges erroneously.
Users are encouraged to attempt again after what seems to be
convergence, starting from the latest point. Although getting stuck at
local maxima can be problematic, the Hooke-Jeeves algorithm is also
attractive because it is simple, fast, does not depend on derivatives,
and is otherwise relatively robust.

When `Method="LBFGS"`, the limited-memory BFGS
(Broyden-Fletcher-Goldfarb-Shanno) algorithm is called in `optim`, once
per iteration.

When `Method="LM"`, the Levenberg-Marquardt algorithm (Levenberg, 1944;
Marquardt, 1963) is used. Also known as the Levenberg-Marquardt
Algorithm (LMA) or the Damped Least-Squares (DLS) method, LM is a trust
region (not to be confused with TR below) quasi-Newton optimization
algorithm that minimizes nonlinear least squares, and has been adapted
here to maximize LP. LM uses partial derivatives and approximates the
Hessian with outer-products. It is suitable for nonlinear optimization
up to a few hundred parameters, but loses its efficiency in larger
problems due to matrix inversion. LM is considered between the
Gauss-Newton algorithm and gradient descent. When far from the solution,
LM moves slowly like gradient descent, but is guaranteed to converge.
When LM is close to the solution, LM becomes a damped Gauss-Newton
method. This was adapted from the `lsqnonlin` algorithm in package
pracma.

When `Method="NM"`, the Nelder-Mead (1965) algorithm is used. This was
adapted from the `nelder_mead` function in package pracma. Nelder-Mead
is a derivative-free, direct search method that is known to become
inefficient in large-dimensional problems. As the dimension increases,
the search direction becomes increasingly orthogonal to the steepest
ascent (usually descent) direction. However, in smaller dimensions, it
is a popular algorithm. At each iteration, three steps are taken to
improve a simplex: reflection, extension, and contraction.

When `Method="NR"`, the Newton-Raphson optimization algorithm, also
known as Newton's Method, is used. Newton-Raphson uses derivatives and a
Hessian matrix. The algorithm is included for its historical
significance, but is known to be problematic when starting values are
far from the targets, and calculating and inverting the Hessian matrix
can be computationally expensive. As programmed here, when the Hessian
is problematic, it tries to use only the derivatives, and when that
fails, a jitter is applied. Newton-Raphson should not be the first
choice of the user, and BFGS should always be preferred.

When `Method="PSO"`, the Standard Particle Swarm Optimization 2007
algorithm is used. A swarm of particles is moved according to velocity,
neighborhood, and the best previous solution. The neighborhood for each
particle is a set of informing particles. PSO is derivative-free. PSO
has been adapted from the `psoptim` function in package pso.

When `Method="Rprop"`, the approximate gradient is taken for each
parameter in each iteration, and its sign is compared to the approximate
gradient in the previous iteration. A weight element in a weight vector
is associated with each approximate gradient. A weight element is
multiplied by 1.2 when the sign does not change, or by 0.5 if the sign
changes. The weight vector is the step size, and is constrained to the
interval \[0.001, 50\], and initial weights are 0.0125. This is the
resilient backpropagation algorithm, which is often denoted as the
"Rprop-" algorithm of Riedmiller (1994).

When `Method="SGD"`, a stochastic gradient descent algorithm is used
that is designed only for big data, and gained popularity after
successful use in the NetFlix competition. This algorithm has special
requirements for the `Model` specification function and the `Data` list.
See the "lucifer tutorial" vignette for more information.

When `Method="SOMA"`, a population of ten particles or individuals moves
in the direction of the best particle, the leader. The leader does not
move in each iteration, and a line-search is used for each non-leader,
up to three times the difference in parameter values between each
non-leader and leader. This algorithm is derivative-free and often
considered in the family of evolution algorithms. Numerous model
evaluations are performed per non-leader per iteration. This algorithm
was adapted from package soma.

When `Method="SPG"`, a Spectral Projected Gradient algorithm is used.
SPG is a non-monotone algorithm that is suitable for high-dimensional
models. The approximate gradient is used, but the Hessian matrix is not.
When used with large models, `CovEst="Hessian"` should be avoided. SPG
has been adapted from the `spg` function in package BB.

When `Method="SR1"`, the Symmetric Rank-One (SR1) algorithm is used. SR1
is a quasi-Newton algorithm, and the Hessian matrix is approximated,
often without being positive-definite. At the posterior modes, the true
Hessian is usually positive-definite, but this is often not the case
during optimization when the parameters have not yet reached the
posterior modes. Other restrictions, including constraints, often result
in the true Hessian being indefinite at the solution. For these reasons,
SR1 often outperforms BFGS. The approximate Hessian is not guaranteed to
converge to the Hessian. When SR1 is used, the approximate Hessian is
not used to calculate the final covariance matrix.

When `Method="TR"`, the Trust Region algorithm of Nocedal and Wright
(1999) is used. The TR algorithm attempts to reach its objective in the
fewest number of iterations, is therefore very efficient, as well as
safe. The efficiency of TR is attractive when model evaluations are
expensive. The Hessian is approximated each iteration, making TR best
suited to models with small to medium dimensions, say up to a few
hundred parameters. TR has been adapted from the `trust` function in
package trust.

## References

Azevedo-Filho, A. and Shachter, R. (1994). "Laplace's Method
Approximations for Probabilistic Inference in Belief Networks with
Continuous Variables". In "Uncertainty in Artificial Intelligence",
Mantaras, R. and Poole, D., Morgan Kauffman, San Francisco, CA, p.
28–36.

Bernardo, J.M. and Smith, A.F.M. (2000). "Bayesian Theory". John Wiley &
Sons: West Sussex, England.

Berndt, E., Hall, B., Hall, R., and Hausman, J. (1974). "Estimation and
Inference in Nonlinear Structural Models". *Annals of Economic and
Social Measurement*, 3, p. 653–665.

Broyden, C.G. (1970). "The Convergence of a Class of Double Rank
Minimization Algorithms: 2. The New Algorithm". Journal of the Institute
of Mathematics and its Applications, 6, p. 76–90.

Fletcher, R. (1970). "A New Approach to Variable Metric Algorithms".
Computer Journal, 13(3), p. 317–322.

Garthwaite, P., Fan, Y., and Sisson, S. (2010). "Adaptive Optimal
Scaling of Metropolis-Hastings Algorithms Using the Robbins-Monro
Process."

Goldfarb, D. (1970). "A Family of Variable Metric Methods Derived by
Variational Means". Mathematics of Computation, 24(109), p. 23–26.

Hestenes, M.R. and Stiefel, E. (1952). "Methods of Conjugate Gradients
for Solving Linear Systems". *Journal of Research of the National Bureau
of Standards*, 49(6), p. 409–436.

Hooke, R. and Jeeves, T.A. (1961). "'Direct Search' Solution of
Numerical and Statistical Problems". *Journal of the Association for
Computing Machinery*, 8(2), p. 212–229.

Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". *Journal of the
American Statistical Association*, 90(430), p. 773–795.

Laplace, P. (1774). "Memoire sur la Probabilite des Causes par les
Evenements". l'Academie Royale des Sciences, 6, 621–656. English
translation by S.M. Stigler in 1986 as "Memoir on the Probability of the
Causes of Events" in Statistical Science, 1(3), 359–378.

Laplace, P. (1814). "Essai Philosophique sur les Probabilites". English
translation in Truscott, F.W. and Emory, F.L. (2007) from (1902) as "A
Philosophical Essay on Probabilities". ISBN 1602063281, translated from
the French 6th ed. (1840).

Levenberg, K. (1944). "A Method for the Solution of Certain Non-Linear
Problems in Least Squares". *Quarterly of Applied Mathematics*, 2, p.
164–168.

Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
Posterior Simulation with the Laplace-Metropolis Estimator". *Journal of
the American Statistical Association*, 92, p. 648–655.

Marquardt, D. (1963). "An Algorithm for Least-Squares Estimation of
Nonlinear Parameters". *SIAM Journal on Applied Mathematics*, 11(2), p.
431–441.

Nelder, J.A. and Mead, R. (1965). "A Simplex Method for Function
Minimization". *The Computer Journal*, 7(4), p. 308–313.

Nocedal, J. and Wright, S.J. (1999). "Numerical Optimization".
Springer-Verlag.

Riedmiller, M. (1994). "Advanced Supervised Learning in Multi-Layer
Perceptrons - From Backpropagation to Adaptive Learning Algorithms".
*Computer Standards and Interfaces*, 16, p. 265–278.

Shanno, D.F. (1970). "Conditioning of quasi-Newton Methods for Function
Minimization". Mathematics of Computation, 24(111), p. 647–650.

Tierney, L. and Kadane, J.B. (1986). "Accurate Approximations for
Posterior Moments and Marginal Densities". *Journal of the American
Statistical Association*, 81(393), p. 82–86.

Tierney, L., Kass, R., and Kadane, J.B. (1989). "Fully Exponential
Laplace Approximations to Expectations and Variances of Nonpositive
Functions". *Journal of the American Statistical Association*, 84(407),
p. 710–716.

Zelinka, I. (2004). "SOMA - Self Organizing Migrating Algorithm". In:
Onwubolu G.C. and Babu, B.V., editors. "New Optimization Techniques in
Engineering". Springer: Berlin, Germany.

## See also

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`BayesianBootstrap`](https://robustecologies.github.io/lucifer/reference/BayesianBootstrap.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`optim`](https://rdrr.io/r/stats/optim.html),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

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
Fit <- LaplaceApproximation(Model, Initial.Values, Data = MyData,
     Iterations = 100, Method = "SPG", CPUs = 1)
print(Fit)
summary(Fit)
plot(Fit, MyData, PDF = FALSE)

# Fit$Covar is scaled (2.38^2/d) and submitted to lucifer as Covar.
# Fit$Summary[,1] is submitted to lucifer as Initial.Values.

Pred <- predict(Fit, Model, MyData, CPUs = 1)
summary(Pred, Discrep = "Chi-Square")
plot(Pred, Style = "Fitted")
} # }
```
