# Variational Bayes methods in lucifer

## Introduction to variational inference

Bayesian inference requires computing the posterior distribution
\\p(\theta \mid y) = p(y \mid \theta) \\ p(\theta) / p(y)\\, where
\\\theta\\ denotes the vector of model parameters and \\y\\ the observed
data. The normalizing constant \\p(y) = \int p(y \mid \theta) \\
p(\theta) \\ d\theta\\ is generally intractable for all but the simplest
conjugate models, and Markov chain Monte Carlo methods obtain samples
from the posterior without evaluating this integral. Variational Bayes
(VB) takes a fundamentally different approach: it recasts posterior
inference as an optimization problem by searching over a family of
tractable distributions \\\mathcal{Q} = \\q(\theta; \lambda) : \lambda
\in \Lambda\\\\ for the member closest to the true posterior
[\[1\]](#ref1).

### Kullback-Leibler divergence

The notion of “closeness” is formalized through the Kullback-Leibler
divergence from \\q\\ to \\p\\:

\\ \mathrm{KL}(q \\\\\\ p) = \int q(\theta; \lambda) \log
\frac{q(\theta; \lambda)}{p(\theta \mid y)} \\ d\theta =
\mathbb{E}\_q\\\left\[\log q(\theta; \lambda)\right\] -
\mathbb{E}\_q\\\left\[\log p(\theta \mid y)\right\]. \\

This quantity is non-negative and equals zero if and only if \\q = p\\
almost everywhere. Minimizing \\\mathrm{KL}(q \\\\\\ p)\\ directly is
infeasible because \\p(\theta \mid y)\\ contains the unknown normalizing
constant \\p(y)\\. However, substituting \\\log p(\theta \mid y) = \log
p(\theta, y) - \log p(y)\\ into the expression and rearranging yields

\\ \log p(y) = \underbrace{\mathbb{E}\_q\\\left\[\log p(\theta,
y)\right\] - \mathbb{E}\_q\\\left\[\log q(\theta;
\lambda)\right\]}\_{\mathcal{L}(\lambda)} + \mathrm{KL}(q \\\\\\ p). \\

Since \\\log p(y)\\ is a constant with respect to \\\lambda\\ and
\\\mathrm{KL}(q \\\\\\ p) \geq 0\\, the quantity
\\\mathcal{L}(\lambda)\\ is a lower bound on the log-evidence.
Maximizing this evidence lower bound (ELBO) is equivalent to minimizing
the KL divergence.

### The ELBO

The ELBO decomposes naturally into two interpretable terms:

\\ \mathcal{L}(\lambda) = \underbrace{\mathbb{E}\_q\\\left\[\log
p(\theta, y)\right\]}\_{\text{expected log-joint}} +
\underbrace{\mathbb{H}\[q\]}\_{\text{entropy of } q} \\

where \\\mathbb{H}\[q\] = -\mathbb{E}\_q\[\log q(\theta; \lambda)\]\\ is
the entropy of the variational distribution. The first term encourages
\\q\\ to place mass where the joint density \\p(\theta, y)\\ is large,
while the entropy term prevents \\q\\ from collapsing to a point mass.
This tension between fit and spread drives the entire variational
inference machinery.

### The mean-field assumption

The most common structural restriction on \\\mathcal{Q}\\ is the
mean-field assumption, which factorizes the variational distribution
over groups of parameters:

\\ q(\theta; \lambda) = \prod\_{j=1}^{J} q_j(\theta_j; \lambda_j). \\

This factorization yields tractable coordinate-ascent updates in
conjugate models but, importantly, ignores all posterior correlations
between groups. When the true posterior exhibits strong dependencies,
the mean-field approximation can severely underestimate marginal
variances and produce overconfident inference [\[3\]](#ref3). The eight
methods implemented in **lucifer** span the spectrum from fully
factorized (mean-field) to correlation-preserving (full-rank),
nonparametric (particle-based), and quasi-Newton trajectory-based
approximations.

## Fixed-form VB: Salimans1 and Salimans2

Salimans and Knowles [\[1\]](#ref1) introduced a fixed-form variational
Bayes framework that recasts the ELBO optimization as a stochastic
linear regression problem. The variational family is a multivariate
Gaussian \\q(\theta) = \mathcal{N}(\theta; m, V)\\, and the task is to
estimate the mean \\m\\ and covariance \\V\\ that maximize the ELBO. The
key insight is that the gradient of the ELBO with respect to the natural
parameters of this Gaussian can be expressed as the coefficients of a
linear regression of the log-joint density on the sufficient statistics
of \\q\\, evaluated at samples drawn from \\q\\ itself.

### The stochastic regression interpretation

At each iteration \\t\\, the algorithm draws a sample \\\theta^\* \sim
q_t(\theta) = \mathcal{N}(m_t, V_t)\\, evaluates the gradient \\g =
\nabla\_\theta \log p(\theta^\*, y)\\ and (for Salimans2) the Hessian
\\H = \nabla^2\_\theta \log p(\theta^\*, y)\\ at \\\theta^\*\\, and
performs stochastic averaging updates on the natural parameters.
Specifically, the precision matrix \\P = V^{-1}\\ and the
precision-weighted mean \\a = P \cdot m\\ are updated according to

\\ a\_{t+1} = (1 - w) \\ a_t + w \\ g, \qquad P\_{t+1} = (1 - w) \\
P_t - w \\ H \\

where the step size \\w = 1/\sqrt{T}\\ with \\T\\ the total number of
iterations. The mean is then recovered as \\m = V \cdot a + z\\, where
\\z\\ is a similarly averaged sample location. Convergence is assessed
in the second half of iterations using Polyak averaging of both \\m\\
and \\V\\.

### Gradient-only vs. gradient + Hessian

**Salimans2** (the default method in
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md))
uses both the gradient and the Hessian of the log-joint density at the
sampled point. Computing the Hessian requires \\O(J^2/2)\\ additional
model evaluations per iteration (via central cross-differences), which
yields an accurate curvature estimate but becomes expensive when the
number of parameters \\J\\ is large.

**Salimans1** avoids the Hessian entirely by approximating the curvature
with the negative outer product of the gradient, \\\hat{H} = -g \\
g^\top\\. This reduces the per-iteration cost to \\O(J + 1)\\ model
evaluations, making it practical for moderate-dimensional problems. The
outer-product approximation is a positive semidefinite estimate of the
negative Hessian, which provides a conservative (often
diagonal-dominant) update to the precision matrix. In practice,
Salimans1 converges more slowly than Salimans2 for smooth,
low-dimensional posteriors, but scales more gracefully when \\J\\
exceeds approximately 20 parameters.

Both methods share the same Polyak-averaged convergence criterion and
accept-reject step, so switching between them requires only changing the
`Method` argument.

## Automatic differentiation variational inference (ADVI)

The ADVI framework of Kucukelbir et al. [\[3\]](#ref3) automates
variational inference for a broad class of models by combining a
Gaussian variational family with the reparameterization trick and
adaptive stochastic optimization. **lucifer** implements both the
mean-field (`ADVI.mf`) and full-rank (`ADVI.fr`) variants.

### The reparameterization trick

The central difficulty in optimizing the ELBO is that the expectation
\\\mathbb{E}\_q\[\cdot\]\\ is taken with respect to a distribution that
itself depends on \\\lambda\\. The reparameterization trick resolves
this by expressing a sample \\\theta \sim q(\theta; \lambda)\\ as a
deterministic function of a fixed noise variable \\\varepsilon \sim
\mathcal{N}(0, I)\\:

\\ \theta = \mu + \sigma \odot \varepsilon \quad \text{(mean-field)},
\qquad \theta = \mu + L \varepsilon \quad \text{(full-rank)} \\

where \\\mu\\ is the variational mean, \\\sigma\\ is the vector of
standard deviations (mean-field), \\L\\ is a lower-triangular Cholesky
factor (full-rank), and \\\odot\\ denotes element-wise multiplication.
This transformation moves the dependence on \\\lambda = (\mu, \sigma)\\
or \\\lambda = (\mu, L)\\ outside the expectation, enabling
straightforward gradient estimation by differentiating through the
sample.

### The softplus parameterization

Since standard deviations and Cholesky diagonal elements must be
positive, ADVI applies the softplus transformation \\\sigma_j = \log(1 +
e^{\omega_j})\\ to unconstrained parameters \\\omega_j\\. The softplus
function is a smooth approximation to the rectifier that avoids the
gradient discontinuity at zero. Its inverse, \\\omega = \log(e^\sigma -
1)\\, is used for initialization. The implementation in **lucifer** uses
a numerically stable branch: for \\\omega \> 20\\,
\\\text{softplus}(\omega) \approx \omega\\ avoids overflow in the
exponential.

### Adam optimizer

Both ADVI variants optimize the variational parameters using the Adam
algorithm [\[7\]](#ref7), an adaptive first-order method that maintains
per-parameter exponential moving averages of the gradient (\\m_t\\) and
squared gradient (\\v_t\\):

\\ m_t = \beta_1 \\ m\_{t-1} + (1 - \beta_1) \\ g_t, \qquad v_t =
\beta_2 \\ v\_{t-1} + (1 - \beta_2) \\ g_t^2 \\

with bias-corrected estimates \\\hat{m}\_t = m_t / (1 - \beta_1^t)\\ and
\\\hat{v}\_t = v_t / (1 - \beta_2^t)\\. The parameter update is
\\\lambda\_{t+1} = \lambda_t + \alpha \\ \hat{m}\_t /
(\sqrt{\hat{v}\_t} + \epsilon)\\ with default learning rate \\\alpha =
0.01\\, \\\beta_1 = 0.9\\, \\\beta_2 = 0.999\\, and \\\epsilon =
10^{-8}\\. Adam’s per-coordinate adaptive scaling is well suited to
variational problems where gradients differ substantially in magnitude
across parameters.

### Mean-field vs. full-rank tradeoffs

**ADVI.mf** approximates the posterior with a fully factorized Gaussian,
\\q(\theta) = \prod_j \mathcal{N}(\theta_j; \mu_j, \sigma_j^2)\\,
requiring only \\2J\\ variational parameters. This is the cheapest ADVI
option and works well when posterior correlations are weak. However, it
systematically underestimates marginal variances of correlated
parameters.

**ADVI.fr** replaces the diagonal covariance with \\\Sigma = L L^\top\\
where \\L\\ is lower-triangular, capturing the full posterior
correlation structure at the cost of \\J + J(J-1)/2\\ additional
variational parameters for the Cholesky factor. This quadratic scaling
means that ADVI.fr becomes expensive for \\J \> 50\\; **lucifer** issues
a warning when this threshold is exceeded. ADVI.fr converges more slowly
per iteration than ADVI.mf (because there are more parameters to
optimize) but can produce substantially better posterior approximations
when correlations are strong.

Convergence for both ADVI variants is assessed by the relative change in
the ELBO between successive iterations, and the algorithm terminates
when this change falls below `Stop.Tolerance`.

## Black-box variational inference (BBVI)

Black-box variational inference, introduced by Ranganath, Gerrish, and
Blei [\[4\]](#ref4), provides a fully generic framework that requires
only the ability to sample from the variational distribution and
evaluate its log-density. Unlike ADVI, BBVI does not require the
reparameterization trick and can therefore be applied to discrete or
non-differentiable variational families, although **lucifer** uses a
Gaussian family for consistency with the other methods.

### The score function estimator

BBVI uses the score function estimator (also called the REINFORCE
estimator) to compute gradients of the ELBO. Starting from the identity

\\ \nabla\_\lambda \mathcal{L}(\lambda) =
\mathbb{E}\_q\\\left\[\left(\log p(\theta, y) - \log q(\theta;
\lambda)\right) \nabla\_\lambda \log q(\theta; \lambda)\right\], \\

the gradient is approximated by drawing \\S\\ samples \\\theta_1,
\ldots, \theta_S \sim q(\theta; \lambda)\\ and computing

\\ \widehat{\nabla}\_\lambda \mathcal{L} = \frac{1}{S} \sum\_{s=1}^{S}
f_s \\ \nabla\_\lambda \log q(\theta_s; \lambda), \qquad f_s = \log
p(\theta_s, y) - \log q(\theta_s; \lambda). \\

The quantity \\f_s\\ is the difference between the log-joint and the
log-variational density, and \\\nabla\_\lambda \log q(\theta_s;
\lambda)\\ is the score function of the variational family. The
**lucifer** implementation uses \\S = 10\\ Monte Carlo samples per
iteration by default.

### Control variates and Rao-Blackwellization

The score function estimator has notoriously high variance, which can
make optimization unstable. BBVI mitigates this through two variance
reduction techniques. First, per-component control variates (baselines)
are subtracted from \\f_s\\. For each variational parameter
\\\lambda_j\\, the optimal baseline is

\\ b_j = \frac{\text{Cov}(f \cdot \nabla_j \log q, \\ \nabla_j \log
q)}{\text{Var}(\nabla_j \log q)} \\

computed from the current batch of samples. This is a form of
Rao-Blackwellization: by exploiting the structure of the score function,
each component of the gradient uses only the terms that are relevant to
that particular variational parameter, reducing the effective noise. The
**lucifer** implementation computes these baselines independently for
the mean parameters \\\mu\\ and the log-scale parameters \\\log
\sigma\\, using the sample covariance and variance within each batch.

Even with control variates, BBVI typically requires more iterations than
reparameterization-based methods (ADVI) to achieve comparable ELBO
values. Its advantage lies in generality: it is the only method in
**lucifer** that does not require differentiating through the model,
making it a fallback when the log-joint density is not smooth.

## Stein variational gradient descent (SVGD)

Stein variational gradient descent, introduced by Liu and Wang
[\[5\]](#ref5), is a nonparametric particle-based method that occupies a
unique niche among the variational algorithms in **lucifer**. Rather
than parameterizing a single distribution, SVGD maintains a set of \\K\\
particles \\\\\theta_k\\\_{k=1}^K\\ and iteratively transports them
toward the posterior by following the direction of steepest descent in
the KL divergence within a reproducing kernel Hilbert space (RKHS).

### Stein’s identity and the RKHS

The foundation of SVGD is Stein’s identity: for a smooth density
\\p(\theta)\\ and a smooth vector-valued function \\\phi(\theta)\\,

\\ \mathbb{E}\_p\\\left\[\mathcal{A}\_p \phi(\theta)\right\] = 0, \qquad
\mathcal{A}\_p \phi(\theta) = \nabla\_\theta \log p(\theta) \cdot
\phi(\theta) + \nabla\_\theta \cdot \phi(\theta) \\

where \\\mathcal{A}\_p\\ is the Stein operator. Liu and Wang showed that
the direction of steepest descent of \\\mathrm{KL}(q \\\\\\ p)\\ in the
unit ball of an RKHS \\\mathcal{H}\\ with kernel \\k(\cdot, \cdot)\\ is

\\ \phi^\*(\theta) = \mathbb{E}\_{q}\\\left\[k(\theta', \theta) \\
\nabla\_{\theta'} \log p(\theta') + \nabla\_{\theta'} k(\theta',
\theta)\right\]. \\

The first term drives each particle toward regions of high posterior
density (the “attractive” force), while the second term is a repulsive
force that prevents particles from collapsing to a single point, thereby
maintaining diversity in the approximation.

### Particle dynamics

In practice, the expectation over \\q\\ is replaced by the empirical
average over the \\K\\ particles, yielding the update

\\ \theta_k \leftarrow \theta_k + \epsilon \\ \hat{\phi}^\*(\theta_k),
\qquad \hat{\phi}^\*(\theta_k) = \frac{1}{K} \sum\_{j=1}^{K}
\left\[k(\theta_j, \theta_k) \\ \nabla\_{\theta_j} \log p(\theta_j, y) +
\nabla\_{\theta_j} k(\theta_j, \theta_k)\right\] \\

where \\\epsilon\\ is the step size. In **lucifer**, each particle is
updated using Adam (with per-particle state) rather than a fixed step
size, which improves convergence stability.

### The median heuristic bandwidth

SVGD uses a radial basis function (RBF) kernel \\k(\theta, \theta') =
\exp\\\left(-\\\theta - \theta'\\^2 / (2h^2)\right)\\, and the bandwidth
\\h\\ is set adaptively at each iteration via the median heuristic:

\\ h = \frac{\text{median}(\\\theta_i - \theta_j\\)}{\sqrt{2 \log(K +
1)}}. \\

This rule adapts the kernel width to the current spread of the
particles. When particles are far apart, \\h\\ is large and the kernel
interactions are long-range; as particles converge toward the posterior,
\\h\\ shrinks and the interactions become local. The \\\sqrt{2 \log(K +
1)}\\ normalization ensures that the effective bandwidth scales
appropriately with the number of particles.

The default number of particles is \\K = 20\\, which can be changed by
passing an integer to the `Covar` argument. SVGD requires \\K \times
(J + 1)\\ model evaluations per iteration (one evaluation plus \\J\\
finite-difference gradient evaluations per particle), making it the most
expensive method per iteration. However, it is the only method in
**lucifer** capable of representing multimodal posteriors, since the
particles can populate distinct modes simultaneously.

## Natural gradient variational inference (NGD)

Natural gradient methods, originating with Amari’s work on information
geometry [\[6\]](#ref6) and developed for variational inference by Khan
and Lin [\[2\]](#ref2), replace the ordinary (Euclidean) gradient with
the natural gradient, which accounts for the curvature of the
variational family’s parameter space.

### The Fisher information metric

The space of probability distributions has a Riemannian geometry induced
by the Fisher information matrix

\\ F(\lambda) = \mathbb{E}\_q\\\left\[\nabla\_\lambda \log q(\theta;
\lambda) \\ \nabla\_\lambda \log q(\theta; \lambda)^\top\right\]. \\

The natural gradient is defined as \\\tilde{\nabla}\_\lambda \mathcal{L}
= F(\lambda)^{-1} \\ \nabla\_\lambda \mathcal{L}\\, and it gives the
steepest ascent direction of \\\mathcal{L}\\ when distance in parameter
space is measured by the KL divergence rather than Euclidean distance.
This has the profound consequence that the natural gradient is invariant
to reparameterization of the variational family: different
parameterizations of the same Gaussian family yield identical updates.

### Natural parameters

For the Gaussian variational family \\q(\theta) = \mathcal{N}(\mu,
\Sigma)\\, the natural parameters are \\\eta_1 = \Sigma^{-1} \mu\\ and
\\\eta_2 = -\frac{1}{2} \Sigma^{-1}\\. A remarkable property of
exponential families is that the natural gradient update in natural
parameters takes an especially simple form. The **lucifer**
implementation updates natural parameters according to

\\ \eta_1^{(t+1)} = (1 - \rho_t) \\ \eta_1^{(t)} + \rho_t \\ g, \qquad
\eta_2^{(t+1)} = (1 - \rho_t) \\ \eta_2^{(t)} - \tfrac{1}{2} \rho_t \\ H
\\

where \\g\\ and \\H\\ are the gradient and Hessian of \\\log p(\theta,
y)\\ evaluated at a sample from \\q\\. The step size follows a
Robbins-Monro schedule \\\rho_t = (t + \tau)^{-\alpha}\\ with \\\tau =
10\\ and \\\alpha = 0.75\\, which satisfies the conditions \\\sum_t
\rho_t = \infty\\ and \\\sum_t \rho_t^2 \< \infty\\ necessary for
convergence of stochastic approximation.

### Relationship to mirror descent and Salimans2

The NGD update in natural parameters is structurally very similar to
Salimans2: both perform stochastic averaging of the gradient and
Hessian, and both maintain a full precision matrix. The key difference
is the step-size schedule. Salimans2 uses a fixed \\w = 1/\sqrt{T}\\
that requires knowing the total number of iterations in advance, while
NGD uses the Robbins-Monro schedule that decays adaptively. From an
optimization-theoretic perspective, the NGD update in natural parameters
is equivalent to mirror descent with the Bregman divergence induced by
the log-partition function of the Gaussian family [\[2\]](#ref2). This
connection provides convergence guarantees under mild regularity
conditions that Salimans2 does not formally enjoy.

In practice, NGD tends to converge faster than Salimans2 in the initial
phase (due to the larger early step sizes) but may require more
fine-tuning of \\\tau\\ and \\\alpha\\ for difficult posteriors. For
small \\J\\, NGD and Salimans2 produce nearly identical results.

## Pathfinder: quasi-Newton variational inference

Pathfinder, introduced by Zhang, Carpenter, Gelman, and Vehtari
[\[10\]](#ref10), takes a fundamentally different approach to
variational inference. Rather than iteratively optimizing variational
parameters to maximize the ELBO, Pathfinder exploits the trajectory of
an L-BFGS optimizer to construct a sequence of Gaussian approximations
and selects the best one. This yields an algorithm that typically
requires 1-2 orders of magnitude fewer gradient evaluations than ADVI
while producing comparable posterior approximation quality.

### The L-BFGS trajectory

The algorithm begins by running an L-BFGS optimizer to maximize the
log-posterior density \\\log p(\theta \mid y)\\. L-BFGS is a
limited-memory quasi-Newton method that maintains a compact
representation of the inverse Hessian using the \\m\\ most recent pairs
of position and gradient differences \\(s_i, y_i)\\ where \\s_i =
\theta\_{i+1} - \theta_i\\ and \\y_i = \nabla f\_{i+1} - \nabla f_i\\
with \\f = -\log p(\theta \mid y)\\. At each iterate \\k\\, the L-BFGS
inverse Hessian approximation \\H_k\\ is reconstructed from these stored
pairs via the compact outer-product form of Nocedal and Wright (2006),
beginning with the initial scaling \\H_k^0 = \gamma_k I\\ where
\\\gamma_k = s\_{k-1}^\top y\_{k-1} / (y\_{k-1}^\top y\_{k-1})\\.

The key insight is that each iterate \\\theta_k\\ together with its
associated inverse Hessian \\H_k\\ defines a Gaussian approximation to
the posterior:

\\ q_k(\theta) = \mathcal{N}(\theta \mid \theta_k, H_k). \\

This approximation is available at no additional cost beyond what the
optimizer already computes.

### ELBO selection

Not all iterates along the optimization path yield equally good
approximations. Early iterates may be far from the mode, while later
iterates may have overly narrow covariance due to the curvature at the
optimum. Pathfinder evaluates the ELBO at each iterate:

\\ \mathcal{L}\_k = \frac{1}{S} \sum\_{s=1}^{S} \log p(\theta_k^{(s)},
y) + \frac{1}{2}\log\det(H_k) + \frac{J}{2}(1 + \log 2\pi) \\

where \\\theta_k^{(s)} \sim \mathcal{N}(\theta_k, H_k)\\ are Monte Carlo
draws, and selects \\k^\* = \arg\max_k \mathcal{L}\_k\\. The first term
estimates the expected log-joint under the approximation, while the
remaining terms give the Gaussian entropy. In **lucifer**, \\S = 10\\
draws are used for ELBO estimation (increased to 25 when \\J \> 50\\).

### Multi-path Pathfinder

A single L-BFGS run may get trapped in a local mode or may not
adequately explore the parameter space. Multi-path Pathfinder addresses
this by running \\M\\ independent single-path instances from jittered
initial values and combining their results. The procedure pools all
draws from the \\M\\ best approximations into a mixture:

\\ q(\theta) = \frac{1}{M} \sum\_{m=1}^{M} \mathcal{N}(\theta \mid
\theta\_{m}^\*, H\_{m}^\*) \\

and computes importance weights \\w_i = \log p(\theta_i, y) - \log
q(\theta_i)\\ at each pooled draw. Pareto smoothed importance sampling
(PSIS) is then applied to smooth extreme weights and produce a reliable
resampled set of posterior draws. The Pareto \\\hat{k}\\ diagnostic from
the PSIS step quantifies the quality of the multi-path approximation:

- \\\hat{k} \< 0.5\\: good; the importance weights are well-behaved
- \\0.5 \leq \hat{k} \< 0.7\\: acceptable; some weight variability but
  generally reliable
- \\0.7 \leq \hat{k} \< 1.0\\: questionable; the approximation may be
  inaccurate
- \\\hat{k} \geq 1.0\\: unreliable; the approximation is poor and MCMC
  should be used

In **lucifer**, the default is \\M = 4\\ paths, controlled by the
`Covar` argument. Setting `Covar = 1` gives single-path Pathfinder
without PSIS.

### When to use Pathfinder

Pathfinder is the recommended starting point for most variational
inference problems. It captures posterior correlations via the L-BFGS
inverse Hessian (unlike mean-field methods), runs faster than ADVI due
to the small number of gradient evaluations required, and provides a
built-in reliability diagnostic through the Pareto \\\hat{k}\\. Its main
limitation is the Gaussian approximation: for strongly non-Gaussian
posteriors (heavy tails, skewness, multimodality), the Pareto
\\\hat{k}\\ will be large, signaling that MCMC methods should be used
instead.

## Worked examples

The following two examples illustrate the practical differences between
the eight variational methods. The first fits a standard Bayesian linear
regression where all methods are expected to succeed, comparing
convergence speed, posterior density recovery, and computation time. The
second constructs a deliberately challenging posterior with strong
parameter correlations, exposing the limitations of mean-field
approximations and demonstrating why full-rank or trajectory-based
methods matter when posterior dependencies are strong.

### Example 1: Bayesian linear regression

#### Model specification

All variational methods in **lucifer** share the same interface: the
user specifies a model function, initial values, and a data list, then
calls
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
with the desired `Method`. The following example uses the `demonsnacks`
data set to fit a Bayesian linear regression of log-calories on
log-sodium, with weakly informative priors on the coefficients and a
half-Cauchy prior on the residual standard deviation.

``` r

library(lucifer)

# Load and prepare data
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[, 10] + 1)))
J <- ncol(X)
for (j in 2:J) X[, j] <- CenterScale(X[, j])

# Parameter setup
mon.names  <- "mu[1]"
parm.names <- as.parm.names(list(beta = rep(0, J), sigma = 0))
pos.beta   <- grep("beta", parm.names)
pos.sigma  <- grep("sigma", parm.names)

PGF <- function(Data) {
  beta  <- rnorm(Data$J)
  sigma <- runif(1)
  return(c(beta, sigma))
}

MyData <- list(
  J = J, PGF = PGF, X = X, mon.names = mon.names,
  parm.names = parm.names, pos.beta = pos.beta,
  pos.sigma = pos.sigma, y = y
)

Model <- function(parm, Data) {
  beta  <- parm[Data$pos.beta]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  beta.prior  <- sum(dnormv(beta, 0, 1000, log = TRUE))
  sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + beta.prior + sigma.prior
  list(LP = LP, Dev = -2 * LL, Monitor = mu[1],
       yhat = rnorm(length(mu), mu, sigma), parm = parm)
}

IV <- c(rep(0, J), 1)
```

The model function returns a list with five named elements: `LP`
(log-posterior), `Dev` (deviance, \\-2 \times\\ log-likelihood),
`Monitor` (monitored quantities), `yhat` (posterior predictive draws),
and `parm` (parameter vector, potentially constrained by
[`interval()`](https://robustecologies.github.io/lucifer/reference/interval.md)).
This interface is shared with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
and the other estimation functions.

#### Fitting all eight methods

The following loop fits every variational method using the same initial
values and tolerance. The `Covar` argument is repurposed for methods
that need it: SVGD uses it as the particle count and Pathfinder as the
number of independent L-BFGS paths. Each result is wrapped in `tryCatch`
so that a failure in one method does not halt the comparison.

``` r

methods <- c("Salimans1", "Salimans2", "ADVI.mf", "ADVI.fr",
             "BBVI", "SVGD", "NGD", "Pathfinder")
results <- list()

for (m in methods) {
  cat("\n===== Fitting:", m, "=====\n")
  covar_arg <- if (m == "SVGD") 20L else if (m == "Pathfinder") 4L else NULL
  iter_arg  <- if (m == "SVGD") 500L else 2000L
  results[[m]] <- tryCatch(
    VariationalBayes(
      Model, IV, Data = MyData, Covar = covar_arg,
      Iterations = iter_arg, Method = m,
      Stop.Tolerance = 1e-5, sir = TRUE, Samples = 1000
    ),
    error = function(e) {
      cat("  Method", m, "failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}

ok_methods <- Filter(function(m) !is.null(results[[m]]), methods)
```

#### Point estimates and timing

A first comparison extracts the posterior means, standard deviations,
iteration count, and elapsed time from each converged fit. On this
well-identified regression all methods should recover similar point
estimates; differences appear primarily in estimated variances and
computation cost.

``` r

pnames <- rownames(results[[ok_methods[1]]]$Summary1)

# OLS reference values
ols_fit  <- lm(y ~ X[, -1, drop = FALSE])
ols_vals <- setNames(c(coef(ols_fit), summary(ols_fit)$sigma), pnames)

comp <- data.frame(Method = ok_methods, stringsAsFactors = FALSE)
for (p in pnames) {
  comp[[paste0(p, " Mean")]] <- sapply(ok_methods, function(m)
    results[[m]]$Summary1[p, "Mean"])
  comp[[paste0(p, " SD")]] <- sapply(ok_methods, function(m)
    results[[m]]$Summary1[p, "SD"])
}
comp[["Iterations"]] <- sapply(ok_methods, function(m) results[[m]]$Iterations)
comp[["Minutes"]]    <- sapply(ok_methods, function(m) results[[m]]$Minutes)

# Insert OLS reference as first row
ref_row <- data.frame(Method = "OLS (reference)", stringsAsFactors = FALSE)
for (p in pnames) {
  ref_row[[paste0(p, " Mean")]] <- ols_vals[p]
  ref_row[[paste0(p, " SD")]] <- NA
}
ref_row[["Iterations"]] <- NA
ref_row[["Minutes"]]    <- NA
comp <- rbind(ref_row, comp)

knitr::kable(comp, digits = 3, align = c("l", rep("r", ncol(comp) - 1)),
  caption = "Posterior summaries and computation time across VB methods.")
```

#### Posterior density comparison

The most informative comparison overlays the marginal posterior
densities from each method. Methods that converge to the same solution
should produce nearly identical marginals; visible discrepancies
indicate either incomplete convergence or structural limitations of the
variational family.

``` r

library(ggplot2)

# Collect posterior draws from all methods (Gaussian approximation when
# SIR posterior is unavailable due to non-convergence)
post_list <- lapply(ok_methods, function(m) {
  df <- as.data.frame(.vb_draws(results[[m]]))
  df$Method <- m
  df
})
post_df <- do.call(rbind, post_list)
pnames_post <- setdiff(names(post_df), "Method")

post_long <- do.call(rbind, lapply(pnames_post, function(p) {
  data.frame(
    Parameter = p,
    Value     = post_df[[p]],
    Method    = post_df$Method,
    stringsAsFactors = FALSE
  )
}))

method_cols <- setNames(contrasting[seq_along(ok_methods)], ok_methods)

ggplot(post_long, aes(x = Value, color = Method)) +
  geom_density(linewidth = 0.6) +
  facet_wrap(~ Parameter, scales = "free") +
  scale_color_manual(values = method_cols) +
  labs(x = "Parameter value", y = "Density",
       title = "Marginal posterior densities across VB methods") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
```

#### Point estimates with credible intervals

A complementary view displays each method’s posterior mean alongside the
approximate 95% credible interval derived from the point-estimate
summary (mean \\\pm\\ 2 SD). Intervals that are notably narrower than
the consensus indicate variance underestimation, a hallmark of
mean-field methods on models with moderate posterior correlation.

``` r

est_df <- do.call(rbind, lapply(ok_methods, function(m) {
  s1 <- results[[m]]$Summary1
  data.frame(
    Method    = m,
    Parameter = rownames(s1),
    Mean      = s1[, "Mean"],
    LB        = s1[, "LB"],
    UB        = s1[, "UB"],
    stringsAsFactors = FALSE
  )
}))

ggplot(est_df, aes(x = Method, y = Mean, ymin = LB, ymax = UB, color = Method)) +
  geom_pointrange(size = 0.4) +
  facet_wrap(~ Parameter, scales = "free_y") +
  coord_flip() +
  scale_color_manual(values = method_cols) +
  labs(y = "Estimate (mean and 95% interval)",
       title = "Point estimates and credible intervals across VB methods") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
```

#### Convergence and computation time

The deviance history reveals how quickly each method finds the posterior
mode. The accompanying timing chart puts the convergence speed in
context: a method that converges in fewer iterations may still be slower
overall if each iteration is expensive (as with SVGD and its
per-particle gradient evaluations, or Salimans2 and NGD with their
Hessian computations).

``` r

dev_data <- do.call(rbind, lapply(ok_methods, function(m) {
  dev <- results[[m]]$Deviance
  data.frame(
    Iteration = seq_along(dev),
    Deviance  = dev,
    Method    = m,
    stringsAsFactors = FALSE
  )
}))

p_conv <- ggplot(dev_data, aes(x = Iteration, y = Deviance, color = Method)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = method_cols) +
  labs(x = "Iteration", y = "Deviance",
       title = "Deviance history") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

time_df <- data.frame(
  Method  = ok_methods,
  Minutes = sapply(ok_methods, function(m) results[[m]]$Minutes),
  stringsAsFactors = FALSE
)

p_time <- ggplot(time_df,
    aes(x = reorder(Method, Minutes), y = Minutes, fill = Method)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = method_cols) +
  labs(x = "", y = "Minutes",
       title = "Computation time") +
  theme_minimal(base_size = 11)

if (requireNamespace("gridExtra", quietly = TRUE))
  gridExtra::grid.arrange(p_conv, p_time, ncol = 2, widths = c(3, 2))
```

In typical runs, Salimans2 and NGD converge most quickly (100-300
iterations), followed by ADVI.mf, Salimans1, and Pathfinder. ADVI.fr
converges more slowly due to the larger parameter count. BBVI exhibits
the noisiest trace because of the variance inherent in the score
function estimator. SVGD takes fewer iterations but each one is
expensive. Because this is a well-identified regression, all methods
should recover similar posterior means; the differences surface in
estimated variances and computation time.

### Example 2: strongly correlated posterior

#### Motivation

The previous example is deliberately easy: the posterior is nearly
Gaussian and parameter correlations are moderate. Many real models
produce posteriors where parameters are tightly coupled, and this second
example constructs precisely that scenario. A regression with collinear
predictors creates a posterior for the regression coefficients that lies
along a narrow ridge rather than filling a roughly spherical region.
Mean-field methods approximate this ridge with an axis-aligned
ellipsoid, severely underestimating marginal variances and missing the
dominant correlation structure entirely. Full-rank and trajectory-based
methods should recover the ridge geometry faithfully.

#### Model specification

The model generates data from a linear regression with two predictors
correlated at \\\rho = 0.95\\ and known error variance \\\sigma = 0.5\\.
Using a known \\\sigma\\ makes the posterior for the three regression
coefficients exactly multivariate normal with an analytical closed form,
providing a ground truth against which every VB method can be
benchmarked without ambiguity.

``` r

set.seed(321)
n2 <- 200
rho_x <- 0.95
x1 <- rnorm(n2)
x2 <- rho_x * x1 + sqrt(1 - rho_x^2) * rnorm(n2)
X2 <- cbind(1, x1, x2)
true_beta2 <- c(1.0, 2.0, -1.5)
sigma_known <- 0.5
y2 <- as.numeric(X2 %*% true_beta2) + rnorm(n2, 0, sigma_known)

# Analytical posterior: beta | y ~ N(beta_hat, sigma^2 * (X'X)^{-1})
XtX <- crossprod(X2)
beta_hat  <- as.numeric(solve(XtX, crossprod(X2, y2)))
post_cov2 <- sigma_known^2 * solve(XtX)

# Model function (known sigma, flat prior)
J2 <- 3L
parm.names2 <- c("intercept", "beta1", "beta2")
mon.names2  <- "LP"

MyData2 <- list(
  parm.names = parm.names2,
  mon.names  = mon.names2,
  X = X2, y = y2, n = n2,
  sigma = sigma_known, J = J2,
  PGF = function(Data) rnorm(Data$J)
)

Model2 <- function(parm, Data) {
  beta <- parm
  mu   <- Data$X %*% beta
  LL   <- sum(dnorm(Data$y, mu, Data$sigma, log = TRUE))
  LP   <- LL  # flat prior
  list(LP = LP, Dev = -2 * LL, Monitor = LP,
       yhat = rnorm(Data$n, mu, Data$sigma), parm = parm)
}

IV2 <- rep(0, J2)
```

With \\\rho_x = 0.95\\ the posterior correlation between `beta1` and
`beta2` is approximately \\-0.95\\: a strongly tilted, elongated ellipse
that exposes the structural limitations of diagonal variational
families.

#### Fitting all methods

``` r

results2 <- list()

for (m in methods) {
  cat("\n===== Fitting:", m, "=====\n")
  covar_arg <- if (m == "SVGD") 20L else if (m == "Pathfinder") 4L else NULL
  iter_arg  <- if (m == "SVGD") 500L else 2000L
  results2[[m]] <- tryCatch(
    VariationalBayes(
      Model2, IV2, Data = MyData2, Covar = covar_arg,
      Iterations = iter_arg, Method = m,
      Stop.Tolerance = 1e-5, sir = TRUE, Samples = 1000
    ),
    error = function(e) {
      cat("  Method", m, "failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}

ok2 <- Filter(function(m) !is.null(results2[[m]]), methods)
```

#### Joint posterior: 2D scatter vs. analytical contour

The critical diagnostic for this problem is the joint distribution of
`beta1` and `beta2`. Each method’s posterior draws are displayed as a
scatter cloud, with the analytical 95% highest-density ellipse overlaid
in red. Methods that capture the correlation produce elongated clouds
aligned with the ellipse; mean-field methods produce nearly circular
clouds that spill far beyond the ellipse along one axis and fall short
along the other.

``` r

library(ggplot2)

# Analytical 95% HDR ellipse for (beta1, beta2)
sub_mean <- beta_hat[2:3]
sub_cov  <- post_cov2[2:3, 2:3]
eig      <- eigen(sub_cov)
angles   <- seq(0, 2 * pi, length.out = 200)
r95      <- sqrt(qchisq(0.95, df = 2))
ell      <- cbind(cos(angles), sin(angles)) %*%
            diag(sqrt(eig$values) * r95) %*% t(eig$vectors)
ellipse_df <- data.frame(
  beta1 = ell[, 1] + sub_mean[1],
  beta2 = ell[, 2] + sub_mean[2]
)

# Collect 2D posterior draws from all methods
scatter_list <- lapply(ok2, function(m) {
  draws <- .vb_draws(results2[[m]])
  data.frame(
    beta1  = draws[, "beta1"],
    beta2  = draws[, "beta2"],
    Method = m,
    stringsAsFactors = FALSE
  )
})
scatter_df <- do.call(rbind, scatter_list)

method_cols2 <- setNames(contrasting[seq_along(ok2)], ok2)

ggplot(scatter_df, aes(x = beta1, y = beta2, color = Method)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_path(data = ellipse_df, aes(x = beta1, y = beta2),
            color = "red", linewidth = 0.8, linetype = "dashed",
            inherit.aes = FALSE) +
  geom_point(data = data.frame(beta1 = sub_mean[1], beta2 = sub_mean[2]),
             aes(x = beta1, y = beta2),
             color = "red", shape = 3, size = 3, stroke = 1.2,
             inherit.aes = FALSE) +
  facet_wrap(~ Method, ncol = 4) +
  scale_color_manual(values = method_cols2) +
  labs(x = expression(beta[1]), y = expression(beta[2]),
       title = "Joint posterior draws vs. analytical 95% ellipse (red dashed)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
```

On this problem the contrast is stark. ADVI.mf and BBVI, both
mean-field, produce nearly circular scatter clouds centered at the
correct mode but with the wrong shape: they underestimate the marginal
variance for each coefficient because the strong negative correlation
between `beta1` and `beta2` is ignored. The full-rank methods
(Salimans2, NGD, ADVI.fr, Pathfinder) produce tilted ellipses that align
closely with the analytical contour. SVGD’s particle cloud typically
follows the ridge but may show more ragged edges due to the finite
particle count.

#### Marginal density comparison

The marginal densities for `beta1` and `beta2` reveal the same pattern
from a different angle. The analytical marginal
\\\mathcal{N}(\hat\beta_j, \[\sigma^2 (X^\top X)^{-1}\]\_{jj})\\ is
plotted as a black dashed curve. Mean-field methods produce marginals
that are too narrow because they attribute no variance to the
correlation component.

``` r

# Analytical marginal densities on a fine grid
grid_b1 <- seq(sub_mean[1] - 4 * sqrt(sub_cov[1, 1]),
               sub_mean[1] + 4 * sqrt(sub_cov[1, 1]), length.out = 300)
grid_b2 <- seq(sub_mean[2] - 4 * sqrt(sub_cov[2, 2]),
               sub_mean[2] + 4 * sqrt(sub_cov[2, 2]), length.out = 300)

analytic_df <- rbind(
  data.frame(Parameter = "beta1", Value = grid_b1,
    Density = dnorm(grid_b1, sub_mean[1], sqrt(sub_cov[1, 1]))),
  data.frame(Parameter = "beta2", Value = grid_b2,
    Density = dnorm(grid_b2, sub_mean[2], sqrt(sub_cov[2, 2])))
)

# VB posterior draws in long format (all methods)
marg_long <- do.call(rbind, lapply(c("beta1", "beta2"), function(p) {
  do.call(rbind, lapply(ok2, function(m) {
    draws <- .vb_draws(results2[[m]])
    data.frame(
      Parameter = p,
      Value     = draws[, p],
      Method    = m,
      stringsAsFactors = FALSE
    )
  }))
}))

ggplot(marg_long, aes(x = Value, color = Method)) +
  geom_density(linewidth = 0.7) +
  geom_line(data = analytic_df, aes(x = Value, y = Density),
            color = "black", linewidth = 1, linetype = "dashed",
            inherit.aes = FALSE) +
  facet_wrap(~ Parameter, scales = "free") +
  scale_color_manual(values = method_cols2) +
  labs(x = "Parameter value", y = "Density",
       title = "Marginal posteriors vs. analytical density (dashed black)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
```

#### Correlation recovery

The table below compares the estimated posterior correlation
\\\text{Cor}(\beta_1, \beta_2)\\ from each method against the analytical
value. Full-rank methods should recover the correlation closely;
mean-field methods by construction estimate zero (or near-zero)
correlation, regardless of the true posterior geometry.

``` r

analytic_cor <- post_cov2[2, 3] / sqrt(post_cov2[2, 2] * post_cov2[3, 3])

cor_df <- data.frame(
  Method = ok2,
  `Estimated correlation` = sapply(ok2, function(m) {
    fit <- results2[[m]]
    if (is.matrix(fit$Covar) && nrow(fit$Covar) >= 3) {
      v <- fit$Covar
      v[2, 3] / sqrt(v[2, 2] * v[3, 3])
    } else NA
  }),
  `Analytical correlation` = analytic_cor,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

knitr::kable(cor_df, digits = 3,
  caption = paste0("Posterior correlation between beta1 and beta2 ",
                   "(analytical: ", round(analytic_cor, 3), ")."))
```

#### Variance recovery

A related diagnostic compares each method’s estimated marginal standard
deviation against the analytical values \\\sqrt{\[\sigma^2 (X^\top
X)^{-1}\]\_{jj}}\\. This table makes the mean-field variance
underestimation concrete: because ADVI.mf and BBVI cannot represent the
correlation, they absorb less total variance into their marginals than
the true posterior contains.

``` r

analytic_sd <- sqrt(diag(post_cov2))

sd_df <- do.call(rbind, lapply(ok2, function(m) {
  fit <- results2[[m]]
  data.frame(
    Method = m,
    Parameter = parm.names2,
    `Estimated SD` = fit$Summary1[, "SD"],
    `Analytical SD` = analytic_sd,
    Ratio = fit$Summary1[, "SD"] / analytic_sd,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}))

knitr::kable(sd_df, digits = 3,
  caption = paste0("Posterior standard deviations: estimated vs. analytical. ",
                   "Ratio < 1 indicates underestimation."))
```

This comparison makes the structural limitation of mean-field methods
concrete. On problems where posterior correlations are strong, ADVI.mf
and BBVI produce point estimates that are correct (the mode is
unaffected by the variational family’s shape) but uncertainty
quantification that is systematically wrong. Pathfinder and the
full-rank methods avoid this failure mode because they maintain a full
covariance representation, either through the L-BFGS inverse Hessian
(Pathfinder), the Cholesky factor (ADVI.fr), or stochastic Hessian
averaging (Salimans2, NGD).

## Method selection guide

The choice of variational algorithm depends on the number of parameters,
the smoothness of the posterior, the importance of posterior
correlations, and computational budget. The table below summarizes the
key properties of each method.

| Method | Variational family | Variational parameters | Model evals per iter | Captures correlations | Multimodal | Recommended J |
|:---|:---|:---|:---|:--:|:--:|:--:|
| Salimans1 | Full Gaussian | J(J+1)/2 + J | O(J) | Yes | No | \< 100 |
| Salimans2 | Full Gaussian | J(J+1)/2 + J | O(J^2) | Yes | No | \< 20 |
| ADVI.mf | Diagonal Gaussian | 2J | O(J) | No | No | Any |
| ADVI.fr | Full Gaussian (Cholesky) | J + J(J+1)/2 | O(J) | Yes | No | \< 50 |
| BBVI | Diagonal Gaussian | 2J | O(S x J) | No | No | Any |
| SVGD | Nonparametric (particles) | K x J (particles) | O(K x J) | Yes (empirical) | Yes | \< 30 |
| NGD | Full Gaussian | J(J+1)/2 + J | O(J^2) | Yes | No | \< 20 |
| Pathfinder | Gaussian mixture (L-BFGS) | M x J (paths) | O(J) | Yes | No | Any |

Comparison of variational inference methods in lucifer. {.table}

The following guidelines summarize the practical recommendations.
**Pathfinder** is the recommended starting point for most problems: it
is fast, captures correlations via the L-BFGS inverse Hessian, scales
well with dimension, and provides a built-in Pareto \\\hat{k}\\
reliability diagnostic. For small problems with fewer than 20 parameters
and a smooth posterior, **Salimans2** or **NGD** should be tried first;
both capture correlations and converge rapidly. For moderate-dimensional
problems (20-100 parameters), **ADVI.mf** or **Salimans1** avoid the
quadratic cost of Hessian computation while still providing good mean
estimates, though ADVI.mf sacrifices correlation information. When
posterior correlations are important and \\J \leq 50\\, **ADVI.fr** is
the best reparameterization-based option. **BBVI** serves as a generic
fallback when other methods struggle, at the cost of higher gradient
variance and slower convergence. **SVGD** is uniquely suited for
multimodal posteriors but is computationally the most demanding and
should be reserved for problems where multimodality is genuinely
expected.

Regardless of the method chosen, the output of
[`VariationalBayes()`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
should be treated as an approximation. For final inference, consider
using the variational posterior as an informed initialization for MCMC
via
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
passing the estimated covariance matrix (scaled by \\2.38^2/J\\) to the
`Covar` argument. This two-stage strategy combines the speed of
variational methods with the asymptotic exactness of MCMC.

## References

**\[1\]** Salimans, T. and Knowles, D.A. (2013). *Fixed-form variational
posterior approximation through stochastic linear regression*. Bayesian
Analysis, 8(4), 837-882. [DOI:
10.1214/13-BA858](https://doi.org/10.1214/13-BA858)

**\[2\]** Khan, M.E. and Lin, W. (2017). *Conjugate-computation
variational inference: converting variational inference in non-conjugate
models to inferences in conjugate models*. Proceedings of the 20th
International Conference on Artificial Intelligence and Statistics
(AISTATS). [PMLR
54:878-887](http://proceedings.mlr.press/v54/khan17a.md)

**\[3\]** Kucukelbir, A., Tran, D., Ranganath, R., Gelman, A., and Blei,
D.M. (2017). *Automatic differentiation variational inference*. Journal
of Machine Learning Research, 18(14), 1-45. [URL:
https://jmlr.org/papers/v18/16-107.html](https://jmlr.org/papers/v18/16-107.html)

**\[4\]** Ranganath, R., Gerrish, S., and Blei, D.M. (2014). *Black box
variational inference*. Proceedings of the 17th International Conference
on Artificial Intelligence and Statistics (AISTATS). [PMLR
33:814-822](http://proceedings.mlr.press/v33/ranganath14.md)

**\[5\]** Liu, Q. and Wang, D. (2016). *Stein variational gradient
descent: a general purpose Bayesian inference algorithm*. Advances in
Neural Information Processing Systems 29 (NeurIPS). [arXiv:
1608.04471](https://arxiv.org/abs/1608.04471)

**\[6\]** Amari, S. (1998). *Natural gradient works efficiently in
learning*. Neural Computation, 10(2), 251-276. [DOI:
10.1162/089976698300017746](https://doi.org/10.1162/089976698300017746)

**\[7\]** Kingma, D.P. and Ba, J. (2015). *Adam: a method for stochastic
optimization*. Proceedings of the 3rd International Conference on
Learning Representations (ICLR). [arXiv:
1412.6980](https://arxiv.org/abs/1412.6980)

**\[8\]** Blei, D.M., Kucukelbir, A., and McAuliffe, J.D. (2017).
*Variational inference: a review for statisticians*. Journal of the
American Statistical Association, 112(518), 859-877. [DOI:
10.1080/01621459.2017.1285773](https://doi.org/10.1080/01621459.2017.1285773)

**\[9\]** Zhang, C., Butepage, J., Kjellstrom, H., and Mandt, S. (2019).
*Advances in variational inference*. IEEE Transactions on Pattern
Analysis and Machine Intelligence, 41(8), 2008-2026. [DOI:
10.1109/TPAMI.2018.2889774](https://doi.org/10.1109/TPAMI.2018.2889774)

**\[10\]** Zhang, L., Carpenter, B., Gelman, A., and Vehtari, A. (2022).
*Pathfinder: parallel quasi-Newton variational inference*. Journal of
Machine Learning Research, 23(306), 1-49. [URL:
https://jmlr.org/papers/v23/21-0889.html](https://jmlr.org/papers/v23/21-0889.html)

**\[11\]** Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J.
(2024). *Pareto smoothed importance sampling*. Journal of Machine
Learning Research, 25(72), 1-58. [URL:
https://jmlr.org/papers/v25/19-556.html](https://jmlr.org/papers/v25/19-556.html)

**\[12\]** Nocedal, J. and Wright, S.J. (2006). *Numerical
optimization*. 2nd ed. Springer. [DOI:
10.1007/978-0-387-40065-5](https://doi.org/10.1007/978-0-387-40065-5)
