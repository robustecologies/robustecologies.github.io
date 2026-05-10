# Coupled MCMC for unbiased estimation

Runs two coupled Markov chains using reflection-maximal coupling of
Gaussian random-walk proposals. When the chains meet at some random time
\\\tau\\, they remain coupled forever, enabling unbiased posterior
expectations via the debiasing formula of Jacob, O'Leary and Atchade
(2020) and convergence diagnostics through the meeting-time
distribution.

## Usage

``` r
coupled_mcmc(
  Model,
  Data,
  Initial.Values = NULL,
  Covar = NULL,
  max_iterations = 10000,
  k = 0,
  m = 0,
  tol = 1e-08,
  verbose = FALSE
)
```

## Arguments

- Model:

  Model function with signature `function(parm, Data)` returning a list
  with at least `LP` (log posterior), `Dev` (deviance), `Monitor`,
  `yhat`, and `parm`.

- Data:

  A named list of data. Must contain `parm.names`.

- Initial.Values:

  Numeric vector of starting parameter values. When `NULL`, a zero
  vector of length `length(Data$parm.names)` is used.

- Covar:

  Proposal covariance matrix. When `NULL` (default), an identity matrix
  scaled by \\(2.38)^2 / d\\ is used, following the optimal Gaussian
  random-walk scaling of Roberts, Gelman and Gilks (1997).

- max_iterations:

  Maximum number of iterations before the algorithm gives up waiting for
  the chains to meet. Default 10000.

- k:

  Lag parameter for the debiasing estimator. The estimator uses samples
  from iteration `k` onward, so larger values discard more burn-in.
  Default 0.

- m:

  Number of additional iterations to run after the chains have met,
  providing extra samples for variance reduction. Must satisfy \\m \geq
  k\\. Default 0.

- tol:

  Tolerance for detecting that two chains have met, i.e., \\\\X_t -
  Y_t\\\_\infty \< \texttt{tol}\\. Default `1e-8`.

- verbose:

  Logical. If `TRUE`, progress messages are printed to the console via
  [`message()`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An object of class `"coupled_mcmc"` containing:

- meeting_time:

  Integer. The iteration at which the chains first met, or `NA` if they
  did not meet within `max_iterations`.

- met:

  Logical. Whether the chains met.

- samples_X:

  Numeric matrix of dimension `(T x d)` storing the full trajectory of
  chain \\X\\.

- samples_Y:

  Numeric matrix of dimension `(T-1 x d)` storing the trajectory of
  chain \\Y\\ (lagged by 1 relative to \\X\\).

- log_posterior_X:

  Numeric vector of log-posterior values for chain \\X\\.

- log_posterior_Y:

  Numeric vector of log-posterior values for chain \\Y\\.

- k:

  The lag parameter used.

- m:

  The post-meeting extension used.

- parm.names:

  Parameter names from `Data$parm.names`.

- call:

  The matched call.

- minutes:

  Wall-clock time in minutes.

## Details

The algorithm runs two Metropolis-Hastings chains \\(X_t)\\ and
\\(Y_t)\\ with a shared proposal covariance \\\Sigma\\. The chain \\Y\\
is initialized one step behind \\X\\ (a lag of 1). At each iteration, a
reflection-maximal coupling of the proposal kernels ensures that the
proposals for both chains agree with the highest possible probability.
Specifically, given current states \\x\\ and \\y\\, a proposal \\z \sim
N(x, \Sigma)\\ is drawn for chain \\X\\. With probability \$\$p =
\min\left(1, \frac{N(z; y, \Sigma)}{N(z; x, \Sigma)}\right)\$\$ the same
proposal \\z\\ is used for chain \\Y\\; otherwise a reflection of the
standardized innovation across the separating hyperplane between \\x\\
and \\y\\ generates the residual proposal for \\Y\\. Each chain
independently accepts or rejects via the standard Metropolis-Hastings
ratio.

When the chains meet at time \\\tau\\ (i.e., \\X\_\tau = Y\_{\tau-1}\\),
they are set equal and evolve identically thereafter. The meeting time
\\\tau\\ is a random variable whose distribution characterizes the
mixing rate of the underlying kernel.

The debiasing estimator for any function \\h\\ is \$\$H_k = h(X_k) +
\sum\_{t=k+1}^{\tau-1} \left\[ h(X_t) - h(Y\_{t-1}) \right\]\$\$ which
is an unbiased estimator of \\E\_\pi\[h(X)\]\\ regardless of the choice
of burn-in \\k\\, provided \\E\[\tau\] \< \infty\\. When \\m \> k\\, the
time-averaged estimator \$\$H\_{k:m} = \frac{1}{m - k + 1}
\sum\_{l=k}^{m} H_l\$\$ has lower variance.

## References

Jacob, P. E., O'Leary, J. and Atchade, Y. F. (2020). Unbiased Markov
chain Monte Carlo methods with couplings. *Journal of the Royal
Statistical Society, Series B*, **82**(3), 543–600.
[doi:10.1111/rssb.12336](https://doi.org/10.1111/rssb.12336)

Roberts, G. O., Gelman, A. and Gilks, W. R. (1997). Weak convergence and
optimal scaling of random walk Metropolis algorithms. *Annals of Applied
Probability*, **7**(1), 110–120.
[doi:10.1214/aoap/1034625254](https://doi.org/10.1214/aoap/1034625254)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
for standard MCMC sampling.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple normal model
N <- 50
y <- rnorm(N, 3, 1)
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
    dcauchy(sigma, 0, 25, log = TRUE)
  yhat <- rep(mu, Data$N)
  Monitor <- LP
  list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
    yhat = yhat, parm = parm)
}

# Single coupled run
fit <- coupled_mcmc(Model, Data, Initial.Values = c(0, 0),
  max_iterations = 5000, k = 100, verbose = TRUE)
print(fit)
summary(fit)
plot(fit)

# Multiple runs for meeting-time distribution
fits <- replicate(50, coupled_mcmc(Model, Data,
  Initial.Values = c(0, 0), max_iterations = 5000, k = 100),
  simplify = FALSE)
plot_meeting_times(fits)
} # }
```
