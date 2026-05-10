# Simulation-based inference

Performs Bayesian inference for simulator-defined models using neural
density estimation. Unlike MCMC or variational methods, SBI does not
require an evaluable likelihood function. Instead, it trains a neural
network on simulated data to approximate the posterior, likelihood, or
likelihood-to-evidence ratio.

## Usage

``` r
SBI(
  simulator,
  prior,
  x_obs,
  method = "NPE",
  summary_fn = NULL,
  log_prior_fn = NULL,
  n_simulations = 10000L,
  network = NULL,
  backend = "auto",
  n_samples = 1000L,
  mcmc_algorithm = "NUTS",
  mcmc_iterations = 5000L,
  mcmc_specs = NULL,
  n_rounds = 1L,
  n_sims_per_round = NULL,
  hidden_layers = c(128L, 128L),
  n_components = 10L,
  activation = "selu",
  learning_rate = 5e-04,
  batch_size = 256L,
  max_epochs = 200L,
  patience = 20L,
  validation_fraction = 0.1,
  max_grad_norm = 10,
  seed = NULL,
  n_cores = 1L,
  pretrained_network = NULL,
  verbose = TRUE
)
```

## Arguments

- simulator:

  A function that takes a parameter vector and returns a simulated data
  vector (or summary statistics). The signature should be
  `function(theta)` where `theta` is a numeric vector of length equal to
  the number of parameters.

- prior:

  A function that takes no arguments and returns a single draw from the
  prior distribution as a numeric vector.

- x_obs:

  Numeric vector of observed summary statistics (or data) that the
  simulator produces.

- method:

  Character: `"NPE"` (default), `"NLE"`, `"NRE"`, `"SNPE"`, `"SNLE"`, or
  `"SNRE"`. Methods prefixed with `"S"` use sequential (multi-round)
  refinement.

- summary_fn:

  Optional function to reduce raw simulator output to summary
  statistics. Applied to both simulated and observed data.

- log_prior_fn:

  Optional function `function(theta)` returning the scalar log-prior
  density at `theta`. Used by NLE/NRE to evaluate the prior in the MCMC
  model closure. When `NULL` (default), an approximate log-prior is
  constructed from prior draws using a multivariate Gaussian with
  empirical support bounds.

- n_simulations:

  Integer: total number of simulations to generate. Default 10000.

- network:

  An
  [`sbi_network`](https://robustecologies.github.io/lucifer/reference/sbi_network.md)
  specification, or `NULL` for automatic architecture selection.

- backend:

  Character: `"auto"` (default), `"mdn"`, or `"torch"`.

- n_samples:

  Integer: number of posterior samples to draw (NPE only). Default 1000.

- mcmc_algorithm:

  Character: MCMC algorithm for NLE/NRE posterior sampling. Any
  algorithm supported by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  may be used. Default `"NUTS"`.

- mcmc_iterations:

  Integer: MCMC iterations for NLE/NRE. Default 5000.

- mcmc_specs:

  Specs list for
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  If `NULL`, defaults are generated automatically.

- n_rounds:

  Integer: number of sequential rounds. Values \> 1 activate sequential
  mode (SNPE/SNLE/SNRE). Default 1.

- n_sims_per_round:

  Integer or `NULL`: simulations per round. If `NULL`, defaults to
  `n_simulations / n_rounds`.

- hidden_layers:

  Integer vector: hidden layer sizes. Default `c(128, 128)`.

- n_components:

  Integer: MDN mixture components. Default 10.

- activation:

  Character: activation function. Default `"selu"`.

- learning_rate:

  Numeric: Adam learning rate. Default 5e-4.

- batch_size:

  Integer: mini-batch size. Default 256.

- max_epochs:

  Integer: maximum training epochs. Default 200.

- patience:

  Integer: early stopping patience. Default 20.

- validation_fraction:

  Numeric: fraction of data held out for validation. Default 0.1.

- max_grad_norm:

  Numeric: maximum gradient norm for gradient clipping during training.
  Default 10.0.

- seed:

  Integer or `NULL`: random seed for reproducibility. Default `NULL` (no
  seed).

- n_cores:

  Integer: number of cores for parallel simulation. Default 1
  (sequential).

- pretrained_network:

  A previously trained SBI network object (from a prior `SBI()` call's
  `$network` element) to use as initialization for training. Default
  `NULL`.

- verbose:

  Logical: print progress. Default `TRUE`.

## Value

An object of class `sbi` containing:

- Posterior:

  Matrix of posterior samples (n_samples x d).

- method:

  The SBI method used.

- network:

  The trained neural network object.

- training:

  List with loss curves and training diagnostics.

- simulations:

  List with theta and x matrices from training.

- mcmc_fit:

  For NLE/NRE: the `demonoid` object from
  [`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  `NULL` for NPE.

- standardization:

  Standardization parameters (means, sds).

- call:

  The matched call.

- Minutes:

  Wall-clock time in minutes.

## Details

Three method families are available, each training a different quantity
from simulated pairs \\(\theta_i, x_i)\\ where \\\theta_i \sim
p(\theta)\\ and \\x_i \sim p(x \mid \theta_i)\\.

**NPE** (neural posterior estimation) directly approximates
\\q\_\phi(\theta \mid x) \approx p(\theta \mid x)\\ using a conditional
density network. After training, posterior samples for any observation
\\x_o\\ are obtained by a single forward pass (amortized inference).
This is the recommended default.

**NLE** (neural likelihood estimation) approximates \\q\_\phi(x \mid
\theta) \approx p(x \mid \theta)\\, then uses MCMC to sample from
\\p(\theta \mid x_o) \propto q\_\phi(x_o \mid \theta) p(\theta)\\. The
trained neural likelihood is wrapped as a `Model(parm, Data)` closure
compatible with
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
giving access to all 62 MCMC algorithms.

**NRE** (neural ratio estimation) trains a classifier to distinguish
joint pairs \\(\theta, x) \sim p(\theta, x)\\ from marginal pairs
\\(\theta, x) \sim p(\theta)p(x)\\. The classifier logit approximates
the log likelihood-to-evidence ratio, which is used as the
log-likelihood in MCMC.

Sequential variants (SNPE, SNLE, SNRE) refine the proposal distribution
over multiple rounds, concentrating simulations near the observed data
for improved sample efficiency.

The default backend uses a C++ Mixture Density Network (MDN) with SELU
activations and full Cholesky-parametrized covariances. For problems
with more than ~20 parameters, the optional torch backend provides
normalizing flows. Both architectures are trained with the Adam
optimizer, gradient clipping, and early stopping.

## References

Deistler, M., Gloeckler, M., Boelts, J., et al. (2025). Simulation-based
inference: a practical guide. *arXiv:2411.17337*.

Papamakarios, G. & Murray, I. (2016). Fast epsilon-free inference of
simulation models with Bayesian conditional density estimation.
*NeurIPS* 29.

Hermans, J., Begy, V. & Louppe, G. (2020). Likelihood-free MCMC with
amortized approximate ratio estimators. *ICML*.

Cranmer, K., Brehmer, J. & Louppe, G. (2020). The frontier of
simulation-based inference. *PNAS* 117(48), 30055–30062.

Bishop, C.M. (1994). Mixture density networks. Aston University
Technical Report NCRG/94/004.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Normal-normal conjugate model (verifiable)
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)
x_obs <- c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0, 2.4, 1.7, 2.2, 2.6)

# NPE (amortized posterior)
fit <- SBI(simulator, prior, x_obs, method = "NPE",
           n_simulations = 5000, n_samples = 2000)
print(fit)
summary(fit)
plot(fit)

# NLE with NUTS (bridges into lucifer's MCMC)
fit_nle <- SBI(simulator, prior, x_obs, method = "NLE",
               mcmc_algorithm = "NUTS", mcmc_iterations = 3000)
plot(fit_nle, type = "training")
} # }
```
