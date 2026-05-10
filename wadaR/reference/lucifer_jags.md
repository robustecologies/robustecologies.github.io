# Fit a JAGS model and return a lucifer demonoid object

Fits a JAGS model using rjags, then converts the result to a `demonoid`
object with full post-processing compatibility. The model code and data
are stored internally for refitting support
([`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)).

## Usage

``` r
lucifer_jags(
  jags_code = NULL,
  jags_file = NULL,
  data,
  inits = NULL,
  parameters.to.save,
  n.chains = 4L,
  n.iter = 10000L,
  n.burnin = floor(n.iter/2),
  n.thin = 1L,
  n.adapt = 1000L,
  log_lik = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- jags_code:

  Character string containing JAGS model code. Exactly one of
  `jags_code` or `jags_file` must be provided.

- jags_file:

  Path to a JAGS model file (`.bug` or `.jags`).

- data:

  Named list of data for the JAGS model.

- inits:

  Initial values. Either a named list, a function returning a list, or a
  list of such (one per chain). If `NULL`, JAGS generates initial values
  automatically.

- parameters.to.save:

  Character vector of parameter names to monitor.

- n.chains:

  Integer number of chains (default 4).

- n.iter:

  Integer total iterations per chain (default 10000).

- n.burnin:

  Integer burn-in iterations (default `floor(n.iter / 2)`).

- n.thin:

  Integer thinning interval (default 1).

- n.adapt:

  Integer number of adaptation iterations (default 1000).

- log_lik:

  An optional \\N \times S\\ matrix of pointwise log-likelihoods. JAGS
  does not compute these automatically; if needed, monitor the
  individual log-density contributions in the model code.

- verbose:

  Logical; print progress messages (default `TRUE`).

- ...:

  Additional arguments passed to
  [`rjags::jags.model()`](https://rdrr.io/pkg/rjags/man/jags.model.html).

## Value

An object of class `demonoid` with `$.bridge$refit_fn` for
cross-validation support.

## Details

The function requires the `rjags` package (and a working JAGS
installation) to be available.

Unlike Stan, JAGS does not automatically compute pointwise
log-likelihoods. To enable
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
include pointwise log-density nodes in the model and monitor them, then
pass the resulting matrix via `log_lik`. Alternatively, compute the
pointwise log-likelihoods from the posterior samples after fitting.

JAGS uses the precision parameterization for the normal distribution
(`dnorm(mu, tau)` where `tau = 1/sigma^2`). This is handled
transparently in the model code; posterior samples are returned in
whichever parameterization the model uses.

## See also

[`as.demonoid`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
[`lucifer_stan`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md),
[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik`](https://robustecologies.github.io/lucifer/reference/log_lik.md)

## Examples

``` r
if (FALSE) { # \dontrun{
jags_code <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu, tau)
  }
  mu ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
"

data_list <- list(N = 100, y = rnorm(100, 3, 2))
fit <- lucifer_jags(
    jags_code = jags_code,
    data = data_list,
    parameters.to.save = c("mu", "sigma"),
    n.iter = 10000
)

# Full lucifer post-processing
print(fit)
summary(fit)
plot(fit)
} # }
```
