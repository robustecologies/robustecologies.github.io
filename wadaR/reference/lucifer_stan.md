# Fit a Stan model and return a lucifer demonoid object

Compiles and samples a Stan model using rstan, then converts the result
to a `demonoid` object with full post-processing compatibility. The
compiled Stan model is stored internally for refitting support
([`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO()`](https://robustecologies.github.io/lucifer/reference/LFO.md)).

## Usage

``` r
lucifer_stan(
  stan_code = NULL,
  stan_file = NULL,
  data,
  pars = NULL,
  chains = 4L,
  iter = 2000L,
  warmup = floor(iter/2),
  thin = 1L,
  seed = NULL,
  log_lik = "log_lik",
  control = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- stan_code:

  Character string containing Stan model code. Exactly one of
  `stan_code` or `stan_file` must be provided.

- stan_file:

  Path to a `.stan` file.

- data:

  Named list of data for the Stan model.

- pars:

  Character vector of parameters to monitor. If `NULL`, all parameters
  are monitored.

- chains:

  Integer number of chains (default 4).

- iter:

  Integer total iterations per chain including warmup (default 2000).

- warmup:

  Integer warmup iterations per chain (default `floor(iter / 2)`).

- thin:

  Integer thinning interval (default 1).

- seed:

  Integer random seed.

- log_lik:

  Character name of the generated quantity containing pointwise
  log-likelihoods (default `"log_lik"`). Set to `NULL` to disable
  automatic extraction.

- control:

  Named list of control parameters for Stan's sampler (e.g.,
  `list(adapt_delta = 0.99, max_treedepth = 15)`).

- verbose:

  Logical; print progress messages (default `TRUE`).

- ...:

  Additional arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

An object of class `demonoid` with `$.bridge$refit_fn` for
cross-validation support.

## Details

The function requires the `rstan` package to be installed. Stan model
compilation can take 30-60 seconds on first use; the compiled model is
cached by rstan for subsequent calls.

To enable
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`WAIC()`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
include a `generated quantities` block in your Stan code that computes
pointwise log-likelihoods:


      generated quantities {
        vector[N] log_lik;
        for (n in 1:N)
          log_lik[n] = normal_lpdf(y[n] | mu, sigma);
      }
      

## See also

[`as.demonoid`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
[`lucifer_jags`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md),
[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik`](https://robustecologies.github.io/lucifer/reference/log_lik.md)

## Examples

``` r
if (FALSE) { # \dontrun{
stan_code <- "
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | mu, sigma);
}
"

data_list <- list(N = 100, y = rnorm(100, 3, 2))
fit <- lucifer_stan(stan_code = stan_code, data = data_list)

# Full lucifer post-processing
print(fit)
summary(fit)
plot(fit)

# Model comparison
LOO(log_lik(fit))
} # }
```
