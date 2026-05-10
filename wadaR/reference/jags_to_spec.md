# Translate JAGS model code to model_spec notation (experimental)

Parses a JAGS model block and attempts to translate stochastic and
deterministic nodes into lucifer's
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
notation. This translation is partial and covers only simple models.

## Usage

``` r
jags_to_spec(jags_code)
```

## Arguments

- jags_code:

  Character string containing JAGS model code.

## Value

A character string of model_spec notation. Warnings are emitted for
untranslatable constructs.

## Details

This function is experimental with significant limitations. It can
translate simple models with standard distributions. It cannot handle:

- Complex indexing with multiple subscripts

- The zeros/ones trick for custom likelihoods

- Mixture models with discrete latent allocations

- Nested for-loops with dependent indices

- Custom link functions or censoring syntax

JAGS uses the precision parameterization for the normal distribution.
The translator converts `dnorm(mu, tau)` to `Normal(mu, 1/sqrt(tau))`,
which may not always produce valid model_spec code (e.g., when tau is a
complex expression).

## See also

[`stan_to_spec`](https://robustecologies.github.io/lucifer/reference/stan_to_spec.md),
[`model_spec`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
[`lucifer_jags`](https://robustecologies.github.io/lucifer/reference/lucifer_jags.md)

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
spec_code <- jags_to_spec(jags_code)
cat(spec_code)
} # }
```
