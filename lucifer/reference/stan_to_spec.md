# Translate Stan model code to model_spec notation (experimental)

Parses the `model` block of a Stan program and attempts to translate
sampling statements (`y ~ distribution(args)`) into lucifer's
[`model_spec()`](https://robustecologies.github.io/lucifer/reference/model_spec.md)
notation. This translation is partial and covers only simple models with
standard distributions.

## Usage

``` r
stan_to_spec(stan_code)
```

## Arguments

- stan_code:

  Character string containing Stan model code.

## Value

A character string of model_spec notation. Warnings are emitted for
untranslatable constructs.

## Details

This function is experimental and has significant limitations. It can
translate simple models with standard distributions and scalar/vector
parameters. It cannot handle:

- Custom functions defined in the `functions` block

- Complex indexing, matrix operations, or ODE solvers

- Vectorized statements without explicit loops

- Truncated distributions or constrained parameters with complex bounds

- Models using `target +=` syntax instead of `y ~ dist()` syntax

The output should be treated as a starting point that requires manual
review and adjustment. Untranslatable constructs are included as
comments in the output.

## See also

[`jags_to_spec`](https://robustecologies.github.io/lucifer/reference/jags_to_spec.md),
[`model_spec`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
[`lucifer_stan`](https://robustecologies.github.io/lucifer/reference/lucifer_stan.md)

## Examples

``` r
if (FALSE) { # \dontrun{
stan_code <- "
data {
  int<lower=0> N;
  vector[N] y;
  matrix[N, 2] X;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  y ~ normal(X * beta, sigma);
}
"
spec_code <- stan_to_spec(stan_code)
cat(spec_code)
} # }
```
