# Frequentist summary table from a lucifer fit

Extracts maximum likelihood estimates, standard errors, test statistics,
p-values, and confidence intervals from any lucifer fit object,
producing a table that mirrors the coefficient output of
[`summary.lm`](https://rdrr.io/r/stats/summary.lm.html). Under flat
priors the posterior mode coincides with the MLE and the posterior
covariance approximates the inverse observed Fisher information, so a
Bayesian fit can be read as a frequentist one.

## Usage

``` r
freq_summary(object, conf_level = 0.95, df.residual = Inf)
```

## Arguments

- object:

  a fitted model object of class `laplace`, `iterquad`, `demonoid`, or
  `data_cloning`.

- conf_level:

  confidence level for the interval, default 0.95.

- df.residual:

  residual degrees of freedom. When finite, a \\t\\-distribution is used
  instead of the normal. Default `Inf` (z-test).

## Value

An object of class `freq_summary` containing a matrix with columns
`Estimate`, `Std.Error`, `z.value` (or `t.value`), `Pr(>|z|)` (or
`Pr(>|t|)`), `CI.lower`, and `CI.upper`.

## Details

For `laplace` and `iterquad` objects the mode and covariance matrix are
extracted directly. For `demonoid` objects (flat-prior MCMC) the
posterior mean of the stationary samples serves as the MLE and the
posterior covariance as the asymptotic variance estimate. For
`data_cloning` objects the MLE and asymptotic covariance are stored by
the data cloning algorithm itself.

When `df.residual` is finite the test statistic follows a \\t\\
distribution with `df.residual` degrees of freedom; otherwise a standard
normal (\\z\\) distribution is used.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
Rubin, D. B. (2013). \*Bayesian Data Analysis\*, 3rd ed. Chapman &
Hall/CRC. ISBN 978-1439840955.

van der Vaart, A. W. (1998). \*Asymptotic Statistics\*. Cambridge
University Press. ISBN 978-0521784504.

## See also

[`wald_test`](https://robustecologies.github.io/lucifer/reference/wald_test.md),
[`lr_test`](https://robustecologies.github.io/lucifer/reference/lr_test.md),
[`confint_compare`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Fit a linear regression with flat priors via Laplace approximation
fit <- LaplaceApproximation(Model, Initial.Values, Data)

## Extract frequentist summary
fs <- freq_summary(fit)
print(fs)
plot(fs)
} # }
```
