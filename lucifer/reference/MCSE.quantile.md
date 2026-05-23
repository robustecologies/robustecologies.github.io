# Monte Carlo standard error for quantiles

Estimates the Monte Carlo standard error of empirical quantiles from
MCMC output using the batch-means method. While the standard
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
estimates uncertainty in the posterior mean, `MCSE.quantile` estimates
uncertainty in posterior quantiles (e.g., credible interval endpoints).

## Usage

``` r
MCSE.quantile(x, q = c(0.025, 0.25, 0.5, 0.75, 0.975), batch.size = "sqrt")
```

## Arguments

- x:

  a numeric vector of posterior samples.

- q:

  a numeric vector of probabilities in \\(0, 1)\\ for which quantile
  MCSEs are desired. Defaults to `c(0.025, 0.25, 0.5, 0.75, 0.975)`.

- batch.size:

  either the batch size as a positive integer, or one of `"sqrt"`
  (default) or `"cuberoot"`, which set the batch size to
  \\\lfloor\sqrt{N}\rfloor\\ or \\\lfloor N^{1/3}\rfloor\\ respectively.

## Value

A matrix with rows corresponding to each requested quantile and two
columns: `est` (the point estimate of the quantile) and `se` (the Monte
Carlo standard error).

## Details

The batch-means estimator for quantile MCSE proceeds by dividing the
\\N\\ posterior samples into \\a\\ non-overlapping batches of size
\\b\\, computing the \\q\\-th quantile within each batch, and estimating
the variance of the quantile estimator as

\$\$\hat{\sigma}^2_q = \frac{b}{a - 1} \sum\_{k=1}^{a}
\left(\hat{Q}\_k(q) - \hat{Q}(q)\right)^2\$\$

where \\\hat{Q}\_k(q)\\ is the \\q\\-th quantile of batch \\k\\ and
\\\hat{Q}(q)\\ is the overall quantile. The MCSE is then
\\\sqrt{\hat{\sigma}^2_q / N}\\.

This approach follows the general framework of Flegal and Jones (2010)
for functional estimation. The batch size should be large enough to
capture the autocorrelation structure; the default `"sqrt"` is a
reasonable choice for most applications.

## References

Flegal, J.M. and Jones, G.L. (2010). "Batch means and spectral variance
estimators in Markov chain Monte Carlo". *The Annals of Statistics*,
38(2), p. 1034–1070. doi:10.1214/09-AOS735.

## See also

[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(5000)
MCSE.quantile(x)
MCSE.quantile(x, q = c(0.05, 0.5, 0.95))
} # }
```
