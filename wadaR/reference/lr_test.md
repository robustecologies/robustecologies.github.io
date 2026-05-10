# Likelihood ratio test

Computes the likelihood ratio statistic \\\Lambda =
2(\hat\ell\_{\mathrm{full}} - \hat\ell\_{\mathrm{reduced}})\\ and tests
it against \\\chi^2(\mathrm{df})\\. The maximized log-likelihood is
extracted from the deviance component of any lucifer fit object
(`Dev = -2 * log-lik`).

## Usage

``` r
lr_test(object_full, object_reduced, df = NULL)
```

## Arguments

- object_full:

  fit object for the full (unrestricted) model.

- object_reduced:

  fit object for the reduced (restricted) model.

- df:

  degrees of freedom for the test. If `NULL`, inferred as the difference
  in parameter counts between the two models.

## Value

An object of class `lr_test` with elements `statistic`, `df`, `p.value`,
and log-likelihoods for both models.

## Details

For `data_cloning` objects the deviance at the final clone size is
divided by the number of clones to recover the un-cloned log-likelihood.
For all other objects the minimum deviance (at the mode) is used.

## References

Wilks, S. S. (1938). The large-sample distribution of the likelihood
ratio for testing composite hypotheses. \*Annals of Mathematical
Statistics\*, 9(1), 60-62. DOI: 10.1214/aoms/1177732360.

## See also

[`wald_test`](https://robustecologies.github.io/lucifer/reference/wald_test.md),
[`score_test`](https://robustecologies.github.io/lucifer/reference/score_test.md),
[`freq_ic`](https://robustecologies.github.io/lucifer/reference/freq_ic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Compare full and reduced models
fit_full    <- LaplaceApproximation(Model_full, IV_full, Data)
fit_reduced <- LaplaceApproximation(Model_reduced, IV_reduced, Data)
lr_test(fit_full, fit_reduced)
} # }
```
