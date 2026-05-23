# Frequentist information criteria from a lucifer fit

Extracts AIC, AICc, and BIC from any lucifer fit object by using the
deviance component (`Dev = -2 * log-lik`) and the parameter count.

## Usage

``` r
freq_ic(object, n = NULL)
```

## Arguments

- object:

  a fitted model object of class `laplace`, `iterquad`, `demonoid`, or
  `data_cloning`.

- n:

  sample size. If `NULL`, extracted from the object or its Data
  component.

## Value

A named numeric vector with elements `AIC`, `AICc`, and `BIC`.

## Details

The criteria are computed as: \$\$\mathrm{AIC} = -2\hat\ell + 2k\$\$
\$\$\mathrm{AICc} = \mathrm{AIC} + \frac{2k(k+1)}{n - k - 1}\$\$
\$\$\mathrm{BIC} = -2\hat\ell + k \log(n)\$\$ where \\\hat\ell\\ is the
maximized log-likelihood, \\k\\ is the number of parameters, and \\n\\
is the sample size.

## References

Akaike, H. (1974). A new look at the statistical model identification.
\*IEEE Transactions on Automatic Control\*, 19(6), 716-723. DOI:
10.1109/TAC.1974.1100705.

Schwarz, G. (1978). Estimating the dimension of a model. \*Annals of
Statistics\*, 6(2), 461-464. DOI: 10.1214/aos/1176344136.

## See also

[`lr_test`](https://robustecologies.github.io/lucifer/reference/lr_test.md),
[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- LaplaceApproximation(Model, Initial.Values, Data)
freq_ic(fit, n = 200)
} # }
```
