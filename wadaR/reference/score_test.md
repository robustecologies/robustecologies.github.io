# Score test (experimental)

Computes the score (Lagrange multiplier) test statistic \\S =
U(\theta_0)^T I(\theta_0)^{-1} U(\theta_0)\\, where \\U(\theta_0)\\ is
the score (gradient of the log-likelihood) evaluated at the restricted
estimate and \\I(\theta_0)\\ is the Fisher information. The statistic is
compared against \\\chi^2(\mathrm{df})\\.

## Usage

``` r
score_test(Model, Data, parm_restricted, df, vcov_method = c("hessian", "opg"))
```

## Arguments

- Model:

  the model function, as used by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  the data list, as used by
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- parm_restricted:

  parameter vector at the restricted MLE.

- df:

  degrees of freedom (number of restrictions). Required.

- vcov_method:

  method for estimating the information matrix: `"hessian"` (default) or
  `"opg"`.

## Value

An object of class `score_test` with elements `statistic`, `df`,
`p.value`.

## Details

This function is marked experimental because the gradient is computed
via finite differences
([`partial`](https://robustecologies.github.io/lucifer/reference/Math.md)),
which introduces numerical noise. For most applications,
[`wald_test`](https://robustecologies.github.io/lucifer/reference/wald_test.md)
or
[`lr_test`](https://robustecologies.github.io/lucifer/reference/lr_test.md)
will be more reliable.

The score test has the advantage that only the restricted model needs to
be fitted; the gradient is evaluated at the restricted MLE without
requiring the full model fit.

## References

Rao, C. R. (1948). Large sample tests of statistical hypotheses
concerning several parameters with applications to problems of
estimation. \*Mathematical Proceedings of the Cambridge Philosophical
Society\*, 44(1), 50-57. DOI: 10.1017/S0305004100023987.

## See also

[`wald_test`](https://robustecologies.github.io/lucifer/reference/wald_test.md),
[`lr_test`](https://robustecologies.github.io/lucifer/reference/lr_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Score test at the restricted MLE
score_test(Model, Data, parm_restricted = c(1.0, 0, 0), df = 2)
} # }
```
