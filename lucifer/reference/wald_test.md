# Wald test for parameter restrictions

Tests the linear hypothesis \\R\theta = r\\ using the Wald statistic \\W
= (R\hat\theta - r)^T \[R \hat{V} R^T\]^{-1} (R\hat\theta - r)\\, which
is asymptotically distributed as \\\chi^2(q)\\ where \\q =
\mathrm{nrow}(R)\\.

## Usage

``` r
wald_test(object, R = NULL, r = NULL, vcov = NULL)
```

## Arguments

- object:

  a fitted model object of class `laplace`, `iterquad`, `demonoid`, or
  `data_cloning`.

- R:

  contrast matrix with `q` rows and `k` columns, where `k` is the number
  of parameters. If `NULL`, tests each parameter individually against
  zero.

- r:

  right-hand side vector of length `q`. Default is a vector of zeros.

- vcov:

  optional covariance matrix override. If `NULL`, the covariance is
  extracted from `object`.

## Value

An object of class `wald_test` with elements `statistic`, `df`,
`p.value`, `hypothesis`, and a print method.

## Details

When `R` and `r` are both `NULL`, each parameter is tested individually
against zero (equivalent to the z-tests in
[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)).

When `R` is a contrast matrix, the function performs a joint test of the
\\q\\ restrictions \\R\theta = r\\. The test statistic is compared
against the chi-squared distribution with \\q\\ degrees of freedom.

## References

Wald, A. (1943). Tests of statistical hypotheses concerning several
parameters when the number of observations is large. \*Transactions of
the American Mathematical Society\*, 54(3), 426-482.

## See also

[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md),
[`lr_test`](https://robustecologies.github.io/lucifer/reference/lr_test.md),
[`score_test`](https://robustecologies.github.io/lucifer/reference/score_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- LaplaceApproximation(Model, Initial.Values, Data)

## Test all parameters individually
wald_test(fit)

## Joint test: beta2 = beta3 = 0
R <- matrix(0, nrow = 2, ncol = length(Initial.Values))
R[1, 2] <- 1
R[2, 3] <- 1
wald_test(fit, R = R)
} # }
```
