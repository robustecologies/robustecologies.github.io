# Posterior predictive check summary for Laplace approximation

This may be used to summarize either new, unobserved instances of
\\\textbf{y}\\ (called \\\textbf{y}^{new}\\) or replicates of
\\\textbf{y}\\ (called \\\textbf{y}^{rep}\\). Either
\\\textbf{y}^{new}\\ or \\\textbf{y}^{rep}\\ is summarized, depending on
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md).

## Usage

``` r
# S3 method for class 'laplace.ppc'
summary(
  object = NULL,
  Categorical = FALSE,
  Rows = NULL,
  Discrep = NULL,
  d = 0,
  Quiet = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `laplace.ppc` is required.

- Categorical:

  Logical. If `TRUE`, then `y` and `yhat` are considered to be
  categorical (such as y=0 or y=1), rather than continuous.

- Rows:

  An optional vector of row numbers, for example `c(1:10)`. All rows
  will be estimated, but only these rows will appear in the summary.

- Discrep:

  A character string indicating a discrepancy test. `Discrep` defaults
  to `NULL`.

- d:

  This is an optional integer to be used with the `Discrep` argument
  above, and it defaults to `d=0`.

- Quiet:

  This logical argument defaults to `FALSE` and will print results to
  the console. When `TRUE`, results are not printed.

- ...:

  Additional arguments are unused.

## Value

This function returns a list with the following components:

- BPIC:

  The Bayesian Predictive Information Criterion (BPIC) was introduced by
  Ando (2007). \\BPIC = Dbar + 2pD\\. The goal is to minimize BPIC.

- Concordance:

  The percentage of records of y within the 95% quantile-based
  probability interval (see
  [`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md))
  of \\\textbf{y}^{rep}\\. Occurs only when \\\textbf{y}\\ is
  continuous.

- Mean.Lift:

  The mean of the record-level lifts, reported only when \\\textbf{y}\\
  is specified as categorical.

- Discrepancy.Statistic:

  Reported only if the `Discrep` argument receives a valid discrepancy
  measure.

- L.criterion:

  The L-criterion (Laud and Ibrahim, 1995) for model and variable
  selection.

- Monitor:

  A \\N \times 5\\ matrix with columns Mean, SD, LB, Median, and UB for
  monitored variables.

- Summary:

  When \\\textbf{y}\\ is continuous, a \\N \times 8\\ matrix with
  columns y, Mean, SD, LB, Median, UB, PQ, and Discrep.

## Details

This function summarizes an object of class `laplace.ppc`, which
consists of posterior predictive checks on either \\\textbf{y}^{new}\\
or \\\textbf{y}^{rep}\\, depending respectively on whether unobserved
instances of \\\textbf{y}\\ or the model sample of \\\textbf{y}\\ was
used in the
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md)
function. The deviance and monitored variables are also summarized.

The purpose of a posterior predictive check is to assess how well (or
poorly) the model fits the data, or to assess discrepancies between the
model and the data. For more information on posterior predictive checks,
see
<https://web.archive.org/web/20150215050702/http://www.bayesian-inference.com/posteriorpredictivechecks>.

When \\\textbf{y}\\ is continuous and known, this function estimates the
predictive concordance between \\\textbf{y}\\ and \\\textbf{y}^{rep}\\
as per Gelfand (1996), and the predictive quantile (PQ), which is for
record-level outlier detection used to calculate Gelfand's predictive
concordance.

When \\\textbf{y}\\ is categorical and known, this function estimates
the record-level lift, which is `p(yhat[i,] = y[i]) / [p(y = j) / n]`,
or the number of correctly predicted samples over the rate of that
category of \\\textbf{y}\\ in vector \\\textbf{y}\\.

A discrepancy measure is an approach to studying discrepancies between
the model and data (Gelman et al., 1996). Supported discrepancy measures
for continuous `y` include: `"Chi-Square"`, `"Chi-Square2"`,
`"Kurtosis"`, `"L-criterion"`, `"MASE"`, `"MSE"`, `"PPL"`,
`"Quadratic Loss"`, `"Quadratic Utility"`, `"RMSE"`, `"Skewness"`,
`"max(yhat[i,]) > max(y)"`, `"mean(yhat[i,]) > mean(y)"`,
`"mean(yhat[i,] > d)"`, `"mean(yhat[i,] > mean(y))"`,
`"min(yhat[i,]) < min(y)"`, `"round(yhat[i,]) = d"`, and
`"sd(yhat[i,]) > sd(y)"`. For categorical `y`: `"p(yhat[i,] != y[i])"`.

## References

Ando, T. (2007). "Bayesian Predictive Information Criterion for the
Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
*Biometrika*, 94(2), p. 443–458.
[doi:10.1093/biomet/asm017](https://doi.org/10.1093/biomet/asm017)

Gelfand, A. (1996). "Model Determination Using Sampling Based Methods".
In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in Markov
Chain Monte Carlo in Practice. Chapman and Hall: Boca Raton, FL. ISBN
0-412-07371-7.

Gelfand, A. and Ghosh, S. (1998). "Model Choice: A Minimum Posterior
Predictive Loss Approach". *Biometrika*, 85, p. 1–11.
[doi:10.1093/biomet/85.1.1](https://doi.org/10.1093/biomet/85.1.1)

Gelman, A., Meng, X.L., and Stern H. (1996). "Posterior Predictive
Assessment of Model Fitness via Realized Discrepancies". *Statistica
Sinica*, 6, p. 733–807.

Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model Selection".
*Journal of the Royal Statistical Society*, B 57, p. 247–262.
[doi:10.1111/j.2517-6161.1995.tb02025.x](https://doi.org/10.1111/j.2517-6161.1995.tb02025.x)

Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A.
(2002). "Bayesian Measures of Model Complexity and Fit (with
Discussion)". *Journal of the Royal Statistical Society*, B 64, p.
583–639.
[doi:10.1111/1467-9868.00353](https://doi.org/10.1111/1467-9868.00353)

## See also

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md),
and
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- LaplaceApproximation(Model, parm, Data, Method = "SPG")
pred <- predict(fit, Model, Data)
summary(pred, Discrep = "Chi-Square")
} # }
```
