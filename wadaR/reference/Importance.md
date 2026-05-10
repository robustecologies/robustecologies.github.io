# Variable importance

The `Importance` function considers variable importance (or predictor
importance) to be the effect that each variable has on replicates
\\\textbf{y}^{rep}\\ (or \\\textbf{Y}^{rep}\\) when the variable is
removed from the model by setting it equal to zero. Variable importance
is assessed via posterior predictive checks, making this a form of
sensitivity analysis useful for model revision, variable selection, and
model interpretation. Currently this function only tests the variable
importance of design matrix \\\textbf{X}\\.

## Usage

``` r
Importance(
  object,
  Model,
  Data,
  Categorical = FALSE,
  Discrep,
  d = 0,
  CPUs = 1,
  Type = "PSOCK"
)
```

## Arguments

- object:

  an object of class `demonoid`, `iterquad`, `laplace`, `pmc`, or `vb`.

- Model:

  the model specification function.

- Data:

  a list of data. The dependent variable must be named either `y` or
  `Y`, and the design matrix `X` must be present.

- Categorical:

  logical; if `TRUE`, `y` and `yhat` are treated as categorical.
  Defaults to `FALSE`.

- Discrep:

  optional discrepancy statistic. See
  [`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md)
  for details.

- d:

  optional integer for use with `Discrep`, defaults to `0`.

- CPUs:

  integer; the number of CPUs for parallel processing. Defaults to `1`.

- Type:

  character; the type of parallel processing, either `"PSOCK"` or
  `"MPI"`.

## Value

An object of class `importance`: a matrix with rows equal to the number
of columns in \\\textbf{X}\\ + 1 (including the full model) and 4
columns (BPIC, Concordance or Mean.Lift, Discrep, and L-criterion). For
non-categorical dependent variables, an attribute `S.L` (the calibration
number of the L-criterion) is attached.

## Details

The full model is first predicted and summarized via the appropriate
`predict` and `summary` methods. Each successive row in the output
corresponds to the removal of a column in design matrix \\\textbf{X}\\
(set to zero), showing the impact on posterior predictive checks. BPIC
(Ando, 2007) is the recommended default criterion for variable
importance, where larger values indicate a more important variable
because removing it worsened the fit.

The
[`plot.importance`](https://robustecologies.github.io/lucifer/reference/plot.importance.md)
function is available to plot the output according to BPIC, predictive
concordance, the selected discrepancy statistic, or the L-criterion.

## References

Ando, T. (2007). "Bayesian Predictive Information Criterion for the
Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
*Biometrika*, 94(2), p. 443–458.

Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model Selection".
*Journal of the Royal Statistical Society*, B 57, p. 247–262.

## See also

[`is.importance`](https://robustecologies.github.io/lucifer/reference/is.class.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`plot.importance`](https://robustecologies.github.io/lucifer/reference/plot.importance.md),
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md),
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`predict.iterquad`](https://robustecologies.github.io/lucifer/reference/predict.iterquad.md),
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# First, update the model with lucifer, creating object Fit.
Importance(Fit, Model, MyData, Discrep = "Chi-Square", CPUs = 1)
} # }
```
