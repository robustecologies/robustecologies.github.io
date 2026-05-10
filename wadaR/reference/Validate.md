# Holdout validation

This function performs holdout validation on an object of class
`demonoid` or `pmc`, given both a modeled and validation data set.

## Usage

``` r
Validate(object, Model, Data, plot = FALSE, PDF = FALSE)
```

## Arguments

- object:

  An object of class `demonoid` or `pmc`.

- Model:

  A model specification function for
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  or
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- Data:

  A list that contains two lists of data, as specified for
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  The first component in the list is the list of modeled data, and the
  second component in the list is the list of validation data.

- plot:

  Logical. When `plot=TRUE`, two plots are displayed. The upper plot
  shows the density of the modeled deviance in black and the density of
  the validation deviance in red. The lower plot shows the density of
  the change in deviance in gray. The `plot` argument defaults to
  `FALSE`.

- PDF:

  Logical. When `PDF=TRUE` (and `plot=TRUE`), the plot is saved as a
  .pdf file. The `PDF` argument defaults to `FALSE`.

## Value

A list with three components. The first two components are also lists.
Each list consists of `y`, `yhat`, and `Deviance`. The third component
is a matrix that reports the expected deviance, pD, and BPIC. The object
is of class `demonoid.val` for
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
or `pmc.val` when associated with
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Details

There are numerous ways to validate a model. In this context, validation
means to assess the predictive performance of a model on out-of-sample
data. If reasonable, leave-one-out cross-validation (LOOCV) via the
conditional predictive ordinate (CPO) should be considered when using
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
or [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).
For more information on CPO, see the accompanying vignettes entitled
"Bayesian Inference" and "Examples". CPO is unavailable when using
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

For
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
it is recommended that the user perform holdout validation by comparing
posterior predictive checks, comparing the differences in the specified
discrepancy measure.

When LOOCV is unreasonable, popular alternatives include k-fold
cross-validation and holdout validation. Although k-fold
cross-validation is not performed explicitly here, the user may
accomplish it with some effort. Of these methods, holdout validation
includes the most bias, but is the most common in applied use, since
only one model is fitted, rather than \\k-1\\ models in k-fold
cross-validation. The `Validate` function performs holdout validation.

For holdout validation, the observed data is sampled randomly into two
data sets of approximately equal size, or three data sets that consists
of two data sets of approximately equal size and a remainder data set.
Of the two data sets approximately equal in size, one is called the
modeled (or training) data set, and the other is called the validation
(or test) data set. The modeled data set is used when updating the
model. After the model is updated, both data sets are predicted in the
`Validate` function, given the model. Predictive loss is estimated for
the validation data set, relative to the modeled data set.

Predictive loss is associated with overfitting, differences between the
model and validation data set, or model misspecification. Bayesian
inference is reputed to be much more robust to overfitting than
frequentist inference.

There are many ways to measure predictive loss, and within each
approach, there are usually numerous possible loss functions. The
log-likelihood of the model is a popular approximate utility function,
and consequently, the deviance of the model is a popular loss function.

A vector of model-level (rather than record-level) deviance samples is
returned with each object of class `demonoid` or `pmc`. The `Validate`
function obtains this vector for each data set, and then calculates the
Bayesian Predictive Information Criterion (BPIC), as per Ando (2007).
BPIC is a variation of the Deviance Information Criterion (DIC) that has
been modified for predictive distributions. For more information on DIC
(Spiegelhalter et al., 2002), see the accompanying vignette entitled
"Bayesian Inference". The goal is to minimize BPIC.

When DIC is applied after the model, such as with a predictive
distribution, it is positively biased, or too small. The bias is due to
the same data \\\textbf{y}\\ being used both to construct the posterior
distributions and to evaluate pD, the penalty term for model complexity.
For example, for validation data set \\\textbf{y}\_{new}\\, BPIC is:

\$\$BPIC = -2\mathrm{log}\[p(\textbf{y}\_{new}\|\textbf{y},\Theta)\] +
2pD\$\$

When `plot=TRUE`, the distributions of the modeled and validation
deviances are plotted above, and the lower plot is the modeled deviance
subtracted from the validation deviance. When positive, this
distribution of the change in deviance is the loss in predictive
deviance associated with moving from the modeled data set to the
validation data set.

After using the `Validate` function, the user is encouraged to perform
posterior predictive checks on each data set via the
[`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md)
or
[`summary.pmc.ppc`](https://robustecologies.github.io/lucifer/reference/summary.pmc.ppc.md)
function.

## References

Ando, T. (2007). "Bayesian Predictive Information Criterion for the
Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
*Biometrika*, 94(2), p. 443–458.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A.
(2002). "Bayesian Measures of Model Complexity and Fit (with
Discussion)". *Journal of the Royal Statistical Society*, B 64, p.
583–639.

## See also

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Given an object called Fit of class demonoid, a Model specification,
## and a modeled data set (MyData.M) and validation data set (MyData.V):
Validate(Fit, Model, Data=list(MyData.M=MyData.M, MyData.V=MyData.V))
} # }
```
