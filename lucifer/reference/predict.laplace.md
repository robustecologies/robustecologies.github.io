# Posterior predictive checks for Laplace approximation

This may be used to predict either new, unobserved instances of
\\\textbf{y}\\ (called \\\textbf{y}^{new}\\) or replicates of
\\\textbf{y}\\ (called \\\textbf{y}^{rep}\\), and then perform posterior
predictive checks. Either \\\textbf{y}^{new}\\ or \\\textbf{y}^{rep}\\
is predicted given an object of class `laplace`, the model
specification, and data. This function requires that posterior samples
were produced with
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).

## Usage

``` r
# S3 method for class 'laplace'
predict(object, Model, Data, CPUs = 1, Type = "PSOCK", ...)
```

## Arguments

- object:

  An object of class `laplace` is required.

- Model:

  The model specification function is required.

- Data:

  A data set in a list is required. The dependent variable is required
  to be named either `y` or `Y`.

- CPUs:

  This argument accepts an integer that specifies the number of central
  processing units (CPUs) of the multicore computer or computer cluster.
  This argument defaults to `CPUs=1`, in which parallel processing does
  not occur.

- Type:

  This argument specifies the type of parallel processing to perform,
  accepting either `Type="PSOCK"` or `Type="MPI"`.

- ...:

  Additional arguments are unused.

## Value

This function returns an object of class `laplace.ppc` (where ppc stands
for posterior predictive checks). The returned object is a list with the
following components:

- y:

  This stores \\\textbf{y}\\, the dependent variable.

- yhat:

  This is a \\N \times S\\ matrix, where \\N\\ is the number of records
  of \\\textbf{y}\\ and \\S\\ is the number of posterior samples.

- Deviance:

  This is a vector of length \\S\\, where \\S\\ is the number of
  independent posterior samples. Samples are obtained with the sampling
  importance resampling algorithm,
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md).

- monitor:

  This is a \\N \times S\\ matrix, where \\N\\ is the number of
  monitored variables and \\S\\ is the number of independent posterior
  samples. Samples are obtained with the sampling importance resampling
  algorithm,
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md).

## Details

Since Laplace Approximation characterizes marginal posterior
distributions with modes and variances, and posterior predictive checks
involve samples, the `predict.laplace` function requires the use of
independent samples of the marginal posterior distributions, provided by
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
when `sir=TRUE`.

The samples of the marginal posterior distributions of the target
distributions (the parameters) are passed along with the data to the
`Model` specification and used to draw samples from the deviance and
monitored variables. At the same time, the fourth component in the
returned list, which is labeled `yhat`, is a vector of expectations of
\\\textbf{y}\\, given the samples, model specification, and data. To
predict \\\textbf{y}^{rep}\\, simply supply the data set used to
estimate the model. To predict \\\textbf{y}^{new}\\, supply a new data
set instead (though for some model specifications, this cannot be done,
and \\\textbf{y}\_{new}\\ must be specified in the `Model` function). If
the new data set does not have \\\textbf{y}\\, then create `y` in the
list and set it equal to something sensible, such as `mean(y)` from the
original data set.

The variable `y` must be a vector. If instead it is matrix `Y`, then it
will be converted to vector `y`. The vectorized length of `y` or `Y`
must be equal to the vectorized length of `yhat`, the fourth component
of the returned list of the `Model` function.

Parallel processing may be performed when the user specifies `CPUs` to
be greater than one, implying that the specified number of CPUs exists
and is available. Parallelization may be performed on a multicore
computer or a computer cluster. Either a Simple Network of Workstations
(SNOW) or Message Passing Interface is used (MPI). With small data sets
and few samples, parallel processing may be slower, due to computer
network communication. With larger data sets and more samples, the user
should experience a faster run-time.

For more information on posterior predictive checks, see
<https://web.archive.org/web/20150215050702/http://www.bayesian-inference.com/posteriorpredictivechecks>.

## See also

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
and [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit Laplace approximation
fit <- LaplaceApproximation(Model, parm, Data, Method = "SPG")

# Predict
pred <- predict(fit, Model, Data)
summary(pred)
} # }
```
