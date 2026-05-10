# Posterior predictive checks

This may be used to predict either new, unobserved instances of
\\\mathbf{y}\\ (called \\\mathbf{y}^{new}\\) or replicates of
\\\mathbf{y}\\ (called \\\mathbf{y}^{rep}\\), and then perform posterior
predictive checks. Either \\\mathbf{y}^{new}\\ or \\\mathbf{y}^{rep}\\
is predicted given an object of class `demonoid`, the model
specification, and data.

## Usage

``` r
# S3 method for class 'demonoid'
predict(object, Model, Data, CPUs = 1, Type = "PSOCK", ...)
```

## Arguments

- object:

  An object of class `demonoid` is required.

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

This function returns an object of class `demonoid.ppc` (where ppc
stands for posterior predictive checks). The returned object is a list
with the following components:

- y:

  This stores the vectorized form of \\\mathbf{y}\\, the dependent
  variable.

- yhat:

  This is a \\N \times S\\ matrix, where \\N\\ is the number of records
  of \\\mathbf{y}\\ and \\S\\ is the number of posterior samples.

- Deviance:

  This is a vector of predictive deviance.

## Details

This function passes each iteration of marginal posterior samples along
with data to `Model`, where the fourth component in the return list is
labeled `yhat`, and is a vector of expectations of \\\mathbf{y}\\, given
the samples, model specification, and data. Stationary samples are used
if detected, otherwise non-stationary samples will be used. To predict
\\\mathbf{y}^{rep}\\, simply supply the data set used to estimate the
model. To predict \\\mathbf{y}^{new}\\, supply a new data set instead
(though for some model specifications, this cannot be done, and
\\\mathbf{y}^{new}\\ must be specified in the `Model` function). If the
new data set does not have \\\mathbf{y}\\, then create `y` in the list
and set it equal to something sensible, such as `mean(y)` from the
original data set.

The variable `y` must be a vector. If instead it is matrix `Y`, then it
will be converted to vector `y`. The vectorized length of `y` or `Y`
must be equal to the vectorized length of `yhat`, the fourth component
of the return list of the `Model` function.

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

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit model
fit <- lucifer(Model, Data, Iterations = 1000, Status = 200, Thinning = 1,
  Algorithm = "CHARM", Specs = NULL)

# Predict
pred <- predict(fit, Model, Data)
summary(pred)
} # }
```
