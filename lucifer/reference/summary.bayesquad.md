# Summarize a `bayesquad` object

Provides an extended summary of a Bayesian quadrature fit, including GP
diagnostics, integral uncertainty quantification, and convergence
assessment.

## Usage

``` r
# S3 method for class 'bayesquad'
summary(object, ...)
```

## Arguments

- object:

  An object of class `bayesquad`.

- ...:

  Additional arguments are unused.

## Value

A list of class `summary.bayesquad` containing:

- Algorithm:

  The algorithm used.

- Converged:

  Logical convergence indicator.

- GP.diagnostics:

  GP fit diagnostics (number of evaluation points, kernel,
  hyperparameters, log marginal likelihood).

- Integral:

  Integral posterior summary (mean, variance, CI).

- Summary1:

  Point-estimate summary.

- Summary2:

  Posterior sample summary (if available).

## Details

Produces a tabular summary of a Bayesian quadrature fit produced by
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md).
Summary of the content is given below. The returned object has class
`summary.<class>` and carries marginal posterior quantiles, effective
sample sizes where applicable, and diagnostic flags. Printing the
summary object yields the human-readable table; subscripting it with
`` `$` `` exposes the underlying numeric matrix.

## References

O'Hagan, A. (1991). Bayes-Hermite quadrature. *Journal of Statistical
Planning and Inference*, 29(3), 245-260.
[doi:10.1016/0378-3758(91)90002-V](https://doi.org/10.1016/0378-3758%2891%2990002-V)

## See also

[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the BayesianQuadrature function for an example.
} # }
```
