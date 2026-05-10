# Posterior predictive checks for `bayesquad`

Generates posterior predictive samples by evaluating the user model at
each posterior sample from a converged Bayesian quadrature fit.

## Usage

``` r
# S3 method for class 'bayesquad'
predict(object, Model, Data, CPUs = 1, Type = "PSOCK", ...)
```

## Arguments

- object:

  An object of class `bayesquad` with `Converged = TRUE`.

- Model:

  The model function used in `BayesianQuadrature`.

- Data:

  The data list used in `BayesianQuadrature`.

- CPUs:

  Number of CPUs for parallel evaluation (default 1).

- Type:

  Cluster type: `"PSOCK"` or `"MPI"`.

- ...:

  Additional arguments are unused.

## Value

An object of class `bayesquad.ppc` containing:

- y:

  Observed response vector.

- yhat:

  Matrix of posterior predictive samples (observations by posterior
  samples).

- Deviance:

  Vector of deviance values from posterior samples.

- monitor:

  Matrix of monitored values.

## Details

Posterior predictive sampling from a Bayesian quadrature fit produced by
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md).
Summary of the content is given below. For each posterior draw, the user
Model is re-evaluated and `yhat` is collected. The resulting posterior
predictive matrix enables posterior predictive checks via
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)
and diagnostic plotting via the `ppc` family of plot methods.

## References

O'Hagan, A. (1991). Bayes-Hermite quadrature. *Journal of Statistical
Planning and Inference*, 29(3), 245-260.
[doi:10.1016/0378-3758(91)90002-V](https://doi.org/10.1016/0378-3758%2891%2990002-V)

## See also

[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md),
[`summary.iterquad.ppc`](https://robustecologies.github.io/lucifer/reference/summary.iterquad.ppc.md),
[`plot.iterquad.ppc`](https://robustecologies.github.io/lucifer/reference/plot.iterquad.ppc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the BayesianQuadrature function for an example.
} # }
```
