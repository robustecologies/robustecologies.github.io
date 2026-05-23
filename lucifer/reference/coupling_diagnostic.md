# L-lag coupling diagnostics for MCMC convergence

Computes non-asymptotic upper bounds on the total variation distance
between the chain's marginal distribution and the target, using coupled
chains with time lag L. Based on Biswas, Jacob and Vanetti (2019).

## Usage

``` r
coupling_diagnostic(
  Model,
  Data,
  Initial.Values = NULL,
  Covar = NULL,
  K = 100,
  L_max = 500,
  max_iterations = 10000
)
```

## Arguments

- Model:

  Model function (see
  [`coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md)).

- Data:

  Data list (see
  [`coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/coupled_mcmc.md)).

- Initial.Values:

  Initial parameter values.

- Covar:

  Proposal covariance matrix.

- K:

  Number of independent coupled chain pairs to run.

- L_max:

  Maximum lag to evaluate. Default 500.

- max_iterations:

  Maximum iterations per coupled pair. Default 10000.

## Value

A list of class `coupling_diagnostic` with components:

- tv_upper_bound:

  Numeric vector of TV upper bounds at lags 1:L_max

- meeting_times:

  Integer vector of meeting times from K pairs

- L_max:

  Maximum lag evaluated

- K:

  Number of coupled pairs

## Details

The TV upper bound at lag L is estimated as: \$\$\bar{U}(L) =
\frac{1}{K} \sum\_{k=1}^{K} \max(0, \lceil (\tau_k - L) / L \rceil)\$\$
where \\\tau_k\\ is the meeting time of the k-th coupled pair. This
provides a computable, non-asymptotic certificate that the chain has
mixed, unlike heuristic diagnostics like R-hat.

## References

Biswas, N., Jacob, P.E. and Vanetti, P. (2019). "Estimating convergence
of Markov chains with L-lag couplings". *NeurIPS* 32.

## See also

[`plot.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/plot.coupled_mcmc.md),
[`plot.coupling_diagnostic`](https://robustecologies.github.io/lucifer/reference/plot.coupling_diagnostic.md),
[`print.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/print.coupled_mcmc.md),
[`print.coupling_diagnostic`](https://robustecologies.github.io/lucifer/reference/print.coupling_diagnostic.md),
[`summary.coupled_mcmc`](https://robustecologies.github.io/lucifer/reference/summary.coupled_mcmc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
Model <- function(parm, Data) {
    LP <- sum(dnorm(parm, 0, 1, log = TRUE))
    list(LP = LP, Dev = -2*LP, Monitor = LP,
         yhat = parm, parm = parm)
}
Data <- list(parm.names = c("x1", "x2"), mon.names = "LP")
diag <- coupling_diagnostic(Model, Data, c(0, 0), K = 50)
print(diag)
plot(diag)
} # }
```
