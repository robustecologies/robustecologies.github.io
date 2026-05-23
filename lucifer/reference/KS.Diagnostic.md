# Kolmogorov-Smirnov convergence diagnostic

The Kolmogorov-Smirnov test is a nonparametric test of stationarity that
has been applied as an MCMC diagnostic (Brooks et al, 2003), such as to
the posterior samples from the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function. The first and last halves of the chain are compared. This test
assumes IID, which is violated in the presence of autocorrelation.

The `KS.Diagnostic` is a univariate diagnostic that is usually applied
to each marginal posterior distribution. A multivariate form is not
included. By chance alone due to multiple independent tests, 5% of the
marginal posterior distributions should appear non-stationary when
stationarity exists. Assessing multivariate convergence is difficult.

## Usage

``` r
KS.Diagnostic(x)
```

## Arguments

- x:

  a vector of posterior samples for which a Kolmogorov-Smirnov test will
  be applied that compares the first and last halves for stationarity.

## Value

A frequentist p-value, where stationarity is indicated when p \> 0.05.

## Details

There are two main approaches to using the Kolmogorov-Smirnov test as an
MCMC diagnostic. There is a version of the test that has been adapted to
account for autocorrelation (and is not included here). Otherwise, the
chain is thinned enough that autocorrelation is not present or is
minimized, in which case the two-sample Kolmogorov-Smirnov test is
applied. The CDFs of both samples are compared. The `ks.test` function
in base R is used.

The advantage of the Kolmogorov-Smirnov test is that it is easier and
faster to calculate. The disadvantages are that autocorrelation biases
results, and the test is generally biased on the conservative side
(indicating stationarity when it should not).

## References

Brooks, S.P., Giudici, P., and Philippe, A. (2003). "Nonparametric
Convergence Assessment for MCMC Model Selection". *Journal of
Computational and Graphical Statistics*. 12(1), p. 1–22.

## See also

[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md),
[`ks.test`](https://rdrr.io/r/stats/ks.test.html), and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(1000)
KS.Diagnostic(x)
} # }
```
