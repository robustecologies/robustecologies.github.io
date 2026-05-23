# Heidelberger and Welch's MCMC convergence diagnostic

Heidelberger and Welch (1981; 1983) proposed a two-part MCMC convergence
diagnostic that calculates a test statistic (based on the Cramer-von
Mises test statistic) to accept or reject the null hypothesis that the
Markov chain is from a stationary distribution.

## Usage

``` r
Heidelberger.Diagnostic(x, eps = 0.1, pvalue = 0.05)
```

## Arguments

- x:

  an object of class `demonoid`. It attempts to use `Posterior2`, but
  when this is missing it uses `Posterior1`.

- eps:

  the target value for the ratio of halfwidth to sample mean.

- pvalue:

  the level of statistical significance.

## Value

An object of class `heidelberger`. This object is a \\J \times 6\\
matrix, and it is intended to be summarized with the
[`print.heidelberger`](https://robustecologies.github.io/lucifer/reference/print.heidelberger.md)
function. The object has \\J\\ rows, each of which corresponds to a
Markov chain. The column names are `stest`, `start`, `pvalue`, `htest`,
`mean`, and `halfwidth`. The `stest` column indicates convergence with a
one, and non-convergence with a zero, regarding the stationarity test.
When non-convergence is indicated, the remaining columns have missing
values. The `start` column indicates the starting iteration, and the
`pvalue` column shows the p-value associated with the first test. The
`htest` column indicates convergence for the halfwidth test. The `mean`
and `halfwidth` columns report the mean and halfwidth.

## Details

The Heidelberg and Welch MCMC convergence diagnostic consists of two
parts. In the first part, a chain of \\N\\ iterations is generated and
an alpha level is defined. The test statistic is calculated on the whole
chain to accept or reject the null hypothesis that the chain is from a
stationary distribution. If the null hypothesis is rejected, the first
10% of the chain is discarded and the test statistic is recalculated.
This process repeats until the null hypothesis is accepted or 50% of the
chain is discarded. If the test still rejects the null hypothesis, then
the chain fails the test and needs to be run longer.

In the second part, if the chain passes the first part, then the portion
of the chain that was not discarded is used to test the halfwidth
criterion. The halfwidth test calculates half the width of the (1 -
alpha)% probability interval (credible interval) around the mean. If the
ratio of the halfwidth and the mean is lower than `eps`, then the chain
passes the halfwidth test. Otherwise, the chain must be updated for more
iterations until sufficient accuracy is obtained. In order to avoid
problems caused by sequential testing, the test should not be repeated
too frequently. Heidelberger and Welch (1981) suggest increasing the run
length by a factor I \> 1.5 each time, so that the estimate has the
same, reasonably large, proportion of new data.

The Heidelberger and Welch MCMC convergence diagnostic conducts multiple
hypothesis tests. The number of potentially wrong results increases with
the number of non-independent hypothesis tests conducted.

The `Heidelberger.Diagnostic` is a univariate diagnostic that is usually
applied to each marginal posterior distribution. A multivariate form is
not included. By chance alone due to multiple independent tests, 5% of
the marginal posterior distributions should appear non-stationary when
stationarity exists. Assessing multivariate convergence is difficult.

## Note

The `Heidelberger.Diagnostic` function was adapted from the
`heidel.diag` function in the coda package.

## References

Heidelberger, P. and Welch, P.D. (1981). "A Spectral Method for
Confidence Interval Generation and Run Length Control in Simulations".
*Comm. ACM.*, 24, p. 233–245.

Heidelberger, P. and Welch, P.D. (1983). "Simulation Run Length Control
in the Presence of an Initial Transient". *Opns Res.*, 31, p. 1109–1144.

Schruben, L.W. (1982). "Detecting Initialization Bias in Simulation
Experiments". *Opns. Res.*, 30, p. 569–590.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`is.stationary`](https://robustecologies.github.io/lucifer/reference/is.stationary.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and
[`print.heidelberger`](https://robustecologies.github.io/lucifer/reference/print.heidelberger.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## After updating with lucifer, do:
# hd <- Heidelberger.Diagnostic(Fit)
# print(hd)
} # }
```
