# Raftery and Lewis's diagnostic

Raftery and Lewis (1992) introduced an MCMC diagnostic that estimates
the number of iterations needed for a given level of precision in
posterior samples, as well as estimating burn-in, when quantiles are the
posterior summaries of interest.

## Usage

``` r
Raftery.Diagnostic(x, q = 0.025, r = 0.005, s = 0.95, eps = 0.001)
```

## Arguments

- x:

  an object of class `demonoid`. It attempts to use `Posterior2`, but
  when this is missing it uses `Posterior1`.

- q:

  the quantile to be estimated.

- r:

  the desired margin of error of the estimate, also called the accuracy.

- s:

  the probability of obtaining an estimate in the interval (q-r, q+r).

- eps:

  the precision required for the estimate of time to convergence.

## Value

An object of class `raftery` that is a list with the following
components:

- params:

  a vector containing the parameters `q`, `r`, and `s`.

- resmatrix:

  a 3-dimensional array containing the results: \\M\\ is the suggested
  burn-in, \\N\\ is the suggested number of iterations, \\Nmin\\ is the
  suggested number of iterations based on zero autocorrelation, and \\I
  = (M+N)/Nmin\\ is the "dependence factor". The dependence factor is
  interpreted as the proportional increase in the number of iterations
  attributable to autocorrelation. Highly autocorrelated chains (\> 5)
  are worrisome, and may be due to influential initial values, parameter
  correlations, or poor mixing.

## Details

In this MCMC diagnostic, a posterior quantile \\q\\ of interest is
specified. Next, an acceptable tolerance \\r\\ is specified for \\q\\,
which means that it is desired to measure \\q\\ with an accuracy of +/-
\\r\\. Finally, the user selects a probability \\s\\, which is the
probability of being within the interval \\(q-r, q+r)\\. The
`Raftery.Diagnostic` then estimates the number \\N\\ of iterations and
the number \\M\\ of burn-in iterations that are necessary to satisfy the
specified conditions regarding quantile \\q\\.

The diagnostic was designed to test a short, initial update, in which
the chains were called pilot chains, and the application was later
suggested for iterative use after any update as a general method for
pursuing convergence (Raftery and Lewis, 1996).

Results of the `Raftery.Diagnostic` differ depending on the chosen
quantile \\q\\. Estimates are conservative, so more iterations are
suggested than necessary.

## Note

The `Raftery.Diagnostic` function was adapted from the `raftery.diag`
function in the coda package, which was adapted from the FORTRAN program
'gibbsit', written by Steven Lewis.

## References

Raftery, A.E. and Lewis, S.M. (1992). "How Many Iterations in the Gibbs
Sampler?" In *Bayesian Statistics*, 4 (J.M. Bernardo, J.O. Berger, A.P.
Dawid and A.F.M. Smith, eds.). Oxford, U.K.: Oxford University Press, p.
763–773.

Raftery, A.E. and Lewis, S.M. (1992). "One Long Run with Diagnostics:
Implementation Strategies for Markov chain Monte Carlo". *Statistical
Science*, 7, p. 493–497.

Raftery, A.E. and Lewis, S.M. (1996). "Implementing MCMC". *In*
Practical Markov Chain Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and
S. Richardson, eds.). Chapman and Hall: Baton Rouge, FL.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`print.raftery`](https://robustecologies.github.io/lucifer/reference/print.raftery.md),
and
[`Thin`](https://robustecologies.github.io/lucifer/reference/Thin.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## After updating with lucifer, do:
# rd <- Raftery.Diagnostic(Fit)
# print(rd)
} # }
```
