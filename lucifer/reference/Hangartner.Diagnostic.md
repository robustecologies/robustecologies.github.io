# Hangartner's convergence diagnostic

Hangartner et al. (2011) proposed a convergence diagnostic for discrete
Markov chains. A simple Pearson's Chi-squared test for two or more
non-overlapping periods of a discrete Markov chain is a reliable
diagnostic of convergence. It does not rely upon the estimation of
spectral density, on suspect normality assumptions, or determining
overdispersion within a small number of outcomes, all of which can be
problematic with discrete measures. A discrete Markov chain is split
into two or more non-overlapping windows. Two windows are recommended,
and results may be sensitive to the number of selected windows, as well
as sample size. As such, a user may try several window configurations
before concluding there is no evidence of non-convergence.

As the number of discrete events in the sample space increases, this
diagnostic becomes less appropriate and standard diagnostics become more
appropriate.

## Usage

``` r
Hangartner.Diagnostic(x, J = 2)
```

## Arguments

- x:

  a vector of marginal posterior samples of a discrete Markov chain,
  such as selected from the output of
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- J:

  the number of windows to be used. Defaults to \\J=2\\.

## Value

An object of class `hangartner`, including the output from a Pearson's
Chi-squared test. A frequentist p-value less than or equal to 0.05 is
usually considered to be indicative of non-convergence.

## Details

Implementation of `Hangartner.Diagnostic`. Refer to the package
vignettes and the cited references for a complete algorithmic and
mathematical description.

## References

Hangartner, D., Gill, J., and Cranmer, S., (2011). "An MCMC Diagnostic
for Purely Discrete Parameters". Paper presented at the annual meeting
of the Southern Political Science Association, Hotel InterContinental,
New Orleans, Louisiana Online.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and
[`TransitionMatrix`](https://robustecologies.github.io/lucifer/reference/Matrices.md).

## Examples

``` r
if (FALSE) { # \dontrun{
N <- 1000
K <- 3
x <- rcat(N, rep(1/K, K))
hd <- Hangartner.Diagnostic(x, J=2)
hd
} # }
```
