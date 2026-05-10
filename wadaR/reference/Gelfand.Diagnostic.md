# Gelfand's convergence diagnostic

Gelfand et al. (1990) proposed a convergence diagnostic for Markov
chains. The `Gelfand.Diagnostic` function is an interpretation of
Gelfand's "thick felt-tip pen" MCMC convergence diagnostic. This
diagnostic plots a series of kernel density plots at \\k\\ intervals of
cumulative samples. Given a vector of \\S\\ samples from a marginal
posterior distribution, \\\theta\\, multiple kernel density lines are
plotted together, where each includes samples from a different interval.
It is assumed that
[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
iterations have been discarded.

Gelfand et al. (1990) assert that convergence is violated when the
plotted lines are farther apart than the width of a thick, felt-tip pen.
This depends on the size of the plot, and, of course, the pen. The
estimated width of a "thick felt-tip pen" is included as a black,
vertical line. The pen in `Gelfand.Diagnostic` is included for
historical reasons. This diagnostic requires numerous samples.

## Usage

``` r
Gelfand.Diagnostic(x, k = 3, pen = FALSE)
```

## Arguments

- x:

  a vector of marginal posterior samples, such as selected from the
  output of
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- k:

  the number \\k\\ of kernel density plots given cumulative intervals of
  samples. Defaults to \\k=3\\.

- pen:

  logical. When `pen=TRUE`, the thick felt-tip pen is included as a
  black, vertical line. Defaults to `FALSE`.

## Value

A plot is returned. The input `x` is returned invisibly.

## Details

Implementation of `Gelfand.Diagnostic`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## References

Gelfand, A.E., Hills, S., Racine-Poon, A., and Smith, A.F.M. (1990).
"Illustration of Bayesian Inference in Normal Data Models Using Gibbs
Sampling". *Journal of the American Statistical Association*, 85, p.
972–985.

## See also

[`burnin`](https://robustecologies.github.io/lucifer/reference/burnin.md)
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(1000)
Gelfand.Diagnostic(x)
} # }
```
