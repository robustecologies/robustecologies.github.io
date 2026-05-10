# The logit and inverse-logit functions

The logit and inverse-logit (also called the logistic function) are
provided.

See Details.

## Usage

``` r
invlogit(x)

logit(p)
```

## Arguments

- x:

  An object containing real values that will be transformed to the
  interval \[0,1\].

- p:

  An object containing probabilities p in the interval \[0,1\] that will
  be transformed to the real line.

## Value

`invlogit` returns probability `p`, and `logit` returns `x`.

See Details.

## Details

The `logit` function is the inverse of the sigmoid or logistic function,
and transforms a continuous value (usually probability \\p\\) in the
interval \[0,1\] to the real line (where it is usually the logarithm of
the odds). The `logit` function is \\\log(p / (1-p))\\.

The `invlogit` function (called either the inverse logit or the logistic
function) transforms a real number (usually the logarithm of the odds)
to a value (usually probability \\p\\) in the interval \[0,1\]. The
`invlogit` function is \\\frac{1}{1 + \exp(-x)}\\.

If \\p\\ is a probability, then \\\frac{p}{1-p}\\ is the corresponding
odds, while the `logit` of \\p\\ is the logarithm of the odds. The
difference between the logits of two probabilities is the logarithm of
the odds ratio. The derivative of probability \\p\\ in a logistic
function (such as `invlogit`) is: \\\frac{d}{dx} = p(1-p)\\.

In the lucifer package, it is common to re-parameterize a model so that
a parameter that should be in an interval can be updated from the real
line by using the `logit` and `invlogit` functions, though the
[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md)
function provides an alternative. For example, consider a parameter
\\\theta\\ that must be in the interval \[0,1\]. The algorithms in
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
are unaware of the desired interval, and may attempt \\\theta\\ outside
of this interval. One solution is to have the algorithms update
`logit(theta)` rather than `theta`. After `logit(theta)` is manipulated
by the algorithm, it is transformed via `invlogit(theta)` in the model
specification function, where \\\theta \in \[0,1\]\\.

Implementation of `logit`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`interval`](https://robustecologies.github.io/lucifer/reference/interval.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`plogis`](https://rdrr.io/r/stats/Logistic.html),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`qlogis`](https://rdrr.io/r/stats/Logistic.html), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- -5:5
p <- invlogit(x)
x <- logit(p)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving logit
} # }
```
