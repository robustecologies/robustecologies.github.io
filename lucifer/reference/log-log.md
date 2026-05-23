# The log-log and complementary log-log functions

The log-log and complementary log-log functions, as well as the inverse
functions, are provided.

See Details.

See Details.

See Details.

## Usage

``` r
loglog(p)

invloglog(x)

cloglog(p)

invcloglog(x)
```

## Arguments

- p:

  A vector of probabilities p in the interval \[0,1\] that will be
  transformed to the real line.

- x:

  A vector of real values that will be transformed to the interval
  \[0,1\].

## Value

`cloglog` returns `x`, `invcloglog` and `invloglog` return probability
`p`, and `loglog` returns `x`.

See Details.

See Details.

See Details.

## Details

The logit and probit links are symmetric, because the probabilities
approach zero or one at the same rate. The log-log and complementary
log-log links are asymmetric. Complementary log-log links approach zero
slowly and one quickly. Log-log links approach zero quickly and one
slowly. Either the log-log or complementary log-log link will tend to
fit better than logistic and probit, and are frequently used when the
probability of an event is small or large. A mixture of the two links,
the log-log and complementary log-log is often used, where each link is
weighted. The reason that logit is so prevalent is because logistic
parameters can be interpreted as odds ratios.

Implementation of `invloglog`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `cloglog`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `invcloglog`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- -5:5
p <- invloglog(x)
x <- loglog(p)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving invloglog
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving cloglog
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving invcloglog
} # }
```
