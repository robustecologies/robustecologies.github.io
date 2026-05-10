# Bernoulli distribution

Density, distribution function, quantile function, and random generation
for the Bernoulli distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dbern(x, prob, log = FALSE)

pbern(q, prob, lower.tail = TRUE, log.p = FALSE)

qbern(p, prob, lower.tail = TRUE, log.p = FALSE)

rbern(n, prob)
```

## Arguments

- x, q:

  vector of quantiles.

- prob:

  probability of success on each trial.

- log, log.p:

  logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\,
  otherwise \\Pr\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dbern` gives the density, `pbern` gives the distribution function,
`qbern` gives the quantile function, and `rbern` generates random
deviates.

See Details.

See Details.

See Details.

See Details.

## Details

The Bernoulli distribution is a binomial distribution with \\n=1\\. The
categorical distribution is the generalization of the Bernoulli
distribution for variables with more than two discrete values. The beta
distribution is the conjugate prior distribution of the Bernoulli
distribution.

Implementation of `dbern`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `pbern`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qbern`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rbern`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`dbinom`](https://rdrr.io/r/stats/Binomial.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
dbern(1, 0.7)
rbern(10, 0.5)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dbern
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pbern
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qbern
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rbern
} # }
```
