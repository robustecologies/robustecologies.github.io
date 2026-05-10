# Categorical distribution

Density, quantile function, and random generation for the categorical
distribution with probabilities parameter \\p\\.

See Details.

See Details.

See Details.

## Usage

``` r
dcat(x, p, log = FALSE)

qcat(pr, p, lower.tail = TRUE, log.pr = FALSE)

rcat(n, p)
```

## Arguments

- x:

  vector of discrete data with \\k\\ discrete categories, of length
  \\n\\. Also accepts an \\n \times k\\ indicator matrix.

- p:

  vector of length \\k\\ or \\n \times k\\ matrix of probabilities.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- pr:

  vector of probabilities, or log-probabilities.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\Pr\[X \le x\]\\.

- log.pr:

  logical; if `TRUE`, probabilities \\pr\\ are given as \\\log(pr)\\.

- n:

  number of observations, which must be a positive integer that has
  length 1.

## Value

`dcat` gives the density, `qcat` gives the quantile function, and `rcat`
generates random deviates.

See Details.

See Details.

See Details.

## Details

Also called the discrete distribution, the categorical distribution
describes the result of a random event that can take on one of \\k\\
possible outcomes, with the probability \\p\\ of each outcome separately
specified. The conjugate prior is the Dirichlet distribution.

Implementation of `dcat`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qcat`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rcat`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`as.indicator.matrix`](https://robustecologies.github.io/lucifer/reference/Matrices.md),
[`ddirichlet`](https://robustecologies.github.io/lucifer/reference/dist.Dirichlet.md),
[`dmultinom`](https://rdrr.io/r/stats/Multinom.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
dcat(x=1, p=c(0.3,0.3,0.4))
rcat(n=10, p=c(0.1,0.3,0.6))
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dcat
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qcat
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rcat
} # }
```
