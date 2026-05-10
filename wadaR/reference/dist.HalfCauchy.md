# Half-Cauchy distribution

Density, distribution function, quantile function, and random generation
for the half-Cauchy distribution.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dhalfcauchy(x, scale = 25, log = FALSE)

phalfcauchy(q, scale = 25)

qhalfcauchy(p, scale = 25)

rhalfcauchy(n, scale = 25)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter \\\alpha\\, which must be positive.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- p:

  vector of probabilities.

- n:

  number of observations, which must be a positive integer that has
  length 1.

## Value

`dhalfcauchy` gives the density, `phalfcauchy` gives the distribution
function, `qhalfcauchy` gives the quantile function, and `rhalfcauchy`
generates random deviates.

See Details.

See Details.

See Details.

See Details.

## Details

The half-Cauchy distribution with scale \\\alpha=25\\ is a recommended,
default, weakly informative prior distribution for a scale parameter.

Implementation of `dhalfcauchy`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `phalfcauchy`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `qhalfcauchy`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rhalfcauchy`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dcauchy`](https://rdrr.io/r/stats/Cauchy.html)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- dhalfcauchy(1, 25)
x <- phalfcauchy(1, 25)
x <- qhalfcauchy(0.5, 25)
x <- rhalfcauchy(1, 25)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhalfcauchy
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving phalfcauchy
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qhalfcauchy
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhalfcauchy
} # }
```
