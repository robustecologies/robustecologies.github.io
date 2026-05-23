# Hyperprior-g prior and Zellner's g-prior

Density and related functions for the hyperprior-g prior and zellner's
g-prior.

See Details.

See Details.

See Details.

## Usage

``` r
dhyperg(g, alpha = 3, log = FALSE)

dzellner(beta, g, sigma, X, log = FALSE)

rzellner(n, g, sigma, X)
```

## Arguments

- g:

  positive scalar hyperparameter.

- alpha:

  positive scale hyperhyperparameter, proper when \\\alpha \> 2\\.
  Default is 3.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- beta:

  regression effects \\\beta\\, a vector of length \\J\\.

- sigma:

  residual standard deviation \\\sigma\\, a positive scalar.

- X:

  full-rank \\N \times J\\ design matrix.

- n:

  number of random deviates to generate.

## Value

`dhyperg` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

## Details

Implementation of `dhyperg`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `dzellner`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rzellner`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhyperg
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dzellner
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rzellner
} # }
```
