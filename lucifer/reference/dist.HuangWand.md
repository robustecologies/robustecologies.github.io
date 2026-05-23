# Huang-Wand distribution

Density and random generation for the Huang-Wand prior distribution for
a covariance matrix.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dhuangwand(x, nu = 2, a, A, log = FALSE)

rhuangwand(nu = 2, a, A)

dhuangwandc(x, nu = 2, a, A, log = FALSE)

rhuangwandc(nu = 2, a, A)
```

## Arguments

- x:

  \\k \times k\\ positive-definite covariance matrix (or Cholesky factor
  for `dhuangwandc`).

- nu:

  scalar degrees of freedom parameter \\\nu\\, default is 2
  (uninformative).

- a:

  positive-only vector of scale parameters of length \\k\\.

- A:

  positive-only vector of scale hyperparameters of length \\k\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dhuangwand` gives the density or other results.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dhuangwand`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rhuangwand`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `dhuangwandc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rhuangwandc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`dhalft`](https://robustecologies.github.io/lucifer/reference/dist.Halft.md),
[`dinvwishart`](https://robustecologies.github.io/lucifer/reference/dist.Inverse.Wishart.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhuangwand
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhuangwand
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dhuangwandc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rhuangwandc
} # }
```
