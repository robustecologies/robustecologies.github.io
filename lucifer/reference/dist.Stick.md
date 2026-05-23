# Truncated stick-breaking prior distribution

Density and related functions for the truncated stick-breaking prior
distribution.

See Details.

See Details.

## Usage

``` r
dStick(theta, gamma, log = FALSE)

rStick(M, gamma)
```

## Arguments

- theta:

  vector of length \\M-1\\, where \\M\\ is the truncated number of
  possible mixture components.

- gamma:

  scalar, usually gamma-distributed.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- M:

  integer equal to one less than the truncated number of possible
  mixture components.

## Value

`dStick` gives the density, and other functions provide related
computations.

See Details.

See Details.

## Details

Implementation of `dStick`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rStick`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dStick
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rStick
} # }
```
