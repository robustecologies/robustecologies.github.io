# Truncated distributions

Density and related functions for the truncated distributions.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
dtrunc(x, spec, a = -Inf, b = Inf, log = FALSE, ...)

extrunc(spec, a = -Inf, b = Inf, ...)

ptrunc(x, spec, a = -Inf, b = Inf, ...)

qtrunc(p, spec, a = -Inf, b = Inf, ...)

rtrunc(n, spec, a = -Inf, b = Inf, ...)

vartrunc(spec, a = -Inf, b = Inf, ...)
```

## Arguments

- x:

  vector to be evaluated.

- spec:

  base name of a probability distribution (e.g., `"norm"`).

- a:

  lower bound of truncation, defaults to negative infinity.

- b:

  upper bound of truncation, defaults to infinity.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- ...:

  additional arguments for the probability distribution specified in
  `spec`.

- p:

  vector of probabilities.

- n:

  number of random draws.

## Value

`dtrunc` gives the density, and other functions provide related
computations.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `dtrunc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `extrunc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `ptrunc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `qtrunc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rtrunc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `vartrunc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dtrunc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving extrunc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving ptrunc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving qtrunc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rtrunc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving vartrunc
} # }
```
