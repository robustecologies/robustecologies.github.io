# Math utility functions

Utility functions for mathematical operations including Gauss-Hermite
quadrature, Hermite polynomials, log-addition, and partial derivatives.

See Details.

See Details.

See Details.

## Usage

``` r
GaussHermiteQuadRule(N)

Hermite(x, N, prob = TRUE)

logadd(x, add = TRUE)

partial(Model, parm, Data, Interval = 1e-06, Method = "simple")
```

## Arguments

- N:

  A positive integer indicating the number of nodes.

- x:

  A numeric vector.

- prob:

  Logical. If `TRUE` (default), uses the probabilist's kernel for the
  Hermite polynomial; otherwise the physicist's kernel.

- add:

  Logical. If `TRUE` (default), \\\log(x+y)\\ is performed; otherwise
  \\\log(x-y)\\.

- Model:

  A model specification function. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- parm:

  A vector of parameters.

- Data:

  A list of data. See
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Interval:

  The interval of numeric differencing.

- Method:

  Either `"simple"` for finite-differencing (default) or `"Richardson"`
  for Richardson extrapolation, which is more accurate but slower.

## Value

`GaussHermiteQuadRule` returns a list with `nodes` and `weights`.
`Hermite` returns a numeric vector. `logadd` returns a scalar. `partial`
returns a vector of partial derivatives.

See Details.

See Details.

See Details.

## Details

`GaussHermiteQuadRule` returns nodes and weights for Gauss-Hermite
quadrature, used by
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md).
`Hermite` returns the probabilist's or physicist's Hermite polynomial.
`logadd` performs addition or subtraction on the log scale, useful for
preventing numerical underflow. `partial` estimates partial derivatives
(the gradient) of the model with respect to its parameters via finite
differencing or Richardson extrapolation.

Implementation of `Hermite`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `logadd`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `partial`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving GaussHermiteQuadRule
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Hermite
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving logadd
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving partial
} # }
```
