# Von Mises distribution

Density, distribution function, and random generation for the von Mises
distribution on the circle.

See Details.

See Details.

See Details.

## Usage

``` r
dvonmises(x, mu = 0, kappa = 1, log = FALSE)

pvonmises(q, mu = 0, kappa = 1)

rvonmises(n, mu = 0, kappa = 1)
```

## Arguments

- x, q:

  vector of angles (in radians).

- mu:

  vector of mean directions (in radians).

- kappa:

  vector of concentration parameters, which must be non-negative.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dvonmises` gives the density, `pvonmises` gives the distribution
function (integrated from \\-\pi\\ to `x`), and `rvonmises` generates
random deviates wrapped to \\\[-\pi, \pi)\\.

See Details.

See Details.

See Details.

## Details

The von Mises distribution is the circular analogue of the normal
distribution on the real line. Its density on the circle \\\[-\pi,
\pi)\\ is \$\$f(x; \mu, \kappa) = \frac{\exp(\kappa \cos(x - \mu))}{2\pi
I_0(\kappa)}\$\$ where \\I_0(\kappa)\\ is the modified Bessel function
of the first kind and order zero. When \\\kappa = 0\\ the distribution
is uniform on the circle; as \\\kappa \to \infty\\ it concentrates at
\\\mu\\. Random generation uses the Best-Fisher algorithm.

Implementation of `dvonmises`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `pvonmises`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rvonmises`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Mardia, K.V. and Jupp, P.E. (2000). *Directional Statistics*. John Wiley
& Sons. ISBN 978-0-471-95333-3.

Best, D.J. and Fisher, N.I. (1979). "Efficient Simulation of the von
Mises Distribution". *Journal of the Royal Statistical Society C*,
28(2), p. 152–157.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- seq(-pi, pi, length.out = 200)
d <- dvonmises(x, mu = 0, kappa = 2)
plot(x, d, type = "l", xlab = expression(theta), ylab = "Density")
r <- rvonmises(500, mu = 0, kappa = 5)
hist(r, breaks = 30, freq = FALSE)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dvonmises
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving pvonmises
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rvonmises
} # }
```
