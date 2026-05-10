# Precision, standard deviation, and variance conversions

Utility functions that convert between precision, standard deviation,
and variance for scalars, vectors, and matrices. Capital letters refer
to matrix operations; lower case letters refer to scalars and vectors.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Usage

``` r
Cov2Prec(Cov)

Prec2Cov(Prec)

prec2sd(prec = 1)

prec2var(prec = 1)

sd2prec(sd = 1)

sd2var(sd = 1)

var2prec(var = 1)

var2sd(var = 1)
```

## Arguments

- Cov:

  A covariance matrix, \\\Sigma\\.

- Prec:

  A precision matrix, \\\Omega\\.

- prec:

  A precision scalar or vector, \\\tau\\.

- sd:

  A standard deviation scalar or vector, \\\sigma\\.

- var:

  A variance scalar or vector, \\\sigma^2\\.

## Value

`Cov2Prec` returns a precision matrix \\\Omega = \Sigma^{-1}\\.
`Prec2Cov` returns a covariance matrix \\\Sigma = \Omega^{-1}\\.
`prec2sd` returns \\\sigma = \sqrt{\tau^{-1}}\\. `prec2var` returns
\\\sigma^2 = \tau^{-1}\\. `sd2prec` returns \\\tau = \sigma^{-2}\\.
`sd2var` returns \\\sigma^2\\. `var2prec` returns \\\tau = 1/\sigma^2\\.
`var2sd` returns \\\sigma = \sqrt{\sigma^2}\\.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Bayesians often use precision rather than variance, where precision is
the inverse of the variance. The modern Bayesian use of precision
developed because it was more straightforward in a normal distribution
to estimate precision \\\tau\\ with a gamma distribution as a conjugate
prior, than to estimate \\\sigma^2\\ with an inverse-gamma distribution
as a conjugate prior. Today, conjugacy is usually considered to be
merely a convenience, and a non-conjugate half-Cauchy prior distribution
is recommended as a weakly informative prior distribution for scale
parameters.

For example, a linear regression may be represented equivalently as
\\\mathbf{y} \sim \mathcal{N}(\mu, \sigma^2)\\, or \\\mathbf{y} \sim
\mathcal{N}(\mu, \tau^{-1})\\, where \\\sigma^2\\ is the variance and
\\\tau\\ is the precision.

Implementation of `Cov2Prec`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `Prec2Cov`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `prec2sd`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `prec2var`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `sd2prec`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `sd2var`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `var2prec`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `var2sd`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`Cov2Cor`](https://robustecologies.github.io/lucifer/reference/Matrices.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
Cov2Prec(matrix(c(1, 0.1, 0.1, 1), 2, 2))
Prec2Cov(matrix(c(1, 0.1, 0.1, 1), 2, 2))
prec2sd(0.5)
prec2var(0.5)
sd2prec(1.4142)
sd2var(1.4142)
var2prec(2)
var2sd(2)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Cov2Prec
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Prec2Cov
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving prec2sd
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving prec2var
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving sd2prec
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving sd2var
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving var2prec
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving var2sd
} # }
```
