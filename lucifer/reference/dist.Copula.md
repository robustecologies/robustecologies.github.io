# Bivariate copula densities and random generation

Density and random generation for Clayton, Gumbel, and Frank bivariate
copulas.

See Details.

See Details.

S3 method: apply `dgumbel()` to objects of class `copula`.

S3 method: apply `rgumbel()` to objects of class `copula`.

See Details.

See Details.

## Usage

``` r
dclayton(u, v, theta, log = FALSE)

rclayton(n, theta = 1)

dgumbel.copula(u, v, theta, log = FALSE)

rgumbel.copula(n, theta = 1)

dfrank(u, v, theta, log = FALSE)

rfrank(n, theta = 1)
```

## Arguments

- u, v:

  vectors of pseudo-observations in \\(0, 1)\\.

- theta:

  vector of copula dependence parameters. For Clayton, \\\theta \> 0\\;
  for Gumbel, \\\theta \geq 1\\; for Frank, \\\theta \neq 0\\.

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

- n:

  number of observations, which must be a positive integer of length 1.

## Value

`dclayton`, `dgumbel.copula`, and `dfrank` give the copula density.
`rclayton`, `rgumbel.copula`, and `rfrank` return an \\n \times 2\\
matrix of pseudo-observations.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Copulas separate dependence structure from marginal distributions. For
bivariate data \\(X_1, X_2)\\ with marginal CDFs \\F_1, F_2\\, the
copula \\C\\ satisfies \\P(X_1 \le x_1, X_2 \le x_2) = C(F_1(x_1),
F_2(x_2))\\. The copula density \\c(u,v) = \partial^2 C / \partial u
\partial v\\ is provided for each family.

The **Clayton** copula captures lower tail dependence with generator
\\C(u,v) = (u^{-\theta} + v^{-\theta} - 1)^{-1/\theta}\\. The **Gumbel**
copula captures upper tail dependence with generator \\C(u,v) =
\exp(-((-\log u)^\theta + (-\log v)^\theta)^{1/\theta})\\. The **Frank**
copula is symmetric (no tail dependence) with generator involving the
Debye function.

Random generation uses conditional inversion for Clayton and Frank, and
the Marshall-Olkin stable variate method for Gumbel.

Implementation of `dclayton`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `rclayton`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `dgumbel.copula`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `rgumbel.copula`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `dfrank`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rfrank`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Nelsen, R.B. (2006). *An Introduction to Copulas*, 2nd ed. Springer.
ISBN 978-0-387-28659-4.

Joe, H. (2015). *Dependence Modeling with Copulas*. Chapman and
Hall/CRC. ISBN 978-1-4665-8323-8.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
u <- runif(100)
v <- runif(100)
d <- dclayton(u, v, theta = 2)
# Generate from Gumbel copula
uv <- rgumbel.copula(500, theta = 3)
plot(uv, pch = ".", main = "Gumbel copula samples")
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dclayton
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rclayton
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dgumbel.copula
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rgumbel.copula
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dfrank
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rfrank
} # }
```
