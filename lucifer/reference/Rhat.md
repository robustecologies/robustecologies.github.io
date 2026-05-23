# Split R-hat convergence diagnostic

Computes the rank-normalized split R-hat statistic of Vehtari, Gelman,
Simpson, Carpenter, and Burkner (2021). Each chain is split in half
before the standard R-hat is computed on the rank-normalized draws,
improving sensitivity to non-stationarity and bias that the classical
R-hat of Gelman and Rubin (1992) may miss.

## Usage

``` r
Rhat(x)
```

## Arguments

- x:

  a numeric vector (single chain, split internally), a matrix where rows
  are iterations and columns are chains (single parameter), a list of
  numeric vectors (one per chain), an object of class `demonoid` (single
  or multi-chain), or a list of `demonoid` objects.

## Value

For a single parameter: a scalar split R-hat value. For a
multi-parameter object (e.g., a `demonoid`): a named numeric vector with
one split R-hat per parameter.

## Details

The split R-hat addresses two weaknesses of the classical R-hat. First,
splitting each chain in half doubles the number of chains used for the
between-chain variance estimate, which makes the diagnostic sensitive to
non-stationarity within individual chains. Second, rank-normalizing the
draws before computing R-hat ensures that the statistic is well-defined
even for distributions with heavy tails or infinite variance.

The rank-normalization proceeds by replacing each draw with its
corresponding standard normal quantile:

\$\$z_r = \Phi^{-1}\\\left(\frac{r - 3/8}{S - 1/4}\right)\$\$

where \\r\\ is the rank of the draw among all \\S\\ draws across all
chains, and \\\Phi^{-1}\\ is the standard normal quantile function. The
fractional offset \\3/8\\ is the Blom (1958) adjustment.

After rank-normalization, the standard R-hat formula is applied to the
transformed, split chains:

\$\$\hat{R} = \sqrt{\frac{\hat{V}}{W}}\$\$

where \\W\\ is the mean within-chain variance and \\\hat{V} =
\frac{n-1}{n}W + \frac{B}{n}\\ with \\B\\ the between-chain variance.

A split R-hat below 1.01 is generally considered adequate for most
applications. Values above 1.05 strongly indicate non-convergence.

When `x` is a single chain (vector) or a single-chain `demonoid`, the
chain is split in half and R-hat is computed on the two halves. For
multi-chain input, each chain is split, yielding \\2M\\ half-chains. A
minimum of 100 draws per half-chain is recommended.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Burkner, P.-C.
(2021). "Rank-normalization, folding, and localization: An improved
R-hat for assessing convergence of MCMC (with discussion)". *Bayesian
Analysis*, 16(2), p. 667–718. doi:10.1214/20-BA1221.

Gelman, A. and Rubin, D.B. (1992). "Inference from Iterative Simulation
using Multiple Sequences". *Statistical Science*, 7, p. 457–511.

## See also

[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`ESS.bulk`](https://robustecologies.github.io/lucifer/reference/ESS.bulk.md),
[`ESS.tail`](https://robustecologies.github.io/lucifer/reference/ESS.tail.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single chain: split internally
x <- cumsum(rnorm(2000))
Rhat(x)

# Multiple chains as matrix columns
chains <- matrix(rnorm(4000), ncol = 4)
Rhat(chains)

# Multi-chain demonoid
fit <- lucifer(Model, MyData, IV,
    Iterations = 5000, Algorithm = "NUTS",
    Specs = list(A = 500, delta = 0.65, epsilon = 0.01, Lmax = 20),
    Chains = 4, CPUs = 4)
Rhat(fit)
} # }
```
