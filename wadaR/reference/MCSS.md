# Monte Carlo sample size

The `MCSS` function estimates the required number of posterior samples,
given the user-specified acceptable error `a`, posterior samples `x`,
and the observed variance (rather than asymptotic variance). Due to the
observed variance, this is a rough estimate.

## Usage

``` r
MCSS(x, a)
```

## Arguments

- x:

  a vector of posterior samples.

- a:

  a positive scalar of acceptable error for the mean of `x`. As
  acceptable error decreases, the required number of samples increases.

## Value

The estimated required number of posterior samples.

## Details

Implementation of `MCSS`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(1000)
MCSS(x, a=0.01)
} # }
```
