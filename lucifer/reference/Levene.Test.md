# Levene's test

A Bayesian form of Levene's test (Levene, 1960) for equality of
variances across groups. This function estimates the test statistic
\\W\\ from observed residuals (\\W^{obs}\\) and from residuals
replicated under a homoskedastic process (\\W^{rep}\\), for each
posterior sample.

## Usage

``` r
Levene.Test(x, Method = "U", G = NULL, Data = NULL)
```

## Arguments

- x:

  an object of class `demonoid.ppc`, `iterquad.ppc`, `laplace.ppc`,
  `pmc.ppc`, or `vb.ppc`.

- Method:

  character; `"U"` (default) for univariate DV `y`, `"C"` to apply the
  test to each column of `Y`, or `"R"` to apply the test to each row of
  `Y`.

- G:

  a grouping indicator with the same dimensions as the DV, or `NULL`
  (default) in which case 4 groups are created by partitioning records
  into quarters.

- Data:

  a list of data required when `Method = "C"` or `Method = "R"`. The DV
  must be named `Y`.

## Value

A plot (or series of plots for multivariate `Y`) showing the density of
\\W^{obs}\\ and \\W^{rep}\\, and a numeric vector of the probability
\\p(W^{obs} \> W^{rep})\\ for each test conducted.

## Details

Levene's test assesses the probability of equality of residual variances
across groups. An advantage over other tests of homoskedasticity is that
Levene's test does not require normality of the residuals. The function
reports \\p(W^{obs} \> W^{rep})\\, where values below 0.025 or above
0.975 indicate heteroskedastic variance.

## References

Levene, H. (1960). "Robust Tests for Equality of Variances". In I.
Olkins et al. (Eds.), *Contributions to Probability and Statistics*, p.
278–292. Stanford University Press.

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# First, update the model with lucifer (or similar), then:
# Pred <- predict(Fit, Model, MyData)
# Levene.Test(Pred)
} # }
```
