# Plot variable importance

This may be used to plot variable importance with BPIC, predictive
concordance, a discrepancy statistic, or the L-criterion regarding an
object of class `importance`.

## Usage

``` r
# S3 method for class 'importance'
plot(x, Style = "BPIC", col = NULL, ...)
```

## Arguments

- x:

  This required argument is an object of class `importance`.

- Style:

  When `Style="BPIC"`, BPIC is shown, and `BPIC` is the default.
  Otherwise, predictive concordance is plotted when
  `Style="Concordance"`, a discrepancy statistic is plotted when
  `Style="Discrep"`, or the L-criterion is plotted when
  `Style="L-criterion"`.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

The x-axis is either BPIC (Ando, 2007), predictive concordance (Gelfand,
1996), a discrepancy statistic (Gelman et al., 1996), or the L-criterion
(Laud and Ibrahim, 1995) of the
[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md)
function (depending on the `Style` argument), and variables are on the
y-axis. A more important variable is associated with a dot that is
plotted farther to the right. For more information on variable
importance, see the
[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md)
function.

## References

Ando, T. (2007). "Bayesian Predictive Information Criterion for the
Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
*Biometrika*, 94(2), p. 443–458.

Gelfand, A. (1996). "Model Determination Using Sampling Based Methods".
In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in Markov
Chain Monte Carlo in Practice. Chapman and Hall: Boca Raton, FL.

Gelman, A., Meng, X.L., and Stern H. (1996). "Posterior Predictive
Assessment of Model Fitness via Realized Discrepancies". *Statistica
Sinica*, 6, p. 733–807.

Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model Selection".
*Journal of the Royal Statistical Society*, B 57, p. 247–262.

## See also

[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.importance
} # }
```
