# Compute model weights via stacking or pseudo-BMA

Given two or more models with pointwise predictive criteria (LOO, WAIC,
K-fold, or LFO), computes optimal combination weights. Stacking (Yao et
al. 2018) maximizes the combined log predictive density and outperforms
naive model selection when models are misspecified. Pseudo-BMA and
pseudo-BMA+ provide faster alternatives based on information criterion
differences.

## Usage

``` r
stacking_weights(
  ...,
  method = c("stacking", "pseudobma", "pseudobma+"),
  criterion = NULL
)
```

## Arguments

- ...:

  Two or more objects of class `lucifer_loo`, `lucifer_waic`,
  `lucifer_kfold`, or `lucifer_lfo`, each containing a `pointwise` data
  frame with an `elpd_*` column.

- method:

  One of `"stacking"` (default), `"pseudobma"`, or `"pseudobma+"`. See
  Details.

- criterion:

  Optional column name override for the elpd column. If `NULL`,
  auto-detected.

## Value

An object of class `lucifer_stacking_weights`, which is a named numeric
vector of model weights summing to 1, with additional attributes storing
the method used.

## Details

**Stacking** (Yao et al. 2018) finds weights \\w\\ that maximize
\\\sum\_{i=1}^N \log \sum\_{k=1}^K w_k \exp(\mathrm{ elpd}\_{ik})\\
subject to \\w_k \ge 0\\ and \\\sum w_k = 1\\. This is solved via
[`stats::constrOptim()`](https://rdrr.io/r/stats/constrOptim.html).

**Pseudo-BMA** assigns weights proportional to
\\\exp(\mathrm{elpd}\_k)\\ where \\\mathrm{elpd}\_k = \sum_i
\mathrm{elpd}\_{ik}\\.

**Pseudo-BMA+** incorporates Bayesian bootstrap uncertainty: bootstrap
weights are drawn from an exponential distribution and applied to the
pointwise elpd differences, yielding a distribution of model rankings.

## References

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). "Using
stacking to average Bayesian predictive distributions (with
discussion)." *Bayesian Analysis*, 13(3), 917–1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

## See also

[`loo_compare`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`Kfold`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO`](https://robustecologies.github.io/lucifer/reference/LFO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
loo1 <- LOO(log_lik_model1)
loo2 <- LOO(log_lik_model2)
loo3 <- LOO(log_lik_model3)

# Stacking weights (recommended)
w <- stacking_weights(m1 = loo1, m2 = loo2, m3 = loo3)
print(w)

# Pseudo-BMA
w2 <- stacking_weights(m1 = loo1, m2 = loo2, method = "pseudobma")
} # }
```
