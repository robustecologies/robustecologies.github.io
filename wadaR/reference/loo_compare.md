# Compare models using pointwise predictive criteria

Compares two or more models fitted to the same data using their
pointwise expected log predictive density (elpd) estimates. Models are
ranked from best to worst, and the standard error of the pairwise
differences accounts for the correlation between pointwise estimates
evaluated on the same observations (Vehtari et al. 2017, Section 5).

## Usage

``` r
loo_compare(..., criterion = NULL)
```

## Arguments

- ...:

  Two or more objects of class `lucifer_loo`, `lucifer_waic`,
  `lucifer_kfold`, or `lucifer_lfo`. Each must contain a `pointwise`
  data frame with a column whose name starts with `elpd_`. All objects
  must have the same number of observations (rows in `pointwise`).

- criterion:

  Optional character string to override automatic detection of the elpd
  column name. If `NULL` (default), the function extracts the first
  column matching `"elpd_*"` from each object's pointwise data frame.

## Value

An object of class `lucifer_loo_comparison`, which is a data frame with
one row per model sorted by elpd (best first), and columns for elpd,
se_elpd, elpd_diff, and se_diff. The difference is always relative to
the best model (first row), so the first row has elpd_diff = 0 and
se_diff = 0.

## Details

The standard error of the elpd difference between two models is computed
as \\\mathrm{se\\diff} = \sqrt{N \\ \mathrm{Var} (\mathrm{elpd}\_{1,i} -
\mathrm{elpd}\_{2,i})}\\, which is always smaller than the naive
\\\sqrt{\mathrm{se}\_1^2 + \mathrm{se}\_2^2}\\ because the pointwise
estimates are positively correlated across models evaluated on the same
data.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). "Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC." *Statistics
and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). "Using
stacking to average Bayesian predictive distributions (with
discussion)." *Bayesian Analysis*, 13(3), 917–1007.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`Kfold`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`LFO`](https://robustecologies.github.io/lucifer/reference/LFO.md),
[`stacking_weights`](https://robustecologies.github.io/lucifer/reference/stacking_weights.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit two models and compute LOO for each
loo1 <- LOO(log_lik_model1)
loo2 <- LOO(log_lik_model2)

# Compare
comp <- loo_compare(model1 = loo1, model2 = loo2)
print(comp)
plot(comp)
} # }
```
