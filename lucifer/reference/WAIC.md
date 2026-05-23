# Widely applicable information criterion

Calculates the Widely Applicable Information Criterion (WAIC), also
known as the Watanabe-Akaike information criterion, of Watanabe (2010).

S3 method: apply `WAIC()` to objects of class `default`.

S3 method: apply `WAIC()` to objects of class `demonoid`.

S3 method: apply `WAIC()` to objects of class `vb`.

S3 method: apply `WAIC()` to objects of class `pmc`.

S3 method: apply `WAIC()` to objects of class `sde_fit`.

## Usage

``` r
WAIC(x, ...)

# Default S3 method
WAIC(x, ...)

# S3 method for class 'demonoid'
WAIC(x, ...)

# S3 method for class 'vb'
WAIC(x, ...)

# S3 method for class 'pmc'
WAIC(x, ...)

# S3 method for class 'sde_fit'
WAIC(x, ...)
```

## Arguments

- x:

  See Details.

- ...:

  See Details.

## Value

An object of class `lucifer_waic`, which is a list with components:

- elpd_waic:

  The expected log pointwise predictive density.

- p_waic:

  The effective number of parameters (pWAIC2).

- waic:

  The WAIC value, \\-2 \\ \mathrm{elpd\\waic}\\.

- se_elpd_waic:

  Standard error of `elpd_waic`.

- se_waic:

  Standard error of `waic`.

- pointwise:

  A data.frame with per-observation `elpd_waic`, `p_waic`, and `waic`.

- WAIC:

  Legacy field: \\-2(lppd - pWAIC)\\.

- lppd:

  Legacy field: the log pointwise predictive density.

- pWAIC:

  Legacy field: effective number of parameters (pWAIC2).

- pWAIC1:

  Legacy field: alternate effective parameter count.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

WAIC is an extension of the Akaike Information Criterion (AIC) that is
more fully Bayesian than the Deviance Information Criterion (DIC).

Like DIC, WAIC estimates the effective number of parameters to adjust
for overfitting. Two adjustments have been proposed. pWAIC1 is similar
to pD in the original DIC. In contrast, pWAIC2 is calculated with
variance more similarly to pV, which Gelman proposed for DIC. Gelman et
al. (2014, p.174) recommends pWAIC2 because its results are closer in
practice to the results of leave-one-out cross-validation (LOO-CV).
pWAIC is considered an approximation to the number of unconstrained and
uninformed parameters, where a parameter counts as 1 when estimated
without constraint or any prior information, 0 if fully constrained or
all information comes from the prior distribution, or an intermediate
number if both the data and prior are informative, which is usually the
case.

Gelman et al. (2014, p. 174) scale the WAIC of Watanabe (2010) by a
factor of 2 so that it is comparable to AIC and DIC. WAIC is then
reported as \\-2(lppd - pWAIC)\\. Gelman et al. (2014) prefer WAIC to
AIC or DIC when feasible, which is less often than AIC or DIC. The
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function requires the model specification function to return the
model-level deviance, which is \\-2(LL)\\, where \\LL\\ is the sum of
the record-level log-likelihood. Therefore, if the user desires to
calculate WAIC, then the record-level log-likelihood must be monitored.

Implementation of `WAIC.default`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `WAIC.demonoid`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `WAIC.vb`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `WAIC.pmc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `WAIC.sde_fit`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and
Rubin, D.B. (2014). "Bayesian Data Analysis, 3rd ed.". CRC Press: Boca
Raton, FL.

Watanabe, S. (2010). "Asymptotic Equivalence of Bayes Cross Validation
and Widely Applicable Information Criterion in Singular Learning
Theory". *Journal of Machine Learning Research*, 11, p. 3571–3594.

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`loo_compare`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

[`LOO.default`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`as.demonoid.default`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LOO.demonoid`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik.demonoid`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md),
[`summary.demonoid`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.md).

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LOO.vb`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik.vb`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md),
[`summary.vb`](https://robustecologies.github.io/lucifer/reference/summary.vb.md).

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`LOO.pmc`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik.pmc`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.pmc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`print.pmc`](https://robustecologies.github.io/lucifer/reference/print.pmc.md),
[`to_mcmc_list.pmc`](https://robustecologies.github.io/lucifer/reference/to_mcmc_list.md).

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md),
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
N <- 10
S <- 1000
LL <- t(rmvn(S, -70 + rnorm(N),
     as.positive.definite(matrix(rnorm(N * N), N, N))))
result <- WAIC(LL)
print(result)
plot(result)
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving WAIC.default
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving WAIC.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving WAIC.vb
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving WAIC.pmc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving WAIC.sde_fit
} # }
```
