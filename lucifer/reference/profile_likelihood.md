# Profile likelihood curves from any lucifer fit

Computes profile log-likelihood curves for each parameter by fixing the
parameter of interest at a grid of values and optimizing over the
remaining nuisance parameters. The resulting curves provide confidence
intervals that account for parameter curvature and asymmetry, unlike the
symmetric Wald intervals from
[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md).

## Usage

``` r
profile_likelihood(
  object,
  Model,
  Data,
  parm = NULL,
  conf_level = 0.95,
  n_grid = 50L
)
```

## Arguments

- object:

  a fitted model object of class `laplace`, `iterquad`, `demonoid`, or
  `data_cloning`.

- Model:

  the Model function used for fitting.

- Data:

  the Data list used for fitting.

- parm:

  character vector of parameter names to profile. If `NULL` (default),
  all parameters are profiled.

- conf_level:

  confidence level for the profile CI, default 0.95.

- n_grid:

  number of grid points per parameter, default 50.

## Value

An object of class `profile_likelihood` with elements:

- profiles:

  named list of per-parameter profiles, each containing `grid`,
  `profile_ll`, `threshold`, and `ci`

- MLE:

  numeric vector of MLE estimates

- SE:

  numeric vector of standard errors

- conf_level:

  confidence level used

- parm.names:

  character vector of profiled parameter names

## Details

For each parameter \\\theta_j\\, the profile log-likelihood is

\$\$\ell_P(\theta_j) = \max\_{\theta\_{-j}} \ell(\theta_j,
\theta\_{-j})\$\$

where \\\theta\_{-j}\\ denotes the nuisance parameters. The profile
likelihood ratio \\R_P(\theta_j) = \ell_P(\theta_j) - \ell(\hat\theta)\\
is plotted against \\\theta_j\\, with a horizontal line at the
\\-\chi^2_1(\alpha)/2\\ threshold determining the profile confidence
interval \[\[1\]\](#ref1).

For univariate models (\\p = 1\\) no nuisance optimization is needed;
the profile equals the log-likelihood evaluated directly on the grid.

The MLE and standard errors are extracted from the fit object via
[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)
internals. The grid spans \\\hat\theta_j \pm 4 \cdot \mathrm{SE}\_j\\
with `n_grid` points. Nuisance optimization uses
`stats::optim(method = "BFGS")` starting at the full MLE.

Profile confidence intervals are more reliable than Wald intervals for
nonlinear models, boundary parameters, and small samples, because they
respect the actual curvature of the likelihood surface rather than
approximating it as quadratic \[\[2\]\](#ref2).

## References

Venzon, D. J. and Moolgavkar, S. H. (1988). A method for computing
profile-likelihood-based confidence intervals. \*Applied Statistics\*,
37(1), 87-94. DOI: [10.2307/2347496](https://doi.org/10.2307/2347496).

Meeker, W. Q. and Escobar, L. A. (1995). Teaching about approximate
confidence regions based on maximum likelihood estimation. \*The
American Statistician\*, 49(1), 48-53. DOI:
[10.1080/00031305.1995.10476112](https://doi.org/10.1080/00031305.1995.10476112).

## See also

[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md),
[`confint_compare`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- LaplaceApproximation(Model, Initial.Values, Data,
                             Iterations = 200, Method = "SPG",
                             CovEst = "Hessian", sir = FALSE)

# Profile all parameters
prof <- profile_likelihood(fit, Model, Data)
print(prof)
plot(prof)

# Profile specific parameters
prof_sub <- profile_likelihood(fit, Model, Data, parm = c("beta0", "beta1"))
plot(prof_sub)
} # }
```
