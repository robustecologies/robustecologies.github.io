# Consort with Lucifer

Evaluates a fitted Bayesian model against family-specific convergence
criteria and returns a structured diagnostic report with actionable
suggestions. Supports all inference families: MCMC, variational Bayes,
Laplace approximation, iterative quadrature, PMC, and SMC.

## Usage

``` r
Consort(object = NULL, verbose = TRUE, ...)
```

## Arguments

- object:

  A fitted model object of class `demonoid`, `vb`, `laplace`,
  `iterquad`, `pmc`, or `smc`.

- verbose:

  Logical. If `TRUE` (default), prints the diagnostic report to the
  console via `print.consort`.

- ...:

  Additional arguments. Currently the only recognised name is
  `thresholds`, an optional list of MCMC convergence thresholds that
  overrides the defaults from `.consort_thresholds()`.

## Value

An object of class `consort` containing:

- family:

  Character. Inference family: "MCMC", "VB", "Laplace", "IQ", "PMC", or
  "SMC".

- algorithm:

  Character. Full algorithm name.

- appeased:

  Logical. Whether all convergence criteria are met.

- conditions:

  Data frame with columns `condition`, `status` ("pass" or "fail"), and
  `detail`.

- diagnostics:

  Family-specific list of computed diagnostic values.

- suggestion:

  List with `action`, `next_specs`, `code`, and `rationale`. NULL if
  appeased.

- fit_summary:

  List with `iterations`, `thinned_samples`, `minutes`, `parameters`,
  and `lp_final`.

- Call:

  The matched call.

## Details

Consort dispatches to a family-specific diagnostic engine that evaluates
the fit against calibrated convergence criteria. For MCMC fits, these
include acceptance rate within the algorithm's recommended range, MCSE
below 6.27% of marginal posterior SD, bulk ESS at least 400, tail ESS at
least 200, split R-hat below 1.01 (or BMK stationarity for single
chains), zero divergent transitions, and the non-adaptive algorithm
requirement. For non-MCMC families, convergence, stability, and
method-specific diagnostics are checked.

When all conditions pass, the fit is declared appeased. Otherwise,
Consort generates a structured suggestion containing a recommended
action, ready-to-use specification list, human-readable code string, and
rationale. The suggestion is informed by the algorithm registry and
follows established escalation paths (e.g., HMC to NUTS, AMM to AMWG,
DEMC to AIES).

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-Normalization, Folding, and Localization: An Improved Rhat
for Assessing Convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

Flegal, J.M., Haran, M. and Jones, G.L. (2008). Markov chain Monte
Carlo: Can We Trust the Third Significant Figure? *Statistical Science*,
23, 250-260. [doi:10.1214/08-STS257](https://doi.org/10.1214/08-STS257)

## See also

[`is.appeased`](https://robustecologies.github.io/lucifer/reference/is.appeased.md),
[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit MCMC model
fit <- lucifer(Model, Data, Initial.Values,
  Iterations = 10000, Status = 1000, Thinning = 10,
  Algorithm = "NUTS",
  Specs = list(A = 500, delta = 0.65, epsilon = 0.01, Lmax = 20))

# Get diagnostic report
cx <- Consort(fit)
print(cx)
summary(cx)
plot(cx)

# If not appeased, use suggestion
if (!cx$appeased) {
  cat(cx$suggestion$code)
}

# Works with VB fits too
vb_fit <- VariationalBayes(Model, Data, parm = IV,
  Iterations = 1000, Method = "Pathfinder", sir = TRUE)
Consort(vb_fit)
} # }
```
