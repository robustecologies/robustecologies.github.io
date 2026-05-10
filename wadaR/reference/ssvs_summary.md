# Summarize spike-and-slab variable selection results

Extracts posterior inclusion probabilities and conditional coefficient
distributions from a fitted model that uses spike-and-slab (SSVS)
priors. Returns an S3 object with `print`, `summary`, and `plot` methods
providing 7 specialized visualization types.

## Usage

``` r
ssvs_summary(
  fit,
  indicators = NULL,
  coefficients = NULL,
  spike_var = 0.01,
  slab_var = NULL,
  prior_prob = 0.5,
  ground_truth = NULL,
  labels = NULL,
  dims = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  or a posterior matrix with named columns.

- indicators:

  Character vector of parameter name patterns (matched via `grep`) or
  integer vector of column indices identifying the binary inclusion
  indicators \\\gamma_j \in \\0, 1\\\\. If `NULL`, binary indicators are
  auto-detected via
  [`detect_indicators`](https://robustecologies.github.io/lucifer/reference/detect_indicators.md).

- coefficients:

  Character vector of parameter name patterns or integer vector of
  column indices identifying the regression coefficients whose inclusion
  is controlled by the indicators. The correspondence between
  `indicators` and `coefficients` is positional: the first indicator
  controls the first coefficient, and so on. If `NULL`, no conditional
  coefficient analysis is performed.

- spike_var:

  Numeric. Variance of the spike (null) component. Default is `0.01`.

- slab_var:

  Numeric or `NULL`. Variance of the slab (active) component. If `NULL`,
  the slab variance is not plotted in the spike-slab density overlay.

- prior_prob:

  Numeric. Prior probability of inclusion \\P(\gamma_j = 1)\\ under the
  Bernoulli prior. Default is `0.5`.

- ground_truth:

  Named numeric vector of true inclusion indicators (0 or 1). Names must
  match the indicator parameter names. Used to annotate plots with the
  ground truth status.

- labels:

  Character vector of human-readable labels for the
  indicators/coefficients. If `NULL`, parameter names are used.

- dims:

  Integer vector giving array dimensions for matrix indicators (e.g.,
  `c(J, J)` for a VAR adjacency matrix or `c(J, J, P)` for a multi-lag
  VAR). When supplied, the heatmap and network plot types reshape the
  indicator vector into this array.

## Value

An object of class `ssvs_summary` with components:

- pip:

  Named numeric vector of posterior inclusion probabilities.

- posterior:

  The full posterior matrix.

- ind_cols:

  Integer indices of the indicator columns.

- coef_cols:

  Integer indices of the coefficient columns (or NULL).

- ind_names:

  Character names of the indicator parameters.

- coef_names:

  Character names of the coefficient parameters (or NULL).

- spike_var:

  Spike variance.

- slab_var:

  Slab variance (or NULL).

- prior_prob:

  Prior inclusion probability.

- ground_truth:

  Ground truth indicators (or NULL).

- labels:

  Display labels.

- dims:

  Array dimensions for matrix indicators (or NULL).

## Details

Spike-and-slab priors assign each coefficient a two-component mixture
prior controlled by a binary indicator \\\gamma_j\\:

\$\$(\beta_j \mid \gamma_j) \sim (1 - \gamma_j)\\\mathcal{N}(0, \nu_0) +
\gamma_j\\\mathcal{N}(0, \nu_1)\$\$

where \\\nu_0 \ll \nu_1\\. The spike component \\\mathcal{N}(0, \nu_0)\\
concentrates mass near zero when the variable is excluded; the slab
component \\\mathcal{N}(0, \nu_1)\\ allows the coefficient to take
non-negligible values when included. The posterior inclusion probability
\\\text{PIP}\_j = P(\gamma_j = 1 \mid y)\\ is the primary output for
variable selection, with the median probability model retaining
variables where \\\text{PIP}\_j \geq 0.5\\ (Barbieri and Berger, 2004).

## References

George, E.I. and McCulloch, R.E. (1993). Variable selection via Gibbs
sampling. *Journal of the American Statistical Association*, 88(423),
881-889.
[doi:10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353)

Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model
selection. *Annals of Statistics*, 32(3), 870-897.
[doi:10.1214/009053604000000238](https://doi.org/10.1214/009053604000000238)

## See also

[`plot.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/plot.ssvs_summary.md),
[`print.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/print.ssvs_summary.md),
[`summary.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/summary.ssvs_summary.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit an SSVS model
fit <- lucifer(Model, Data, IV,
               Iterations = 10000, Algorithm = "Gibbs")

# Summarize SSVS results
ss <- ssvs_summary(fit,
                   indicators = "gamma",
                   coefficients = "beta",
                   spike_var = 0.01, slab_var = 100,
                   prior_prob = 0.5,
                   ground_truth = c("gamma[1]" = 1, "gamma[2]" = 0))

# Inspect
print(ss)
summary(ss)

# Visualize
plot(ss)
plot(ss, type = "spike_slab")
plot(ss, type = "heatmap", dims = c(3, 3))
} # }
```
