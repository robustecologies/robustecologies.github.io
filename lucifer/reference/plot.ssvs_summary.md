# Plot spike-and-slab variable selection results

Produces 7 specialized visualizations for spike-and-slab (SSVS) variable
selection models:

- `"pip"`:

  Posterior inclusion probability bar chart with prior reference line
  and decision threshold (default).

- `"spike_slab"`:

  Spike and slab prior densities overlaid with the marginal posterior
  density for each coefficient.

- `"conditional"`:

  Posterior density of each coefficient split by inclusion status
  (\\\gamma = 0\\ vs \\\gamma = 1\\).

- `"heatmap"`:

  Matrix heatmap of inclusion probabilities, suitable for VAR/VARMA
  models where indicators have array structure. Requires `dims` in the
  `ssvs_summary` object.

- `"trajectory"`:

  Running posterior inclusion probability across MCMC iterations,
  showing convergence of the selection decision.

- `"network"`:

  Directed network graph where edge opacity and width are proportional
  to \\P(\Gamma\_{ij} = 1)\\. Requires `dims` with length 2 or 3.

- `"fdr"`:

  Bayesian false discovery rate (FDR) and false omission rate (FOR)
  curves as functions of the inclusion threshold. Requires
  `ground_truth`.

## Usage

``` r
# S3 method for class 'ssvs_summary'
plot(
  x,
  type = c("pip", "spike_slab", "conditional", "heatmap", "trajectory", "network", "fdr"),
  threshold = 0.5,
  col = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `ssvs_summary`.

- type:

  Character. Plot type (see above). Default is `"pip"`.

- threshold:

  Numeric. Inclusion probability threshold for the median model. Default
  is `0.5`.

- col:

  Optional character vector of colors.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the ggplot object(s).

## Details

Implementation of `plot.ssvs_summary`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`print.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/print.ssvs_summary.md),
[`summary.ssvs_summary`](https://robustecologies.github.io/lucifer/reference/summary.ssvs_summary.md).

## Examples

``` r
if (FALSE) { # \dontrun{
ss <- ssvs_summary(fit, indicators = "gamma", coefficients = "beta",
                   spike_var = 0.01, slab_var = 100)
plot(ss)                          # PIP bar chart
plot(ss, type = "spike_slab")     # Spike-slab prior vs posterior
plot(ss, type = "conditional")    # Conditional densities
plot(ss, type = "trajectory")     # Running PIP across iterations
plot(ss, type = "heatmap")        # Matrix heatmap
plot(ss, type = "network")        # Directed network
plot(ss, type = "fdr")            # FDR/FOR curves
} # }
```
