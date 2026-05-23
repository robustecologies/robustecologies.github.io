# Plot method for product_space objects

Produces diagnostic and summary plots for a product space model
selection analysis. The default type (`"all"`) displays a combined panel
with all five visualization types.

## Usage

``` r
# S3 method for class 'product_space'
plot(x, type = "all", col = NULL, ...)
```

## Arguments

- x:

  An object of class `product_space`.

- type:

  Character: `"all"` (default), `"probabilities"`, `"trace"`,
  `"transition"`, `"bayes_factors"`, `"posteriors"`, or `"complexity"`.

- col:

  Optional character vector of hex color strings. Overrides the default
  palette.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns the plot object(s).

## Details

The available plot types are:

- `"all"`:

  Combined panel: top row shows probabilities, transition matrix, Bayes
  factors, and posteriors in a 2x2 grid; bottom row shows the full-width
  model indicator trace.

- `"probabilities"`:

  Bar chart of corrected posterior model probabilities with a horizontal
  reference line at 1/K (uniform). The preferred model is highlighted;
  model space complexity metrics are annotated in the subtitle.

- `"trace"`:

  Trace plot of the model indicator M across MCMC iterations, showing
  model switching patterns.

- `"transition"`:

  Heatmap of the K x K transition matrix p(M\[t+1\] = j \| M\[t\] = i).

- `"bayes_factors"`:

  Forest plot of log Bayes factors relative to the best model, with
  Jeffreys scale reference regions and an interpretive legend.

- `"posteriors"`:

  Density plots of per-model posterior distributions for each parameter,
  faceted by model.

- `"complexity"`:

  Horizontal bullet chart of four model space complexity metrics
  (Shannon normalized entropy, effective K, Simpson diversity, Gini
  effective K), all rescaled to \\\[0,1\]\\. Background zones indicate
  decisive, moderate, and ambiguous model selection regimes. Per-metric
  raw values and theoretical ranges are annotated.

## See also

[`ProductSpace`](https://robustecologies.github.io/lucifer/reference/ProductSpace.md),
[`print.product_space`](https://robustecologies.github.io/lucifer/reference/print.product_space.md),
[`summary.product_space`](https://robustecologies.github.io/lucifer/reference/summary.product_space.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.product_space
} # }
```
