# Plot an arena object

Produces publication-quality visualizations from the universal
benchmarking comparison. Eight plot types are available, each
highlighting a different aspect of the cross-method comparison.

## Usage

``` r
# S3 method for class 'arena'
plot(x, type = "efficiency", Parms = NULL, ...)
```

## Arguments

- x:

  An object of class `arena`.

- type:

  Character. The type of plot to produce. One of `"efficiency"`
  (default), `"accuracy"`, `"pareto"`, `"convergence"`, `"radar"`,
  `"heatmap"`, `"forest"`, or `"pairwise"`.

- Parms:

  Integer vector or character vector of parameter indices or names for
  plots that display per-parameter results. Default `NULL` selects up to
  9 parameters.

- ...:

  Additional arguments passed to individual plot functions.

## Value

Invisibly returns the ggplot object(s).

## Details

The available plot types are:

- `"efficiency"`:

  Horizontal bar chart of methods ordered by ESS per second, color-coded
  by inference category. Converged methods are shown with solid bars;
  non-converged with dashed outlines.

- `"accuracy"`:

  Marginal posterior density overlays per parameter, with each method as
  a colored line and the reference as a thick black dashed line.

- `"pareto"`:

  Scatterplot of log10(minutes) vs. accuracy (Wasserstein distance to
  reference), with the Pareto frontier highlighted and method labels.

- `"convergence"`:

  Convergence trace plots per method (LP trace for MCMC, ELBO for VB,
  ESS history for SMC).

- `"radar"`:

  Spider chart with one polygon per method, axes representing normalized
  efficiency, accuracy, speed, reliability, and time.

- `"heatmap"`:

  Rank heatmap with methods as rows and metrics as columns, cell color
  by rank, annotated with rank numbers.

- `"forest"`:

  Forest plot showing posterior mean and 95% credible interval per
  method, for selected parameters.

- `"pairwise"`:

  Lower-triangle heatmap of pairwise Wasserstein distance between
  methods.

## See also

[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`print.arena`](https://robustecologies.github.io/lucifer/reference/print.arena.md),
[`summary.arena`](https://robustecologies.github.io/lucifer/reference/summary.arena.md)

## Examples

``` r
if (FALSE) { # \dontrun{
arena_result <- Arena(list(NUTS = fit1, VB = fit2, Laplace = fit3))

plot(arena_result)
plot(arena_result, type = "accuracy", Parms = 1:4)
plot(arena_result, type = "pareto")
plot(arena_result, type = "heatmap")
plot(arena_result, type = "forest", Parms = c("mu", "sigma"))
plot(arena_result, type = "pairwise")
} # }
```
