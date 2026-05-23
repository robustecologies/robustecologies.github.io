# Plot method for an ensemble_sim object

Visualises ensemble simulation results with six complementary plot
types: `"fan"` (default) shows quantile ribbons of the trajectories per
state; `"spaghetti"` overlays individual sample paths; `"terminal"`
shows the empirical terminal distribution; `"mean_sd"` decomposes into
mean and standard deviation panels; `"extinction"` plots the cumulative
extinction fraction over time; `"convergence"` tracks the running mean
as replicates accumulate.

## Usage

``` r
# S3 method for class 'ensemble_sim'
plot(
  x,
  type = c("fan", "spaghetti", "terminal", "mean_sd", "extinction", "convergence"),
  state = NULL,
  max_traces = 50,
  title = NULL,
  ...
)
```

## Arguments

- x:

  An `ensemble_sim` object

- type:

  Plot type: `"fan"` (default), `"spaghetti"`, `"terminal"`,
  `"mean_sd"`, `"extinction"`, or `"convergence"`.

- state:

  Which state variable(s) to plot. Defaults to all states.

- max_traces:

  Maximum number of individual trajectories to overlay in spaghetti
  plots (default: 50).

- title:

  Optional plot title. If `NULL` (the default), a descriptive title is
  generated automatically.

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
