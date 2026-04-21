# Plot method for ensemble_sim objects

Visualizes ensemble simulation results. Six plot types are available via
the `type` argument.

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
