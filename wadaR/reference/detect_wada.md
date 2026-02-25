# Detect Wada basins using multiple methods

A unified interface to test for Wada basins using one or more detection
methods. Returns a comprehensive report of the results.

## Usage

``` r
detect_wada(
  basins,
  system_func = NULL,
  attractors = NULL,
  x_range = NULL,
  y_range = NULL,
  methods = c("merging", "grid"),
  grid_params = list(),
  straddle_params = list(),
  verbose = TRUE
)
```

## Arguments

- basins:

  Matrix of basin assignments or wada_basins object.

- system_func:

  Optional function defining the dynamical system.

- attractors:

  Optional list of attractor specifications.

- x_range:

  Optional range of x coordinates.

- y_range:

  Optional range of y coordinates.

- methods:

  Character vector specifying which methods to use. Options: "grid",
  "merging", "straddle", "all".

- grid_params:

  List of parameters for the grid method.

- straddle_params:

  List of parameters for the saddle-straddle method.

- verbose:

  Logical. Show progress and diagnostics.

## Value

A list of class "wada_analysis" containing results from all methods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using the forced damped pendulum
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum$system, c(-pi, pi), c(-3, 3),
                         resolution = 200, attractors = pendulum$attractors)

# Test with all methods
result <- detect_wada(basins, system_func = pendulum$system,
                      attractors = pendulum$attractors, methods = "all")
} # }
```
