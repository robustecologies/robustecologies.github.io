# Detect Wada basins using multiple methods

Unified dispatcher that runs one or more Wada-detection algorithms on a
basin structure and returns a consolidated `wada_analysis` object with a
per-method classification, the underlying diagnostics, and a majority
vote across the methods that ran successfully.

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

  Matrix of basin assignments or `wada_basins` object.

- system_func:

  Optional function defining the dynamical system, required only by the
  saddle-straddle method.

- attractors:

  Optional list of attractor specifications.

- x_range:

  Optional range of x coordinates.

- y_range:

  Optional range of y coordinates.

- methods:

  Character vector specifying which methods to use. Options: `"grid"`,
  `"merging"`, `"straddle"`, `"all"`.

- grid_params:

  List of parameters forwarded to
  [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md).

- straddle_params:

  List of parameters forwarded to
  [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md).

- verbose:

  Logical. Show progress and diagnostics.

## Value

A list of class `wada_analysis` containing the results from each method
invoked, the consensus classification, and the diagnostic metrics that
drove it.

## Details

This wrapper exists so that downstream code can run the three published
methods on the same basin object without juggling their differing input
contracts. The grid method (Daza et al. 2015) uses iterative box
refinement; the merging method (Daza et al. 2018) compares slim
boundaries via the Hausdorff distance; the saddle-straddle method
(Battelino, Grebogi, Ott & Yorke 1988; Wagemakers et al. 2020) tracks
the chaotic saddle by bisection. The merging and grid methods accept a
bare `wada_basins` object; the saddle-straddle method additionally
requires the system function and the attractor list because it
integrates trajectories on demand. The consensus rule is a strict
majority across the methods that returned a non-NA classification; if
fewer than two methods agree, the consensus is reported as inconclusive
rather than forced.

## References

Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
Testing for basins of Wada. *Scientific Reports*, 5, 16579.
[doi:10.1038/srep16579](https://doi.org/10.1038/srep16579)

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining when
a basin is Wada: The merging method. *Scientific Reports*, 8, 9954.
[doi:10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2020). The
saddle-straddle method to test for Wada basins. *Communications in
Nonlinear Science and Numerical Simulation*, 84, 105167.
[doi:10.1016/j.cnsns.2020.105167](https://doi.org/10.1016/j.cnsns.2020.105167)

Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., & Sanjuan, M.
A. F. (2016). Basin entropy: a new tool to analyze uncertainty in
dynamical systems. *Scientific Reports*, 6, 31416.
[doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

## See also

[`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md),
[`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md),
[`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md),
[`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md),
[`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
[`print.wada_analysis()`](https://robustecologies.github.io/wadaR/reference/print.wada_analysis.md),
[`summary.wada_analysis()`](https://robustecologies.github.io/wadaR/reference/summary.wada_analysis.md),
[`plot.wada_analysis()`](https://robustecologies.github.io/wadaR/reference/plot.wada_analysis.md).

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
print(result)
summary(result)
plot(result)
} # }
```
