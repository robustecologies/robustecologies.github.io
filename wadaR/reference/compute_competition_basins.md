# Compute basins of attraction for multispecies competition

Computes the basins of attraction for a Huisman-Weissing multispecies
competition system by varying initial conditions of two species and
determining which outcome (winning coalition) results.

## Usage

``` r
compute_competition_basins(
  system,
  x_range = c(0, 2),
  y_range = c(0, 2),
  resolution = 200,
  t_max = NULL,
  dt = 0.01,
  extinction_threshold = 0.01,
  dominance_threshold = 1,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- system:

  A multispecies competition system object created by
  [`multispecies_competition`](https://robustecologies.github.io/wadaR/reference/multispecies_competition.md).

- x_range:

  Numeric vector of length 2. Range of initial abundance for the first
  projected species (default: N2).

- y_range:

  Numeric vector of length 2. Range of initial abundance for the second
  projected species (default: N4).

- resolution:

  Integer. Number of grid points per dimension.

- t_max:

  Numeric. Maximum simulation time (days). For 5-species scenarios,
  default is 2500. For 8-species, default is 2000.

- dt:

  Numeric. Integration time step (default 0.01 days).

- extinction_threshold:

  Numeric. Abundance below which a species is considered extinct
  (default 0.01).

- dominance_threshold:

  Numeric. For 8-species scenario, minimum abundance to be considered
  dominant (default 1.0).

- n_cores:

  Integer. Number of CPU cores to use for parallel computation. Default
  is `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`
  (leaving one core free for system tasks). Set to 1 to disable
  parallelization.

- verbose:

  Logical. Show progress information.

## Value

A list of class `"wada_basins"` containing:

- basins:

  Integer matrix of outcome assignments

- x_grid:

  Vector of x coordinates (N2 initial conditions)

- y_grid:

  Vector of y coordinates (N4 initial conditions)

- n_attractors:

  Number of distinct outcomes

- unclassified:

  Number of points with undetermined outcome

- x_range:

  Range of x values

- y_range:

  Range of y values

- resolution:

  Grid resolution

- scenario:

  The competition scenario ("5species" or "8species")

- attractors:

  List of outcome descriptions

## Details

Unlike physical dynamical systems where attractors are fixed points in
phase space, the "attractors" in competition systems are *outcomes*:
different stable coalitions of species that can emerge depending on
initial conditions.

**2D projection strategy:**

The Huisman-Weissing system has \\n\\ species dimensions plus \\k\\
resource dimensions. To visualize basins, we project to 2D by varying
the initial conditions of two key species (typically N2 and N4,
following the original paper) while holding others fixed.

**5-species scenario:**

Two possible outcomes based on which species coalition survives:

- Basin 1: Species {1, 2, 3} survive

- Basin 2: Species {1, 4, 5} survive

The boundary between these outcomes is fractal but not Wada (only 2
basins).

**8-species scenario:**

Three outcomes based on which invading species (6, 7, or 8) dominates:

- Basin 1: Species 6 dominates

- Basin 2: Species 7 dominates

- Basin 3: Species 8 dominates

With 3+ basins, the system can exhibit Wada basin boundaries.

**Computational complexity warning:**

Competition basin computation is expensive: O(N^2 \\\times\\ T/dt) where
N is resolution and T is simulation time. For a 200x200 grid with
t_max=2000 and dt=0.01, this means 8 billion ODE steps. Recommendations:

- Use resolution 50-100 for exploration

- Increase to 200-400 for publication-quality figures

- The 8-species scenario is slower due to the invasion phase

## Note

**User abort**: Press `Esc` during computation to abort the simulation.
The function checks for user interrupts periodically and will stop
gracefully.

## References

Huisman, J., & Weissing, F. J. (2001). Fundamental unpredictability in
multispecies competition. *American Naturalist*, 157(5), 488-494.
[doi:10.1086/319929](https://doi.org/10.1086/319929)

## See also

[`multispecies_competition`](https://robustecologies.github.io/wadaR/reference/multispecies_competition.md)
for creating competition systems,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for testing Wada property,
[`plot.wada_basins`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: 5-species competition basins                               #
# ===================================================================== #
hw5 <- multispecies_competition(scenario = "5species")
basins5 <- compute_competition_basins(hw5, c(0, 2), c(0, 2),
                                      resolution = 200)
plot(basins5, title = "5-species competition basins",
     colors = c("#E41A1C", "#377EB8"))

# ===================================================================== #
# Example 2: 8-species competition with Wada potential                  #
# ===================================================================== #
hw8 <- multispecies_competition(scenario = "8species")
basins8 <- compute_competition_basins(hw8, c(0, 2), c(0, 2),
                                      resolution = 200)
plot(basins8, title = "8-species competition basins")

# Test for Wada property
result <- wada_merging_method(basins8)
print(result)

# ===================================================================== #
# Example 3: Basin entropy analysis                                     #
# ===================================================================== #
entropy <- basin_entropy(basins8)
print(entropy)
} # }
```
