# Detect attractors automatically from trajectory analysis

Analyzes trajectories from multiple initial conditions to automatically
identify attractors (fixed points, limit cycles, or strange attractors)
and estimate appropriate convergence radii. Uses C++/OpenMP
parallelization for high performance with interrupt support.

## Usage

``` r
detect_attractors(
  cpp_dynamics,
  params = list(),
  dim = 2L,
  type = c("ode", "map"),
  x_range = c(-2, 2),
  y_range = c(-2, 2),
  n_trajectories = 200L,
  t_max = 100,
  t_transient = 50,
  dt = 0.01,
  min_attractors = 3L,
  max_attractors = 10L,
  cluster_method = c("hclust", "kmeans", "dbscan"),
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- cpp_dynamics:

  Character string containing C++ code for the dynamics. Same format as
  for
  [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md):
  `deriv[0] = ...; deriv[1] = ...;`

- params:

  Named list of scalar numeric parameters.

- dim:

  Integer. Dimension of the state space (default 2).

- type:

  Character. System type: "ode" (default) or "map".

- x_range:

  Numeric vector of length 2. Range for initial x coordinates.

- y_range:

  Numeric vector of length 2. Range for initial y coordinates.

- n_trajectories:

  Integer. Number of trajectories to sample.

- t_max:

  Numeric. Integration time for each trajectory.

- t_transient:

  Numeric. Transient time to discard before recording.

- dt:

  Numeric. Integration time step (for ODEs).

- min_attractors:

  Integer. Minimum number of attractors to detect.

- max_attractors:

  Integer. Maximum number of attractors to consider.

- cluster_method:

  Character. Clustering method: "hclust" (hierarchical), "kmeans", or
  "dbscan" (requires dbscan package).

- n_cores:

  Integer. Number of CPU cores to use for parallel computation. Default
  is `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`.

- verbose:

  Logical. Show progress information.

## Value

A list of class `"attractor_detection"` containing:

- attractors:

  List of detected attractor specifications

- n_attractors:

  Number of attractors found

- attractor_types:

  Character vector of attractor types

- cluster_info:

  Clustering details for diagnostics

- endpoints:

  Matrix of trajectory endpoints

## Details

The detection algorithm:

1.  Compiles C++ code for parallel trajectory integration

2.  Samples random initial conditions from specified ranges

3.  Integrates trajectories using OpenMP parallelization

4.  Discards transient behavior

5.  Clusters trajectory endpoints to identify attractors

6.  Classifies attractors and estimates convergence radii

**Clustering methods:**

- `"hclust"`: Hierarchical clustering (default, most robust)

- `"kmeans"`: K-means clustering

- `"dbscan"`: Density-based clustering (requires dbscan package)

## See also

[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md),
[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ============================================================
# Example 1: Bistable system (two fixed points)
# ============================================================
# x'' + 0.3*x' + x*(x^2 - 1) = 0
# Fixed points at x = -1, 0, +1 (origin is unstable)

detected <- detect_attractors(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = -0.3 * v - x * (x*x - 1);
    ',
    params = list(),
    dim = 2,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    n_trajectories = 200,
    t_max = 50,
    t_transient = 20,
    min_attractors = 2,
    max_attractors = 5,
    cluster_method = "hclust",
    verbose = TRUE
)

# Print detection results
print(detected)

# Visualize trajectory endpoints and detected attractors
plot(detected)

# ============================================================
# Example 2: Triple-well system (three fixed points)
# ============================================================
# x'' + 0.5*x' + x*(x^2 - 4) = 0
# Fixed points at x = -2, 0, +2

detected_triple <- detect_attractors(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -0.5 * state[1] - state[0] * (state[0]*state[0] - 4);
    ',
    params = list(),
    x_range = c(-3, 3),
    y_range = c(-3, 3),
    n_trajectories = 300,
    t_max = 100,
    t_transient = 50,
    min_attractors = 3,
    verbose = TRUE
)

print(detected_triple)

# Use detected attractors in a compiled system
triple_well <- compiled_system(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -0.5 * state[1] - state[0] * (state[0]*state[0] - 4);
    ',
    attractors = detected_triple$attractors,
    params = list(),
    name = "Triple-well system"
)

basins <- compute_basins(triple_well, c(-3, 3), c(-3, 3), resolution = 200)
plot(basins)

# ============================================================
# Example 3: Using different clustering methods
# ============================================================

# K-means clustering
detected_km <- detect_attractors(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -0.3 * state[1] - state[0] * (state[0]*state[0] - 1);
    ',
    params = list(),
    n_trajectories = 200,
    cluster_method = "kmeans",
    min_attractors = 2,
    verbose = TRUE
)

# DBSCAN (density-based, requires dbscan package)
detected_db <- detect_attractors(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -0.3 * state[1] - state[0] * (state[0]*state[0] - 1);
    ',
    params = list(),
    n_trajectories = 200,
    cluster_method = "dbscan",
    min_attractors = 2,
    verbose = TRUE
)
} # }
```
