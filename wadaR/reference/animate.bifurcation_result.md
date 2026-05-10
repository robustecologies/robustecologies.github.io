# Animate bifurcation results

Creates an animated GIF showing basin evolution across parameter values.
Uses parallel frame rendering for high performance.

## Usage

``` r
# S3 method for class 'bifurcation_result'
animate(
  x,
  filename = NULL,
  fps = 5,
  width = 600,
  height = 600,
  palette = "turbo",
  n_cores = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A bifurcation_result object.

- filename:

  Character. Output filename. If `NULL` (default), saves to a temporary
  file and returns the path.

- fps:

  Integer. Frames per second (default 5).

- width, height:

  Integer. Animation dimensions in pixels.

- palette:

  Character. Color palette name (default "turbo").

- n_cores:

  Integer. Number of CPU cores for parallel frame rendering. Default is
  `NULL`, which uses `parallel::detectCores(logical = FALSE) - 1`.

- verbose:

  Logical. Show progress information (default `TRUE`).

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns the output file path.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create animation
animate(bif_result, filename = "basin_evolution.gif", fps = 10)

# Parallel rendering (auto-detected cores)
animate(bif_result, filename = "basins.gif")
} # }
```
