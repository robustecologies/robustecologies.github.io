# Create an escape channel specification

Defines an escape channel for systems where trajectories can escape to
infinity in specific directions (like the Henon-Heiles system).

## Usage

``` r
attractor_exit(angle, label = NULL)
```

## Arguments

- angle:

  Numeric. Angle in radians defining the escape direction. Measured
  counterclockwise from the positive x-axis.

- label:

  Character. Optional human-readable label for the exit.

## Value

A list with class `"attractor_spec"` containing:

- angle:

  Escape direction angle

- label:

  Exit label (or NULL)

- type:

  "exit"

## See also

[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)

## Examples

``` r
# Three exits at 120-degree intervals (like Henon-Heiles)
exit1 <- attractor_exit(angle = pi/2, label = "Top")
exit2 <- attractor_exit(angle = pi/2 + 2*pi/3, label = "Bottom-left")
exit3 <- attractor_exit(angle = pi/2 + 4*pi/3, label = "Bottom-right")

print(exit1)
#> $angle
#> [1] 1.570796
#> 
#> $label
#> [1] "Top"
#> 
#> $type
#> [1] "exit"
#> 
#> attr(,"class")
#> [1] "attractor_spec" "list"          

# The Henon-Heiles system has three symmetric escape channels
# where trajectories can escape to infinity above the critical energy
all_exits <- list(exit1, exit2, exit3)
sapply(all_exits, function(e) e$angle * 180 / pi)  # Angles in degrees
#> [1]  90 210 330
```
