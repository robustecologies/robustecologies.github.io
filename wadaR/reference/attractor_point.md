# Create a point attractor specification

Defines a fixed-point or periodic attractor by its center location and
convergence radius. A trajectory is considered to have converged when it
enters the ball of the specified radius around the center.

## Usage

``` r
attractor_point(center, radius, label = NULL)
```

## Arguments

- center:

  Numeric vector. Coordinates of the attractor center in phase space.
  Length should match the projection dimension (typically 2 for
  plotting).

- radius:

  Numeric. Convergence radius. Trajectories entering this distance from
  the center are classified as converging to this attractor.

- label:

  Character. Optional human-readable label for the attractor.

## Value

A list with class `"attractor_spec"` containing:

- center:

  Attractor center coordinates

- radius:

  Convergence radius

- label:

  Attractor label (or NULL)

- type:

  "point"

## See also

[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md),
[`attractor_cycle`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md),
[`attractor_exit`](https://robustecologies.github.io/wadaR/reference/attractor_exit.md),
[`attractor_outcome`](https://robustecologies.github.io/wadaR/reference/attractor_outcome.md)

## Examples

``` r
# Fixed point at origin with radius 0.5
attr1 <- attractor_point(center = c(0, 0), radius = 0.5, label = "Origin")
print(attr1)
#> $center
#> [1] 0 0
#> 
#> $radius
#> [1] 0.5
#> 
#> $label
#> [1] "Origin"
#> 
#> $type
#> [1] "point"
#> 
#> attr(,"class")
#> [1] "attractor_spec" "list"          

# Period-1 orbit of forced pendulum at 2*pi
attr2 <- attractor_point(center = c(2*pi, 0), radius = 0.5, label = "+2pi")

# Period-1 orbit at -2*pi
attr3 <- attractor_point(center = c(-2*pi, 0), radius = 0.5, label = "-2pi")

# Use in compiled system
if (FALSE) { # \dontrun{
pendulum <- compiled_system(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -damping * state[1] - sin(state[0]) + forcing * cos(t);
    ',
    attractors = list(attr1, attr2, attr3),
    params = list(damping = 0.2, forcing = 1.66),
    name = "Forced pendulum"
)
} # }
```
