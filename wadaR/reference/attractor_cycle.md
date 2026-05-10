# Create a limit cycle attractor specification

Defines a limit cycle attractor. For basin computation, convergence is
detected when trajectories enter a neighborhood of the cycle's
approximate center.

## Usage

``` r
attractor_cycle(center, radius, period = NULL, label = NULL)
```

## Arguments

- center:

  Numeric vector. Approximate center of the limit cycle in phase space.

- radius:

  Numeric. Convergence radius for detecting trajectories that have
  locked onto the cycle.

- period:

  Numeric. Optional period of the limit cycle. Not used for basin
  computation but useful for documentation.

- label:

  Character. Optional human-readable label.

## Value

A list with class `"attractor_spec"` containing:

- center:

  Cycle center coordinates

- radius:

  Convergence radius

- period:

  Cycle period (or NULL)

- label:

  Attractor label (or NULL)

- type:

  "cycle"

## See also

[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)

## Examples

``` r
# Limit cycle centered at (1, 0) with period 2*pi
cycle <- attractor_cycle(
    center = c(1, 0),
    radius = 0.2,
    period = 2*pi,
    label = "Limit cycle"
)
print(cycle)
#> $center
#> [1] 1 0
#> 
#> $radius
#> [1] 0.2
#> 
#> $period
#> [1] 6.283185
#> 
#> $label
#> [1] "Limit cycle"
#> 
#> $type
#> [1] "cycle"
#> 
#> attr(,"class")
#> [1] "attractor_spec" "list"          

# Multiple limit cycles for Van der Pol oscillator
large_orbit <- attractor_cycle(c(2, 0), 0.5, label = "Large amplitude")
small_orbit <- attractor_cycle(c(-2, 0), 0.5, label = "Small amplitude")

# Use in compiled system
if (FALSE) { # \dontrun{
vdp <- compiled_system(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = mu * (1.0 - x*x) * v - x;
    ',
    attractors = list(
        large_orbit,
        small_orbit,
        attractor_point(c(0, 0), 0.3, "Unstable origin")
    ),
    params = list(mu = 1.0),
    name = "Van der Pol oscillator"
)
} # }
```
