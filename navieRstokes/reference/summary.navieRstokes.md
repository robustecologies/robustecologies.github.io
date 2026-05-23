# Summary method for navieRstokes objects

Provides an extended summary of a Navier-Stokes simulation result with
detailed diagnostics, convergence information, and mass conservation
metrics.

## Usage

``` r
# S3 method for class 'navieRstokes'
summary(object, ...)
```

## Arguments

- object:

  Object of class 'navieRstokes' (output from simulate_navier_stokes)

- ...:

  Additional arguments (not currently used)

## Value

An object of class 'summary.navieRstokes' containing:

- parameters:

  Simulation parameters

- diagnostics:

  Diagnostic metrics

- time_range:

  Time span of simulation

- n_saves:

  Number of saved frames

- velocity_stats:

  Statistics on velocity field

- convergence_info:

  Pressure solver convergence statistics

## See also

[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md),
[`print.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/print.navieRstokes.md),
[`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md),
[`flow.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/flow.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 500, nu = 0.01,
  initial_condition = function(x, y) list(u = 0, v = 0),
  boundary_condition = list(
    type = "dirichlet",
    values = list(u_top = 1, u_bottom = 0, u_left = 0, u_right = 0,
                  v_top = 0, v_bottom = 0, v_left = 0, v_right = 0)
  )
)
summary(result)
} # }
```
