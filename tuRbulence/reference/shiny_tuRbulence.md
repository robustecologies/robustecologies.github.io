# Interactive Shiny application for tuRbulence

Launches an interactive Shiny dashboard for exploring the dynamical
systems implemented in the tuRbulence package. The application provides
a unified interface for simulating all seven systems, visualizing
attractors, and computing dynamical invariants such as Lyapunov
exponents and embedding dimensions.

## Usage

``` r
shiny_tuRbulence(launch.browser = TRUE)
```

## Arguments

- launch.browser:

  Logical; if TRUE (default), opens the app in the system's default web
  browser. If FALSE, returns the app object.

## Value

If `launch.browser = FALSE`, returns a Shiny app object. Otherwise,
launches the app and returns NULL invisibly.

## Details

The application includes the following sections:

**Home:** Overview of package features with quick start guide.

**Classical systems:** Interactive simulation of Lorenz, Roessler, and
Lorenz-84 systems with 3D attractor visualization.

**Geophysical models:** Von Karman turbulent flow and Charney-DeVore
atmospheric blocking model with specialized analysis tools.

**Climate models:** Stommel thermohaline circulation and Hasselmann
stochastic climate model.

**Analysis tools:** Cao's method for embedding dimension and Lyapunov
exponent estimation using Wolf and Rosenstein algorithms.

**Theory:** Mathematical background for each system with equations.

## See also

[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md),
[`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md),
[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md),
[`stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Launch the application
shiny_tuRbulence()

# Get app object without launching
app <- shiny_tuRbulence(launch.browser = FALSE)
} # }
```
