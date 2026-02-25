# Forced damped pendulum system

Creates a forced damped pendulum system, a canonical example for
studying Wada basins of attraction in nonlinear dynamics.

## Usage

``` r
forced_damped_pendulum(damping = 0.2, forcing = 1.66)
```

## Arguments

- damping:

  Damping coefficient \\\gamma\\ (default 0.2).

- forcing:

  Forcing amplitude \\F\\. Key values:

  - \\F = 1.66\\: Full Wada basins (classic example)

  - \\F = 1.71\\: Partially Wada basins

  - \\F \> 2.0\\: Chaotic dynamics with complex basin structure

## Value

A list (system object) containing:

- system:

  Function(state, t) defining the ODE right-hand side

- attractors:

  List of 3 attractors with center and radius

- params:

  List with damping and forcing values

- description:

  Human-readable description

## Details

The forced damped pendulum is governed by the second-order ODE:
\$\$\ddot{x} + \gamma \dot{x} + \sin(x) = F \cos(t)\$\$

This system is written as a first-order system in phase space \\(x,
v)\\: \$\$\dot{x} = v\$\$ \$\$\dot{v} = -\gamma v - \sin(x) + F
\cos(t)\$\$

The pendulum exhibits multistability with three coexisting period-1
attractors located approximately at:

- Attractor 1: \\(x, v) \approx (0, 0)\\

- Attractor 2: \\(x, v) \approx (2\pi, 0)\\

- Attractor 3: \\(x, v) \approx (-2\pi, 0)\\

For \\\gamma = 0.2\\ and \\F \approx 1.66\\, the basins exhibit the Wada
property: every boundary point of one basin is also a boundary point of
the other two basins. This remarkable topological property implies that
arbitrarily small uncertainties in initial conditions can lead to any of
the three attractors.

## References

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and basin
cells. *Physica D: Nonlinear Phenomena*, 90(3), 242-261.
[doi:10.1016/0167-2789(95)00249-9](https://doi.org/10.1016/0167-2789%2895%2900249-9)

Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., & Sanjuan, M.
A. F. (2016). Basin entropy: A new tool to analyze uncertainty in
dynamical systems. *Scientific Reports*, 6, 31416.
[doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins of attraction,
[`henon_heiles_system`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md)
for another Wada-exhibiting system,
[`wada_grid_method`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
for testing Wada property.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Classic Wada basins (F = 1.66)                             #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66, damping = 0.2)
print(pendulum$description)

basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 500)
plot(basins, title = "Wada basins (F = 1.66)")

# ===================================================================== #
# Example 2: Partially Wada basins (F = 1.71)                           #
# ===================================================================== #
pendulum_partial <- forced_damped_pendulum(forcing = 1.71)
basins_partial <- compute_basins(pendulum_partial, c(-pi, pi), c(-3, 3),
                                 resolution = 500)
plot(basins_partial, title = "Partially Wada (F = 1.71)")

# ===================================================================== #
# Example 3: Highly chaotic regime (F = 2.3)                            #
# ===================================================================== #
pendulum_chaotic <- forced_damped_pendulum(forcing = 2.3)
basins_chaotic <- compute_basins(pendulum_chaotic, c(-pi, pi), c(-3, 3),
                                 resolution = 500, t_max = 150)
plot(basins_chaotic, title = "Chaotic regime (F = 2.3)")

# ===================================================================== #
# Example 4: Test for Wada property                                     #
# ===================================================================== #
result <- wada_grid_method(basins)
print(result)  # Shows if basins are Wada
plot(result, basins = basins)
} # }
```
