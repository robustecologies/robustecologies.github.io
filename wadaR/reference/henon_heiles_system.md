# Henon-Heiles Hamiltonian system

Creates a Henon-Heiles system for escape basin analysis. This
conservative Hamiltonian system exhibits chaotic scattering and Wada
basins when particles can escape through three symmetric channels.

## Usage

``` r
henon_heiles_system(energy = 0.2)
```

## Arguments

- energy:

  Total energy \\E\\. Must be above the critical escape energy \\E_c =
  1/6 \approx 0.1667\\ for trajectories to escape.

## Value

A list (system object) containing:

- system:

  Function(state, t) defining the 4D ODE right-hand side

- exits:

  List of 3 exit specifications with angle and label

- params:

  List with energy and critical_energy values

- description:

  Human-readable description

## Details

The Henon-Heiles system is a paradigmatic model in nonlinear dynamics,
originally introduced to study stellar motion in galaxies. The
Hamiltonian is: \$\$H = \frac{1}{2}(p_x^2 + p_y^2) + \frac{1}{2}(x^2 +
y^2) + x^2 y - \frac{y^3}{3}\$\$

The corresponding equations of motion are: \$\$\dot{x} = p_x, \quad
\dot{y} = p_y\$\$ \$\$\dot{p_x} = -x - 2xy, \quad \dot{p_y} = -y - x^2 +
y^2\$\$

The potential has a triangular symmetry with three saddle points at the
critical energy \\E_c = 1/6\\. Above this energy, particles can escape
through three channels (exits) located at 120-degree intervals:

- Exit 1: Top (\\\theta = \pi/2\\)

- Exit 2: Bottom-left (\\\theta = 7\pi/6\\)

- Exit 3: Bottom-right (\\\theta = 11\pi/6\\)

The escape basins exhibit the Wada property: every point on the boundary
of one exit basin is simultaneously on the boundary of all three basins.
This makes final-state prediction fundamentally uncertain for
trajectories starting near the boundary.

## References

Henon, M., & Heiles, C. (1964). The applicability of the third integral
of motion: Some numerical experiments. *The Astronomical Journal*, 69,
73. [doi:10.1086/109234](https://doi.org/10.1086/109234)

Aguirre, J., Vallejo, J. C., & Sanjuan, M. A. F. (2001). Wada basins and
chaotic invariant sets in the Henon-Heiles system. *Physical Review E*,
64(6), 066208.
[doi:10.1103/PhysRevE.64.066208](https://doi.org/10.1103/PhysRevE.64.066208)

Seoane, J. M., & Sanjuan, M. A. F. (2013). New developments in classical
chaotic scattering. *Reports on Progress in Physics*, 76(1), 016001.
[doi:10.1088/0034-4885/76/1/016001](https://doi.org/10.1088/0034-4885/76/1/016001)

## See also

[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing escape basins,
[`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
for another Wada-exhibiting system,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for testing Wada property.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Escape basins at E = 0.2                                   #
# ===================================================================== #
hh <- henon_heiles_system(energy = 0.2)
print(hh$description)

basins <- compute_basins(hh, c(-0.4, 0.4), c(-0.4, 0.4),
                         resolution = 500, t_max = 200)
plot(basins, title = "Henon-Heiles escape basins (E = 0.2)",
     colors = c("#E41A1C", "#377EB8", "#4DAF4A"))

# ===================================================================== #
# Example 2: Higher energy (more chaotic)                               #
# ===================================================================== #
hh_high <- henon_heiles_system(energy = 0.25)
basins_high <- compute_basins(hh_high, c(-0.5, 0.5), c(-0.5, 0.5),
                              resolution = 500, t_max = 300)
plot(basins_high, title = "Henon-Heiles (E = 0.25)")

# ===================================================================== #
# Example 3: Near critical energy (barely escaping)                     #
# ===================================================================== #
hh_critical <- henon_heiles_system(energy = 0.17)  # Just above Ec
basins_crit <- compute_basins(hh_critical, c(-0.3, 0.3), c(-0.3, 0.3),
                              resolution = 500, t_max = 500)
plot(basins_crit, title = "Near critical energy (E = 0.17)")

# ===================================================================== #
# Example 4: Test for Wada property                                     #
# ===================================================================== #
result <- wada_merging_method(basins)
print(result)
plot(result)
} # }
```
