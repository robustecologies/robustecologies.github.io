# Huisman-Weissing multispecies competition system

Creates a Huisman-Weissing resource competition system for basin
analysis. This ecological model demonstrates that competition for
resources can produce fractal basin boundaries and Wada-like
unpredictability in species outcomes.

## Usage

``` r
multispecies_competition(
  scenario = c("5species", "8species"),
  dilution = 0.25,
  supply = c(10, 10, 10),
  t_invade = 1000,
  N_invade = 0.1
)
```

## Arguments

- scenario:

  Character. Competition scenario:

  - `"5species"`: 5 species competing for 3 resources. Produces 2
    possible outcomes (binary basin structure).

  - `"8species"`: 8 species with invasion. Produces 3+ possible
    outcomes, suitable for Wada detection.

- dilution:

  Numeric. Resource dilution rate \\D\\ (default 0.25 d\\^{-1}\\).

- supply:

  Numeric vector. Resource supply concentrations \\S_j\\ (default
  `c(10, 10, 10)` \\\mu\\mol L\\^{-1}\\).

- t_invade:

  Numeric. For 8-species scenario, time when species 6-8 invade (default
  1000 days).

- N_invade:

  Numeric. Initial abundance of invading species (default 0.1).

## Value

A list (system object) containing:

- system:

  "huisman_weissing" string identifier

- params:

  List of all model parameters (K, C, r, m, D, S, etc.)

- attractors:

  List of outcome-based attractor specifications

- projection:

  Named vector indicating which species to vary for 2D basin computation
  (default: N2 and N4)

- scenario:

  The scenario name

- type:

  "competition"

- description:

  Human-readable description

## Details

The Huisman-Weissing model describes \\n\\ species competing for \\k\\
resources in a chemostat. The dynamics are governed by:

**Species dynamics:** \$\$\frac{dN_i}{dt} = N_i \left( \mu_i(R_1,
\ldots, R_k) - m_i \right)\$\$

**Resource dynamics:** \$\$\frac{dR_j}{dt} = D(S_j - R_j) -
\sum\_{i=1}^{n} c\_{ji} \mu_i N_i\$\$

**Growth rate (Liebig's law with Monod kinetics):** \$\$\mu_i =
\min\_{j} \left( r_i \frac{R_j}{K\_{ji} + R_j} \right)\$\$

where:

- \\N_i\\: abundance of species \\i\\

- \\R_j\\: concentration of resource \\j\\

- \\\mu_i\\: specific growth rate of species \\i\\

- \\m_i\\: mortality rate of species \\i\\

- \\D\\: dilution rate

- \\S_j\\: supply concentration of resource \\j\\

- \\K\_{ji}\\: half-saturation constant of species \\i\\ for resource
  \\j\\

- \\c\_{ji}\\: resource content (consumption coefficient)

- \\r_i\\: maximum growth rate

The key finding of Huisman & Weissing is that when \\n \> k\\ (more
species than resources), the system can exhibit chaotic dynamics and
fractal basin boundaries. Small changes in initial conditions can lead
to completely different community compositions.

**5-species scenario:**

Five species compete for 3 resources, producing two possible outcomes:

- Outcome 1: Species {1, 2, 3} survive

- Outcome 2: Species {1, 4, 5} survive

The basin boundary between these outcomes is fractal.

**8-species scenario:**

Eight species compete, with species 6, 7, 8 invading at \\t = 1000\\
days. This produces three distinct outcomes based on which invader
dominates:

- Outcome 1: Species 6 dominates

- Outcome 2: Species 7 dominates

- Outcome 3: Species 8 dominates

With 3+ outcomes, the system can exhibit Wada basin boundaries.

## References

Huisman, J., & Weissing, F. J. (2001). Fundamental unpredictability in
multispecies competition. *American Naturalist*, 157(5), 488-494.
[doi:10.1086/319929](https://doi.org/10.1086/319929)

Huisman, J., & Weissing, F. J. (1999). Biodiversity of plankton by
species oscillations and chaos. *Nature*, 402(6760), 407-410.
[doi:10.1038/46540](https://doi.org/10.1038/46540)

## See also

[`compute_competition_basins`](https://robustecologies.github.io/wadaR/reference/compute_competition_basins.md)
for computing basins,
[`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
for comparison with physical systems,
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for Wada detection.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: 5-species competition (fractal boundary, not Wada)         #
# ===================================================================== #
hw5 <- multispecies_competition(scenario = "5species")
print(hw5$description)

# Compute basins varying initial N2 and N4
basins5 <- compute_competition_basins(hw5, c(0, 2), c(0, 2),
                                      resolution = 200)
plot(basins5, title = "5-species competition basins")

# ===================================================================== #
# Example 2: 8-species competition (potential Wada basins)              #
# ===================================================================== #
hw8 <- multispecies_competition(scenario = "8species")
basins8 <- compute_competition_basins(hw8, c(0, 2), c(0, 2),
                                      resolution = 200)
plot(basins8, title = "8-species competition basins")

# Test for Wada property
result <- wada_merging_method(basins8)
print(result)

# ===================================================================== #
# Example 3: Custom parameters                                          #
# ===================================================================== #
hw_custom <- multispecies_competition(scenario = "5species",
                                      dilution = 0.3,
                                      supply = c(12, 8, 10))
} # }
```
