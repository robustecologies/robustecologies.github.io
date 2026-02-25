# Fractal basins in multispecies competition

``` r
library(wadaR)
library(ggplot2)
```

## Introduction

This vignette demonstrates the application of Wada basin detection
methods to ecological systems, specifically the Huisman-Weissing
multispecies competition model. This model reveals that competition for
resources can produce fractal basin boundaries and fundamental
unpredictability in species outcomes.

### Ecological motivation

Classical competition theory (competitive exclusion principle) predicts
that \\n\\ species competing for \\k\\ resources should result in at
most \\k\\ coexisting species. However, Huisman & Weissing (1999, 2001)
demonstrated that with oscillatory and chaotic dynamics:

- More than \\k\\ species can coexist on \\k\\ resources
- The final community composition may be fundamentally unpredictable
- Basin boundaries between different outcomes can be fractal

These findings have profound implications for:

- **Harmful algal blooms**: Cannot reliably predict which phytoplankton
  species will dominate
- **Ecosystem management**: Small perturbations may have large,
  unpredictable effects
- **Biodiversity maintenance**: Chaotic dynamics may explain high
  diversity in plankton communities

------------------------------------------------------------------------

## Mathematical model

### Model formulation

The Huisman-Weissing model describes \\n\\ species competing for \\k\\
abiotic resources in a chemostat-type environment.

#### State variables

- \\N_i(t)\\: abundance of species \\i\\ at time \\t\\
- \\R_j(t)\\: concentration of resource \\j\\ at time \\t\\

#### Species dynamics

\\\frac{dN_i}{dt} = N_i \left\[ \mu_i(R_1, \ldots, R_k) - m_i \right\]\\

where \\\mu_i\\ is the specific growth rate and \\m_i\\ is the mortality
rate.

#### Resource dynamics

\\\frac{dR_j}{dt} = D(S_j - R_j) - \sum\_{i=1}^{n} c\_{ji} \mu_i N_i\\

where:

- \\D\\: dilution (turnover) rate
- \\S_j\\: supply concentration of resource \\j\\
- \\c\_{ji}\\: resource content (consumption coefficient)

#### Growth kinetics

The specific growth rate follows **Monod kinetics** combined with
**Liebig’s law of the minimum**:

\\\mu_i = \min\_{j=1,\ldots,k} \left( \frac{r_i R_j}{K\_{ji} + R_j}
\right)\\

where:

- \\r_i\\: maximum specific growth rate
- \\K\_{ji}\\: half-saturation constant for resource \\j\\ of species
  \\i\\

The \\\min(\cdot)\\ operator implements Liebig’s law: growth is limited
by the most limiting resource.

### Parameter values

Following Huisman & Weissing (2001), we use biologically realistic
values for phytoplankton:

| Parameter | Value | Unit                  | Description         |
|-----------|-------|-----------------------|---------------------|
| \\r_i\\   | 1.0   | d\\^{-1}\\            | Maximum growth rate |
| \\m_i\\   | 0.25  | d\\^{-1}\\            | Mortality rate      |
| \\D\\     | 0.25  | d\\^{-1}\\            | Dilution rate       |
| \\S_j\\   | 10    | \\\mu\\mol L\\^{-1}\\ | Resource supply     |

------------------------------------------------------------------------

## Competition scenarios

The `wadaR` package implements two scenarios from the original paper.

### 5-species scenario

Five species compete for 3 resources, creating intransitive competitive
relationships (rock-paper-scissors dynamics):

``` r
hw5 <- multispecies_competition(scenario = "5species")
print(hw5)
```

This scenario produces **two possible outcomes**:

- **Outcome 1**: Species {1, 2, 3} survive (species 4, 5 go extinct)
- **Outcome 2**: Species {1, 4, 5} survive (species 2, 3 go extinct)

The boundary between these outcomes is fractal, but since there are only
2 basins, the Wada property does not apply (Wada requires \\\geq 3\\
basins).

### 8-species scenario

Eight species compete, with species 6, 7, and 8 invading at \\t = 1000\\
days:

``` r
hw8 <- multispecies_competition(scenario = "8species")
print(hw8)
```

This scenario produces **three possible outcomes** based on which
invading species dominates:

- **Outcome 1**: Species 6 dominates
- **Outcome 2**: Species 7 dominates
- **Outcome 3**: Species 8 dominates

With 3+ basins, this system can potentially exhibit Wada basin
boundaries.

------------------------------------------------------------------------

## Time series dynamics

Before computing basins, it’s instructive to examine the temporal
dynamics.

### Chaotic dynamics in 5-species competition

``` r
hw5 <- multispecies_competition(scenario = "5species")

# Simulate dynamics
result <- simulate_competition(hw5, t_max = 4000)

# Plot time series
matplot(result$time, result$N, type = "l", lty = 1, lwd = 1.5,
        xlab = "Time (days)", ylab = "Species abundance",
        main = "5-species chaotic competition",
        col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
legend("topright", paste0("Sp", 1:5), col = c("#E41A1C", "#377EB8", "#4DAF4A",
       "#984EA3", "#FF7F00"), lty = 1, lwd = 1.5, cex = 0.8)
```

The system exhibits **transient chaos**: the dynamics are chaotic for a
finite time before eventually settling on one of the two limit cycles.
The duration of transient chaos depends sensitively on initial
conditions.

### 8-species invasion dynamics

``` r
hw8 <- multispecies_competition(scenario = "8species")

# Custom initial conditions
N0 <- c(0.1, 0.5, 0.1, 0.8, 0.1, 0, 0, 0)  # Species 6-8 start at 0

# Simulate with invasion
result8 <- simulate_competition(hw8, N0 = N0, t_max = 2000)

# Plot
matplot(result8$time, result8$N, type = "l", lty = 1, lwd = 1.5,
        xlab = "Time (days)", ylab = "Species abundance",
        main = "8-species competition with invasion")
abline(v = 1000, lty = 2, col = "gray50")
text(1050, max(result8$N) * 0.9, "Invasion", adj = 0)
```

------------------------------------------------------------------------

## Basin computation

### Computational complexity

Basin computation for competition models is computationally expensive.
The complexity is:

\\O(N^2 \times T/\Delta t \times n\_{vars})\\

where:

- \\N\\: grid resolution
- \\T\\: simulation time (2000-2500 days)
- \\\Delta t\\: time step (0.01 days)
- \\n\_{vars}\\: number of state variables (species + resources)

For a 200×200 grid with the 8-species scenario, this means
approximately:

- 40,000 initial conditions
- 200,000 integration steps per condition
- 11 state variables (8 species + 3 resources)
- \\\approx 8.8 \times 10^{10}\\ RK4 evaluations

**Recommendation**: Start with resolution 50-100 for exploration,
increase to 200-400 for final figures.

### 5-species basins

``` r
hw5 <- multispecies_competition(scenario = "5species")

# Compute basins (low resolution for exploration)
basins5 <- compute_competition_basins(
  hw5,
  x_range = c(0, 2),
  y_range = c(0, 2),
  resolution = 200,
  t_max = 2500,
  verbose = TRUE
)

# Plot
plot(basins5,
     title = "5-species competition basins",
     colors = c("#E41A1C", "#377EB8"),
     xlab = expression(N[2](0)),
     ylab = expression(N[4](0)))
```

The fractal boundary is clearly visible. Small changes in initial
conditions near the boundary lead to completely different community
compositions.

### 8-species basins

``` r
hw8 <- multispecies_competition(scenario = "8species")

# Compute basins
basins8 <- compute_competition_basins(
  hw8,
  x_range = c(0, 2),
  y_range = c(0, 2),
  resolution = 200,
  t_max = 2000,
  verbose = TRUE
)

# Plot
plot(basins8,
     title = "8-species competition basins",
     colors = c("#E41A1C", "#377EB8", "#4DAF4A"),
     xlab = expression(N[2](0)),
     ylab = expression(N[4](0)))
```

------------------------------------------------------------------------

## Wada analysis

### Testing for Wada property

With 3 basins in the 8-species scenario, we can test for the Wada
property using the merging method:

``` r
# Test for Wada using merging method
result <- wada_merging_method(basins8, verbose = TRUE)
print(result)
```

``` r
plot(result, basins = basins8)
```

### Basin entropy

The basin entropy quantifies the unpredictability of the system:

``` r
entropy <- basin_entropy(basins8, box_size = 10)
print(entropy)
```

``` r
plot(entropy, basins = basins8)
```

Higher entropy values indicate regions where small measurement
uncertainties lead to large uncertainties in the predicted outcome.

------------------------------------------------------------------------

## Ecological interpretation

### Competitive exclusion vs. chaos

The classical competitive exclusion principle states that \\n\\ species
competing for \\k\\ resources can support at most \\k\\ coexisting
species at equilibrium. However, the Huisman-Weissing model demonstrates
that:

1.  **Oscillations extend coexistence**: Non-equilibrium dynamics (limit
    cycles) can maintain more than \\k\\ species
2.  **Chaos enables supersaturation**: Chaotic dynamics can maintain 5+
    species on 3 resources
3.  **Outcomes become unpredictable**: Basin boundaries can be fractal,
    making prediction fundamentally impossible

### Implications for ecosystem management

The fractal basin structure has important practical implications:

- **Prediction limits**: Even with perfect knowledge of model
  parameters, small uncertainties in initial conditions can make
  long-term predictions impossible
- **Tipping points**: The system may appear stable but be arbitrarily
  close to a basin boundary
- **Intervention design**: Management actions must account for the
  fractal structure of basins

### The paradox of the plankton

Hutchinson (1961) posed the “paradox of the plankton”: why do so many
phytoplankton species coexist in seemingly homogeneous environments with
limited resources? The Huisman-Weissing model provides one resolution:

> “Competitive chaos enables the indefinite coexistence of many species
> on few resources, thus providing a solution to the paradox of the
> plankton.” - Huisman & Weissing (1999)

------------------------------------------------------------------------

## Computational details

### Numerical methods

The system is integrated using the **4th-order Runge-Kutta method** with
fixed time step \\\Delta t = 0.01\\ days:

\\\mathbf{k}\_1 = \mathbf{f}(\mathbf{y}\_n)\\ \\\mathbf{k}\_2 =
\mathbf{f}(\mathbf{y}\_n + \frac{\Delta t}{2}\mathbf{k}\_1)\\
\\\mathbf{k}\_3 = \mathbf{f}(\mathbf{y}\_n + \frac{\Delta
t}{2}\mathbf{k}\_2)\\ \\\mathbf{k}\_4 = \mathbf{f}(\mathbf{y}\_n +
\Delta t \cdot \mathbf{k}\_3)\\ \\\mathbf{y}\_{n+1} = \mathbf{y}\_n +
\frac{\Delta t}{6}(\mathbf{k}\_1 + 2\mathbf{k}\_2 + 2\mathbf{k}\_3 +
\mathbf{k}\_4)\\

The implementation uses OpenMP parallelization for efficient computation
across multiple CPU cores.

### Winner determination

For the **5-species scenario**, a species group “wins” if:

- All species in that group have abundance \> \\\theta\_{extinct}\\
  (default 0.01)
- Species not in that group have abundance \< \\\theta\_{extinct}\\

For the **8-species scenario**, the winner is determined by which
invading species (6, 7, or 8) has the highest abundance above a
dominance threshold (default 1.0).

------------------------------------------------------------------------

## References

1.  Huisman, J. & Weissing, F.J. (2001). Fundamental unpredictability in
    multispecies competition. *American Naturalist* 157: 488-494. DOI:
    10.1086/319929

2.  Huisman, J. & Weissing, F.J. (1999). Biodiversity of plankton by
    species oscillations and chaos. *Nature* 402: 407-410. DOI:
    10.1038/46540

3.  Hutchinson, G.E. (1961). The paradox of the plankton. *American
    Naturalist* 95: 137-145.

4.  Tilman, D. (1982). *Resource competition and community structure*.
    Princeton University Press.

5.  Monod, J. (1950). La technique de culture continue, theorie et
    applications. *Annales de l’Institut Pasteur* 79: 390-410.

6.  Grebogi, C., Ott, E. & Yorke, J.A. (1987). Chaos, strange
    attractors, and fractal basin boundaries in nonlinear dynamics.
    *Science* 238: 632-638.
