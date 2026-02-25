# Wada basins: theory and detection methods

``` r
library(wadaR)
library(ggplot2)
```

## Introduction to Wada basins

### Historical background

The concept of Wada basins originates from topology. In 1917, Kunizo
Yoneyama described a construction that divides a region of the plane
into three or more connected sets sharing a common boundary, attributing
the original procedure to his advisor Takeo Wada. These structures,
known as “Wada lakes,” were initially studied within a purely
topological context.

The connection to dynamical systems was established by Kennedy and Yorke
(1991), who demonstrated that basins of attraction in chaotic systems
can exhibit the Wada property. The forced damped pendulum became the
paradigmatic example of a dynamical system with Wada basins.

### Mathematical definitions

#### Basin of attraction

For a dynamical system with attractors \\A_1, \ldots, A\_{N_A}\\, the
**basin of attraction** \\B_i\\ of attractor \\A_i\\ is the set of all
initial conditions that asymptotically converge to \\A_i\\:

\\B_i = \\x \in \Omega : \lim\_{t \to \infty} \phi_t(x) \in A_i\\\\

where \\\phi_t\\ is the flow of the dynamical system.

#### Basin boundary

The **boundary** \\\partial B_i\\ of basin \\B_i\\ is the set of points
\\p\\ such that every neighborhood of \\p\\ contains points from both
\\B_i\\ and its complement:

\\\partial B_i = \\p : \forall \varepsilon \> 0, \\ B(p, \varepsilon)
\cap B_i \neq \emptyset \text{ and } B(p, \varepsilon) \cap B_i^c \neq
\emptyset\\\\

#### Wada point

A point \\p\\ is a **Wada point** if it is on the boundary of all basins
simultaneously:

\\p \text{ is Wada} \iff \forall i \in \\1, \ldots, N_A\\, \\ p \in
\partial B_i\\

#### Wada basin

A collection of basins \\\\B_1, \ldots, B\_{N_A}\\\\ with \\N_A \geq 3\\
is **Wada** if every boundary point is a Wada point:

\\\partial B_1 = \partial B_2 = \cdots = \partial B\_{N_A}\\

**Theorem (Kuratowski):** If a boundary separates three or more
connected regions in the plane, then the boundary must be an
indecomposable continuum.

### Implications for predictability

Wada basins represent an extreme form of unpredictability: a small
perturbation of an initial condition on the Wada boundary can lead to
convergence to *any* of the system’s attractors. This has profound
implications for:

- **Weather and climate prediction**: Small uncertainties in atmospheric
  initial conditions can lead to drastically different outcomes
- **Ecological system dynamics**: Population models with multiple stable
  states may exhibit Wada boundaries
- **Engineering control systems**: Understanding basin structure is
  crucial for robust control design
- **Celestial mechanics**: The three-body problem and spacecraft
  trajectories can exhibit Wada basins

------------------------------------------------------------------------

## Example dynamical systems

The `wadaR` package provides several example systems known to exhibit
Wada basins.

### The forced damped pendulum

The forced damped pendulum is the canonical example for studying Wada
basins. It is governed by the second-order ODE:

\\\ddot{x} + \gamma \dot{x} + \sin(x) = F \cos(t)\\

where \\\gamma\\ is the damping coefficient and \\F\\ is the forcing
amplitude. This system exhibits three coexisting period-1 attractors
when \\\gamma = 0.2\\ and \\F \approx 1.66\\.

``` r
# Create the forced damped pendulum system
pendulum <- forced_damped_pendulum(forcing = 1.66, damping = 0.2)
print(pendulum$description)
```

The system has three attractors located approximately at \\(0, 0)\\,
\\(2\pi, 0)\\, and \\(-2\pi, 0)\\.

#### Computing basins of attraction

``` r
# Compute basins at moderate resolution for the vignette
basins_pendulum <- compute_basins(
    pendulum,
    x_range = c(-pi, pi),
    y_range = c(-3, 3),
    resolution = 300,
    verbose = FALSE
)

print(basins_pendulum)
```

#### Visualizing the basins

``` r
# Custom color palette
pendulum_colors <- c("#E63946", "#457B9D", "#2A9D8F")

plot(basins_pendulum,
     title = "Forced damped pendulum (F = 1.66)",
     colors = pendulum_colors) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    )
```

#### Varying the forcing amplitude

Different forcing values produce different basin structures:

- \\F = 1.66\\: Full Wada basins (classic example)
- \\F = 1.71\\: Partially Wada basins
- \\F \> 2.0\\: Highly chaotic dynamics

``` r
# Partially Wada regime
pendulum_partial <- forced_damped_pendulum(forcing = 1.71)
basins_partial <- compute_basins(
    pendulum_partial,
    x_range = c(-pi, pi),
    y_range = c(-3, 3),
    resolution = 300,
    verbose = FALSE
)

plot(basins_partial,
     title = "Partially Wada basins (F = 1.71)",
     colors = c("#FF6B6B", "#4ECDC4", "#45B7D1")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
```

### The Henon-Heiles system

The Henon-Heiles system is a paradigmatic model in nonlinear dynamics,
originally introduced to study stellar motion in galaxies. The
Hamiltonian is:

\\H = \frac{1}{2}(p_x^2 + p_y^2) + \frac{1}{2}(x^2 + y^2) + x^2 y -
\frac{y^3}{3}\\

The potential has triangular symmetry with three saddle points at the
critical energy \\E_c = 1/6\\. Above this energy, particles can escape
through three channels located at 120-degree intervals.

``` r
# Create Henon-Heiles system above critical energy
hh <- henon_heiles_system(energy = 0.2)
print(hh$description)
```

#### Computing escape basins

``` r
basins_hh <- compute_basins(
    hh,
    x_range = c(-0.4, 0.4),
    y_range = c(-0.4, 0.4),
    resolution = 300,
    t_max = 200,
    verbose = FALSE
)
```

``` r
# Henon-Heiles color palette (vibrant)
hh_colors <- c("#F72585", "#7209B7", "#3A0CA3")

plot(basins_hh,
     title = "Henon-Heiles escape basins (E = 0.2)",
     colors = hh_colors) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
```

### Newton fractals

The Newton-Raphson method for solving \\z^n - 1 = 0\\ produces basins of
attraction for each root. For \\n \geq 3\\, these basins exhibit the
Wada property—this was proven mathematically by Przytycki (1989).

The Newton iteration is: \\z\_{k+1} = z_k - \frac{z_k^n - 1}{n
z_k^{n-1}}\\

The \\n\\ roots of unity are equally spaced on the unit circle:
\\\zeta_k = e^{2\pi i k / n}, \quad k = 0, 1, \ldots, n-1\\

#### Cubic Newton fractal (\\z^3 - 1 = 0\\)

``` r
basins_newton3 <- compute_newton_basins(
    n_roots = 3,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 500
)

# Newton fractal colors (complementary)
newton_colors <- c("#264653", "#E9C46A", "#E76F51")

plot(basins_newton3,
     title = expression(paste("Newton fractal: ", z^3 - 1 == 0)),
     colors = newton_colors) +
    theme_void() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
    )
```

#### Quintic Newton fractal (\\z^5 - 1 = 0\\)

``` r
basins_newton5 <- compute_newton_basins(
    n_roots = 5,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 500
)

# Five-color palette
newton5_colors <- c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6")

plot(basins_newton5,
     title = expression(paste("Newton fractal: ", z^5 - 1 == 0)),
     colors = newton5_colors) +
    theme_void() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
    )
```

#### Zooming into the boundary

The fractal nature of the Wada boundary becomes apparent when zooming
in:

``` r
basins_zoom <- compute_newton_basins(
    n_roots = 3,
    x_range = c(-0.5, 0.5),
    y_range = c(0.5, 1.5),
    resolution = 500
)

plot(basins_zoom,
     title = "Newton fractal: zoomed boundary",
     colors = c("#780000", "#C1121F", "#FDF0D5")) +
    theme_void() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
    )
```

------------------------------------------------------------------------

## Wada detection methods

The `wadaR` package implements three complementary methods for testing
the Wada property, each based on different mathematical principles.

### The grid method

#### Theoretical foundation

The grid method, introduced by Daza et al. (2015), is based on a key
topological observation: for Wada basins, between any two points of
different colors (basins), a third color can always be found at
sufficiently fine resolution.

**Algorithm overview:**

1.  Classify each grid cell by the number of distinct basins in its
    8-connected neighborhood
2.  For cells bordering exactly 2 basins, use bisection refinement to
    search for a third basin
3.  Compute the proportion \\W_m\\ of boundary cells bordering exactly
    \\m\\ basins

**Wada criterion:** The basins are Wada if \\W\_{N_A} \approx 1\\,
meaning all boundary cells border all \\N_A\\ basins.

#### Applying the grid method

``` r
# Apply grid method to forced damped pendulum
result_grid <- wada_grid_method(basins_pendulum, verbose = FALSE)
print(result_grid)
```

#### Visualizing boundary classification

``` r
plot(result_grid, basins = basins_pendulum) +
    scale_fill_manual(
        values = c("white", "#FEE0D2", "#FC9272", "#DE2D26"),
        name = "# basins",
        labels = c("Interior", "2", "3", "")
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
```

``` r
plot(result_grid, basins = basins_pendulum, show_wada_only = TRUE) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
```

#### Grid method on Newton fractal

The Newton fractal is mathematically proven to be Wada, making it an
ideal test case:

``` r
result_grid_newton <- wada_grid_method(basins_newton3, verbose = FALSE)
print(result_grid_newton)
```

As expected, \\W_3 \approx 1\\ confirms the Wada property.

### The merging method

#### Theoretical foundation

The merging method, introduced by Daza et al. (2018), exploits a
fundamental property of Wada basins: **merging basins does not change
their common boundary**.

For each basin \\B_i\\, define the merged basin \\M_i = \bigcup\_{j \neq
i} B_j\\. The boundary between \\B_i\\ and \\M_i\\ is called the *slim
boundary*. For Wada basins, all slim boundaries are identical:

\\\partial(B_1, M_1) = \partial(B_2, M_2) = \cdots = \partial(B\_{N_A},
M\_{N_A})\\

The method compares slim boundaries using the **Hausdorff distance**:

\\d_H(X, Y) = \max\left\\\sup\_{x \in X} d(x, Y), \\ \sup\_{y \in Y}
d(y, X)\right\\\\

**Wada criterion:** If all pairwise Hausdorff distances between slim
boundaries are small relative to the phase space size, the basins are
Wada.

#### Applying the merging method

``` r
# Apply merging method to forced damped pendulum
result_merging <- wada_merging_method(basins_pendulum, verbose = FALSE)
print(result_merging)
```

#### Visualizing slim boundaries

``` r
plot(result_merging, show_all = TRUE) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right"
    )
```

#### Comparing specific boundary pairs

``` r
plot(result_merging, show_all = FALSE, boundary1 = 1, boundary2 = 2) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
```

#### Merging method on Henon-Heiles system

``` r
result_merging_hh <- wada_merging_method(basins_hh, verbose = FALSE)
print(result_merging_hh)
```

``` r
plot(result_merging_hh) +
    scale_color_manual(values = c("#F72585", "#7209B7", "#3A0CA3")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
```

### The saddle-straddle method

#### Theoretical foundation

The saddle-straddle method, introduced by Wagemakers et al. (2020),
probes the dynamical structure of basin boundaries. For dissipative
systems, boundaries contain **chaotic saddles**—invariant sets where
trajectories exhibit transient chaotic behavior before escaping to an
attractor.

**Key insight:** Connected Wada basins are separated by a single chaotic
saddle. If the boundaries have the Wada property, then all merged basin
pairs share the same chaotic saddle:

\\S_1 = S_2 = \cdots = S\_{N_A}\\

**The straddle algorithm:**

1.  Find a segment straddling the boundary (endpoints in different
    basins)
2.  Refine using bisection until segment length \\\< \varepsilon\\
3.  Record the midpoint (approximately on the saddle)
4.  Iterate both endpoints forward (segment expands along unstable
    manifold)
5.  Repeat refinement and collection

**Wada criterion:** If the Hausdorff distance between all saddle
approximations is much smaller than their diameter, the basins are Wada.

#### Applying the saddle-straddle method

``` r
# Apply saddle-straddle method
# Using fewer points for vignette speed
result_straddle <- wada_straddle_method(
    pendulum,
    x_range = c(-pi, pi),
    y_range = c(-3, 3),
    n_points = 3000,
    verbose = FALSE
)
print(result_straddle)
```

#### Visualizing chaotic saddles

``` r
plot(result_straddle) +
    scale_color_manual(values = c("#FF595E", "#FFCA3A", "#8AC926")) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right"
    )
```

#### Saddles overlaid on basins

``` r
plot(result_straddle, basins = basins_pendulum) +
    scale_color_manual(values = c("#FFFFFF", "#FFD700", "#00FFFF")) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right"
    )
```

#### Individual saddle panels

``` r
plot(result_straddle, overlay = FALSE) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
```

------------------------------------------------------------------------

## Comparing detection methods

### Running all methods

``` r
# Summary comparison
cat("=== Wada detection results for forced damped pendulum (F = 1.66) ===\n\n")

cat("Grid method:\n")
cat(sprintf("  W_3 = %.4f\n", result_grid$W_m[3]))
cat(sprintf("  Is Wada: %s\n\n", result_grid$is_wada))

cat("Merging method:\n")
cat(sprintf("  Relative difference = %.4f\n", result_merging$relative_diff))
cat(sprintf("  Is Wada: %s\n\n", result_merging$is_wada))

cat("Saddle-straddle method:\n")
cat(sprintf("  Max Hausdorff / diameter = %.6f\n",
            result_straddle$max_distance / max(result_straddle$diameters)))
cat(sprintf("  Is Wada: %s\n", result_straddle$is_wada))
```

### Comparing full vs. partial Wada

``` r
# Compare F = 1.66 (full Wada) vs F = 1.71 (partial Wada)
grid_partial <- wada_grid_method(basins_partial, verbose = FALSE)
merging_partial <- wada_merging_method(basins_partial, verbose = FALSE)

comparison <- data.frame(
    Method = c("Grid", "Merging"),
    `F = 1.66 (Wada)` = c(
        sprintf("W_3 = %.3f", result_grid$W_m[3]),
        sprintf("rel_diff = %.3f", result_merging$relative_diff)
    ),
    `F = 1.71 (Partial)` = c(
        sprintf("W_3 = %.3f", grid_partial$W_m[3]),
        sprintf("rel_diff = %.3f", merging_partial$relative_diff)
    ),
    check.names = FALSE
)

knitr::kable(comparison, caption = "Comparison of detection methods for full vs. partial Wada basins")
```

### Method recommendations

| Method              | Best for           | Speed  | Accuracy |
|---------------------|--------------------|--------|----------|
| **Grid**            | Final confirmation | Medium | High     |
| **Merging**         | Quick screening    | Fast   | Medium   |
| **Saddle-straddle** | Dynamical insight  | Slow   | High     |

**Recommended workflow:**

1.  Start with the **merging method** for quick screening (no dynamics
    integration required)
2.  Confirm with the **grid method** for higher accuracy
3.  Use **saddle-straddle** when you need to understand the chaotic
    structure

------------------------------------------------------------------------

## Auxiliary functions

### Extracting basin boundaries

The
[`get_boundary()`](https://robustecologies.github.io/wadaR/reference/get_boundary.md)
function extracts boundary points from a basin matrix:

``` r
boundary <- get_boundary(basins_pendulum)
head(boundary)
```

``` r
ggplot(boundary, aes(x = x, y = y, color = factor(basin))) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = pendulum_colors, name = "Basin") +
    coord_fixed() +
    labs(title = "Basin boundaries", x = "x", y = "y") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
```

### Merging basins

The
[`merge_basins()`](https://robustecologies.github.io/wadaR/reference/merge_basins.md)
function creates a two-color basin map by merging all basins except one:

``` r
merged <- merge_basins(basins_pendulum$basins, keep_basin = 1)

# Convert to data frame for plotting
merged_df <- expand.grid(
    x = basins_pendulum$x_grid,
    y = basins_pendulum$y_grid
)
merged_df$basin <- as.vector(merged)

ggplot(merged_df, aes(x = x, y = y, fill = factor(basin))) +
    geom_raster() +
    scale_fill_manual(
        values = c("#E63946", "#457B9D"),
        name = "Basin",
        labels = c("B_1", "Others")
    ) +
    coord_fixed() +
    labs(title = "Merged basins: B_1 vs. M_1", x = "x", y = "y") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
```

### Computing Hausdorff distance

The
[`hausdorff_distance()`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md)
function computes the Hausdorff metric between point sets:

``` r
# Example: distance between two boundary segments
X <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
Y <- matrix(c(0.1, 0, 1.1, 0, 0.1, 1), ncol = 2, byrow = TRUE)

d <- hausdorff_distance(X, Y)
cat(sprintf("Hausdorff distance: %.4f\n", d))
```

------------------------------------------------------------------------

## Mathematical foundations

### Fractal dimension of Wada boundaries

Wada boundaries are typically fractals with dimension \\D\\ satisfying
\\1 \< D \< 2\\. The box-counting dimension can be estimated as:

\\D = \lim\_{\varepsilon \to 0} \frac{\log
N(\varepsilon)}{\log(1/\varepsilon)}\\

where \\N(\varepsilon)\\ is the number of \\\varepsilon\\-boxes covering
the boundary.

### Connection to invariant manifolds

For dissipative systems, the Wada boundary is related to the stable
manifold \\W^s(P)\\ of an unstable periodic orbit \\P\\:

**Theorem (Nusse-Yorke):** If \\P\\ satisfies certain accessibility
conditions, then \\W^s(P)\\ is dense in the Wada boundary.

### Basin entropy

The basin entropy quantifies unpredictability:

\\S_b = \lim\_{\varepsilon \to 0} \frac{1}{N(\varepsilon)}
\sum\_{i=1}^{N(\varepsilon)} H(p_i)\\

where \\H(p_i) = -\sum_j p\_{ij} \log p\_{ij}\\ is the Shannon entropy
of basin probabilities in box \\i\\.

For Wada basins, \\S_b\\ approaches its maximum value \\S_b^{\max} =
\log N_A\\.

``` r
# Compute basin entropy for the pendulum system
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3),
                         resolution = 200, verbose = FALSE)

# Compute entropy with different box sizes
e1 <- basin_entropy(basins, box_size = 10)
e2 <- basin_entropy(basins, box_size = 20)

cat(sprintf("Box size 10: S_b = %.3f (normalized = %.3f)\n",
            e1$S_b, e1$S_normalized))
cat(sprintf("Box size 20: S_b = %.3f (normalized = %.3f)\n",
            e2$S_b, e2$S_normalized))
cat(sprintf("Maximum entropy S_max = log2(3) = %.3f\n", log2(3)))

# Visualize entropy distribution
plot(e1) + labs(title = "Local basin entropy")
```

------------------------------------------------------------------------

## Computational complexity

Understanding the computational cost of each function helps choose
appropriate parameters:

| Function | Complexity | Notes |
|----|----|----|
| [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md) | \\O(N^2 \cdot T/\Delta t)\\ | Each of \\N^2\\ grid points requires ODE integration |
| [`compute_newton_basins()`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md) | \\O(N^2 \cdot k)\\ | Much faster, no ODE, just Newton iteration |
| [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md) | \\O(B \cdot 2^d)\\ | \\B\\ boundary boxes, \\d\\ = `max_refinements` |
| [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md) | \\O(N_A \cdot B^2)\\ | Hausdorff distance is \\O(B^2)\\ without k-d tree |
| [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md) | \\O(N_A \cdot P \cdot T/\Delta t)\\ | Most expensive: full ODE integration |
| [`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md) | \\O(N^2 / s^2)\\ | Fast, just counting over \\s \times s\\ boxes |

**Recommendations:**

- Use resolution 100-200 for exploration, 400-500 for publication
- [`compute_newton_basins()`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md)
  is fastest for mathematical examples
- Start with
  [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
  for quick Wada screening
- Use
  [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
  only when dynamical insight is needed

------------------------------------------------------------------------

## References

1.  Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
    Nonlinear Phenomena*, 51(1-3), 213-225.
    [doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789(91)90234-Z)

2.  Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and basin
    cells. *Physica D: Nonlinear Phenomena*, 90(3), 242-261.
    [doi:10.1016/0167-2789(95)00249-9](https://doi.org/10.1016/0167-2789(95)00249-9)

3.  Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
    Testing for basins of Wada. *Scientific Reports*, 5, 16579.
    [doi:10.1038/srep16579](https://doi.org/10.1038/srep16579)

4.  Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining
    when a basin is Wada: The merging method. *Scientific Reports*,
    8, 9954.
    [doi:10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

5.  Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2020). The
    saddle-straddle method to test for Wada basins. *Communications in
    Nonlinear Science and Numerical Simulation*, 84, 105167.
    [doi:10.1016/j.cnsns.2020.105167](https://doi.org/10.1016/j.cnsns.2020.105167)

6.  Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2021). How to detect
    Wada basins. *Discrete and Continuous Dynamical Systems - B*, 26(1),
    717-739.
    [doi:10.3934/dcdsb.2020330](https://doi.org/10.3934/dcdsb.2020330)

7.  Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., &
    Sanjuan, M. A. F. (2016). Basin entropy: A new tool to analyze
    uncertainty in dynamical systems. *Scientific Reports*, 6, 31416.
    [doi:10.1038/srep31416](https://doi.org/10.1038/srep31416)

8.  Aguirre, J., Vallejo, J. C., & Sanjuan, M. A. F. (2001). Wada basins
    and chaotic invariant sets in the Henon-Heiles system. *Physical
    Review E*, 64(6), 066208.
    [doi:10.1103/PhysRevE.64.066208](https://doi.org/10.1103/PhysRevE.64.066208)

9.  Henon, M., & Heiles, C. (1964). The applicability of the third
    integral of motion: Some numerical experiments. *The Astronomical
    Journal*, 69, 73.
    [doi:10.1086/109234](https://doi.org/10.1086/109234)

10. Przytycki, F. (1989). Remarks on the simple connectedness of basins
    of sinks for iterations of rational maps. *Banach Center
    Publications*, 23(1), 229-235.
    [doi:10.4064/-23-1-229-235](https://doi.org/10.4064/-23-1-229-235)

11. Curry, J. H., Garnett, L., & Sullivan, D. (1983). On the iteration
    of a rational function: Computer experiments with Newton’s method.
    *Communications in Mathematical Physics*, 91(2), 267-277.
    [doi:10.1007/BF01211162](https://doi.org/10.1007/BF01211162)

12. Seoane, J. M., & Sanjuan, M. A. F. (2013). New developments in
    classical chaotic scattering. *Reports on Progress in Physics*,
    76(1), 016001.
    [doi:10.1088/0034-4885/76/1/016001](https://doi.org/10.1088/0034-4885/76/1/016001)

13. Aguirre, J., Viana, R. L., & Sanjuan, M. A. F. (2009). Fractal
    structures in nonlinear dynamics. *Reviews of Modern Physics*,
    81(1), 333-386.
    [doi:10.1103/RevModPhys.81.333](https://doi.org/10.1103/RevModPhys.81.333)

------------------------------------------------------------------------

## Session information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8    
#>  [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8   
#>  [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Madrid
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.1 wadaR_1.0.0  
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6       jsonlite_2.0.0     dplyr_1.1.4        compiler_4.5.2    
#>  [5] tidyselect_1.2.1   Rcpp_1.1.1         dichromat_2.0-0.1  jquerylib_0.1.4   
#>  [9] systemfonts_1.2.1  scales_1.4.0       textshaping_1.0.3  yaml_2.3.12       
#> [13] fastmap_1.2.0      R6_2.6.1           generics_0.1.4     knitr_1.51        
#> [17] htmlwidgets_1.6.4  tibble_3.3.0       desc_1.4.3         bslib_0.9.0       
#> [21] pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.1.7        cachem_1.1.0      
#> [25] xfun_0.55          fs_1.6.6           sass_0.4.10        S7_0.2.1          
#> [29] otel_0.2.0         cli_3.6.5          withr_3.0.2        pkgdown_2.2.0     
#> [33] magrittr_2.0.4     digest_0.6.39      grid_4.5.2         rstudioapi_0.17.1 
#> [37] lifecycle_1.0.5    vctrs_0.6.5        evaluate_1.0.5     glue_1.8.0        
#> [41] farver_2.1.2       ragg_1.2.7         rmarkdown_2.30     tools_4.5.2       
#> [45] pkgconfig_2.0.3    htmltools_0.5.9
```
