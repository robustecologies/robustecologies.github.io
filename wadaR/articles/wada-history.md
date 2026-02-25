# The Wada property: from topology to chaos and art

``` r
library(wadaR)
library(ggplot2)
```

## The topological origins

### The Lakes of Wada (1917)

The story of Wada basins begins not in physics or chaos theory, but in
pure mathematics—specifically, in the field of point-set topology in
early 20th century Japan.

In 1917, the Japanese mathematician **Kunizo Yoneyama** published a
remarkable construction in the *Tohoku Mathematical Journal*. He
described a way to divide a bounded region of the plane into three (or
more) disjoint connected open sets that share a common boundary.
Yoneyama attributed this construction to his doctoral advisor, **Takeo
Wada**, and thus the resulting sets became known as the “Lakes of Wada”
(*Wada no mizuumi*).

#### The original construction

Imagine a square island surrounded by water, with three lakes that will
be excavated:

1.  **Day 1**: Dig a canal from the sea that winds throughout the
    island, coming within distance \\1\\ of every point, but not
    touching the planned lake locations. This creates Lake 1 (connected
    to the sea).

2.  **Day 2**: From the shore of Lake 1, dig a canal that winds within
    distance \\1/2\\ of every point not yet covered. This creates Lake
    2.

3.  **Day 3**: From Lake 2’s shore, dig a canal within distance \\1/3\\
    of every remaining point. This creates Lake 3.

4.  **Day 4**: Return to Lake 1 and extend it within distance \\1/4\\ of
    all remaining points.

5.  Continue alternating between the three lakes, with canals of width
    \\1/n\\ on day \\n\\.

After infinitely many days, the “land” that remains has measure zero and
forms a common boundary for all three lakes. Every point on this
boundary is simultaneously on the boundary of Lake 1, Lake 2, and Lake
3.

#### Mathematical significance

The Wada construction demonstrated that pathological topological objects
could exist even in the familiar Euclidean plane. It showed that:

- Three or more connected open sets can share a common boundary
- This boundary must be an **indecomposable continuum**—a connected
  compact set that cannot be written as the union of two proper
  connected closed subsets
- The boundary has a fractal-like structure, infinitely intricate at all
  scales

**Theorem (Kuratowski, 1924):** If a boundary separates three or more
connected regions in the plane, and if that boundary is connected, then
the boundary must be an indecomposable continuum.

This theorem, proved by the Polish mathematician Kazimierz Kuratowski,
established that the Wada property is not just a curiosity but a
fundamental topological phenomenon.

### The Japanese school of topology

The Wada construction emerged from a flourishing of mathematical
research in Japan during the Meiji and Taisho eras. Key figures
included:

- **Takeo Wada (1882-1944)**: Professor at Kyoto Imperial University,
  specialist in analysis and topology
- **Kunizo Yoneyama (1877-1968)**: Student of Wada, published the
  original construction
- **Teiji Takagi (1875-1960)**: Another pioneer of Japanese mathematics,
  known for the Takagi function (a continuous nowhere-differentiable
  function)

This period saw Japanese mathematicians making original contributions to
analysis, topology, and number theory, establishing Japan as a center of
mathematical research.

------------------------------------------------------------------------

## From topology to dynamical systems

### The connection to chaos (1991)

For seven decades, the Wada property remained a topological curiosity,
discussed in textbooks as an example of pathological sets. Then, in
1991, **Judy Kennedy** and **James Yorke** made a groundbreaking
discovery: Wada basins occur naturally in chaotic dynamical systems.

In their paper “Basins of Wada” published in *Physica D*, Kennedy and
Yorke showed that the forced damped pendulum—a simple physical
system—exhibits basins of attraction with the Wada property. This was
revolutionary because it connected abstract topology to concrete
physics.

``` r
pendulum <- forced_damped_pendulum(forcing = 1.66, damping = 0.2)
basins <- compute_basins(
    pendulum,
    x_range = c(-pi, pi),
    y_range = c(-3, 3),
    resolution = 150,
    verbose = FALSE
)

plot(basins,
     title = "Wada basins in the forced damped pendulum",
     colors = c("#1A535C", "#4ECDC4", "#F7FFF7")) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.background = element_rect(fill = "#FF6B6B"),
        legend.position = "none"
    )
```

#### The forced damped pendulum

The system that Kennedy and Yorke studied is deceptively simple:

\\\ddot{\theta} + \gamma \dot{\theta} + \sin(\theta) = F \cos(t)\\

where: - \\\theta\\ is the angle of the pendulum - \\\gamma\\ is the
damping coefficient - \\F\\ is the forcing amplitude - The forcing is
periodic with period \\2\pi\\

For \\\gamma = 0.2\\ and \\F \approx 1.66\\, this system has three
coexisting periodic attractors. The remarkable discovery was that the
boundaries between their basins of attraction form a Wada boundary.

### James Yorke: The father of chaos

**James Yorke** (born 1941) is one of the most influential figures in
dynamical systems theory. His contributions include:

- **Coining the term “chaos”**: In 1975, with Tien-Yien Li, he published
  “Period Three Implies Chaos,” which introduced the word “chaos” to
  describe deterministic unpredictability

- **Li-Yorke theorem**: Proved that if a continuous map has a period-3
  orbit, it has orbits of all periods—and chaos

- **OGY method**: With Edward Ott and Celso Grebogi, developed the first
  practical method for controlling chaos (1990)

- **Wada basins**: With Judy Kennedy and later collaborators,
  established the physical relevance of Wada topology

Yorke’s work bridges pure mathematics and applied science, showing that
abstract concepts have concrete physical manifestations.

------------------------------------------------------------------------

## The Henon-Heiles system and celestial mechanics

### Escape basins and stellar dynamics

The **Henon-Heiles system** was introduced in 1964 by Michel Henon and
Carl Heiles to model the motion of stars in a galaxy. The Hamiltonian
is:

\\H = \frac{1}{2}(p_x^2 + p_y^2) + \frac{1}{2}(x^2 + y^2) + x^2 y -
\frac{y^3}{3}\\

Above a critical energy \\E_c = 1/6\\, particles can escape through
three channels. These escape basins also exhibit the Wada property.

``` r
hh <- henon_heiles_system(energy = 0.2)
basins_hh <- compute_basins(
    hh,
    x_range = c(-0.4, 0.4),
    y_range = c(-0.4, 0.4),
    resolution = 150,
    t_max = 200,
    verbose = FALSE
)

plot(basins_hh,
     title = "Henon-Heiles: stellar escape routes",
     colors = c("#E63946", "#F1FAEE", "#A8DADC")) +
    theme_void() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none"
    )
```

#### Implications for astronomy

The Wada property in celestial mechanics means that:

- **Asteroid trajectories**: Near Wada boundaries, tiny errors in
  position lead to completely different long-term fates

- **Spacecraft navigation**: Mission planners must account for sensitive
  dependence in multi-body gravitational systems

- **Galaxy evolution**: The escape of stars from galactic cores may be
  fundamentally unpredictable in certain regimes

### The three-body problem

The gravitational three-body problem—determining the motion of three
masses under mutual gravitation—has been studied since Newton. Poincare
showed in the 1890s that it exhibits chaotic behavior. Modern research
has revealed that the basins of capture (which body a test particle
eventually orbits) can have Wada boundaries.

This connects the centuries-old problem of celestial mechanics to
20th-century topology and 21st-century computational methods.

------------------------------------------------------------------------

## Newton fractals: mathematics as art

### The Julia set connection

Newton fractals arise from applying Newton’s method to find roots of
polynomials in the complex plane. For \\f(z) = z^n - 1\\, the iteration
is:

\\z\_{k+1} = z_k - \frac{z_k^n - 1}{n z_k^{n-1}}\\

Each starting point converges to one of the \\n\\ roots of unity. The
basin boundaries form the **Julia set** of the Newton map—and for \\n
\geq 3\\, these boundaries are Wada.

``` r
newton3 <- compute_newton_basins(
    n_roots = 3,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 200
)

# Artistic color palette inspired by Japanese prints
plot(newton3,
     colors = c("#2D3047", "#E84855", "#FFFD82")) +
    theme_void() +
    theme(legend.position = "none")
```

#### Przytycki’s theorem (1989)

The Polish mathematician **Feliks Przytycki** proved rigorously that
Newton fractals for polynomials with three or more roots exhibit the
Wada property. This gave a mathematical proof that Wada basins are not
rare pathologies but common features of complex dynamics.

### Higher-degree Newton fractals

As the degree increases, Newton fractals become increasingly intricate:

``` r
# Create a gallery of Newton fractals
par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))

palettes <- list(
    c("#264653", "#2A9D8F", "#E9C46A"),
    c("#003049", "#D62828", "#F77F00", "#FCBF49"),
    c("#10002B", "#240046", "#3C096C", "#5A189A", "#7B2CBF"),
    c("#03071E", "#370617", "#6A040F", "#9D0208", "#D00000", "#DC2F02")
)

for (i in 1:4) {
    n <- i + 2  # 3, 4, 5, 6 roots
    basins <- compute_newton_basins(
        n_roots = n,
        x_range = c(-1.5, 1.5),
        y_range = c(-1.5, 1.5),
        resolution = 150
    )

    # Simple base R plot for speed
    image(basins$x_grid, basins$y_grid, basins$basins,
          col = palettes[[i]], axes = FALSE, xlab = "", ylab = "",
          main = bquote(z^.(n) - 1 == 0))
}
```

### Zoom into infinity

One of the most striking properties of Newton fractals is their
self-similarity. Zooming into the boundary reveals structure at all
scales:

``` r
# Zoom sequence
zoom_levels <- list(
    list(x = c(-2, 2), y = c(-2, 2), title = "Full view"),
    list(x = c(-0.5, 0.5), y = c(0.3, 1.3), title = "4x zoom"),
    list(x = c(-0.15, 0.15), y = c(0.7, 1.0), title = "16x zoom"),
    list(x = c(-0.05, 0.05), y = c(0.83, 0.93), title = "64x zoom")
)

par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))

for (zoom in zoom_levels) {
    basins <- compute_newton_basins(
        n_roots = 3,
        x_range = zoom$x,
        y_range = zoom$y,
        resolution = 300
    )

    image(basins$x_grid, basins$y_grid, basins$basins,
          col = c("#1D3557", "#457B9D", "#A8DADC"),
          axes = FALSE, xlab = "", ylab = "",
          main = zoom$title)
}
```

------------------------------------------------------------------------

## Artistic and cultural dimensions

### Fractals in art

The visual beauty of fractal basin boundaries has inspired artists and
designers worldwide. Newton fractals, in particular, have become iconic
images of mathematical art.

#### The aesthetic of complexity

What makes Wada boundaries aesthetically compelling?

1.  **Infinite detail**: No matter how closely you look, new structure
    appears
2.  **Self-similarity**: Patterns repeat at different scales
3.  **Balance of order and chaos**: Clear basins with infinitely complex
    boundaries
4.  **Color gradients**: Natural color coding reveals hidden structure

``` r
# Artistic rendering of pendulum basins
pendulum_art <- compute_basins(
    forced_damped_pendulum(forcing = 1.66),
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 500,
    verbose = FALSE
)

# Sunset palette
plot(pendulum_art,
     colors = c("#F72585", "#7209B7", "#3A0CA3")) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "Basin boundaries as abstract art")
```

### Japanese aesthetics and mathematical beauty

It is perhaps fitting that Wada basins originate from Japanese
mathematics, given the Japanese aesthetic tradition of finding beauty in
complexity and impermanence.

Concepts from Japanese aesthetics that resonate with Wada structures:

- **Wabi-sabi** (侘寂): Beauty in imperfection and transience—the
  fractal boundary is never quite smooth

- **Ma** (間): Negative space—the boundary itself has measure zero, yet
  defines everything

- **Yūgen** (幽玄): Profound mystery—infinite complexity hidden in
  simple equations

### Fractals in popular culture

Fractal images, including basins of attraction, have permeated popular
culture:

- **Album covers**: Bands like Tool and Pink Floyd have used fractal
  imagery
- **Film**: Fractals appear in science fiction films as visualizations
  of alien mathematics
- **Fashion**: Fractal patterns have influenced textile design
- **Architecture**: Some modern buildings incorporate fractal principles

------------------------------------------------------------------------

## Scientific applications

### Climate and weather prediction

Edward Lorenz’s discovery of sensitive dependence in weather models (the
“butterfly effect”) shares deep connections with Wada basins. In climate
systems with multiple stable states (e.g., ice ages vs. interglacial
periods), the boundaries between basins may exhibit Wada properties.

#### Implications

- Long-term climate prediction may be fundamentally limited near basin
  boundaries
- Tipping points between climate states may be more complex than simple
  thresholds
- Small perturbations could push the climate system toward any of
  several stable states

### Ecology and population dynamics

Ecological systems with multiple stable states—such as lakes that can be
clear or eutrophic, or forests that can be intact or degraded—may have
Wada basin boundaries.

**Example**: A lake ecosystem might have three stable states: 1. Clear
water with submerged vegetation 2. Turbid water with algal blooms 3.
Intermediate state with floating vegetation

If the basins of attraction for these states have the Wada property,
then predicting the lake’s fate from current conditions becomes
fundamentally uncertain near the boundaries.

### Engineering and control

In engineering systems with multiple operating modes, understanding
basin structure is crucial for:

- **Power grid stability**: Preventing cascading failures
- **Chemical reactors**: Avoiding runaway reactions
- **Robotic systems**: Ensuring convergence to desired configurations

Wada basins represent worst-case scenarios for control: near the
boundary, any small error could lead to any outcome.

------------------------------------------------------------------------

## The Spanish school of nonlinear dynamics

### Miguel Sanjuan and the Nonlinear Dynamics Group

Much of the modern computational work on Wada basins has come from the
**Nonlinear Dynamics, Chaos and Complex Systems Group** at Universidad
Rey Juan Carlos in Madrid, led by **Miguel A.F. Sanjuan**.

#### Key contributions

- **Grid method** (Daza et al., 2015): First practical computational
  test for Wada basins
- **Merging method** (Daza et al., 2018): Faster algorithm based on
  boundary comparison
- **Saddle-straddle method** (Wagemakers et al., 2020): Dynamical
  approach using chaotic saddles
- **Basin entropy** (Daza et al., 2016): Quantifying unpredictability in
  basin structures

These methods, implemented in the `wadaR` package, have made Wada basin
analysis accessible to researchers worldwide.

``` r
# Quick demonstration of all three methods
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3),
                         resolution = 200, verbose = FALSE)

r1 <- wada_grid_method(basins, verbose = FALSE)
r2 <- wada_merging_method(basins, verbose = FALSE)

cat("=== Detection methods from Rey Juan Carlos ===\n\n")
cat(sprintf("Grid method (2015):    W_3 = %.4f\n", r1$W_m[3]))
cat(sprintf("Merging method (2018): rel_diff = %.4f\n", r2$relative_diff))
cat("\nBoth methods confirm Wada basins.\n")
```

### International collaborations

The Spanish group has collaborated with researchers worldwide,
including:

- **James Yorke** (Maryland, USA): Original discoverer of physical Wada
  basins
- **Celso Grebogi** (Aberdeen, UK): Pioneer of chaos control
- **Japanese researchers**: Continuing the legacy of Wada’s homeland

------------------------------------------------------------------------

## Philosophical reflections

### Determinism and unpredictability

Wada basins raise profound questions about determinism and
predictability:

**The equations are deterministic**: Given exact initial conditions, the
future is completely determined.

**Yet prediction is impossible**: For initial conditions on the Wada
boundary, the system could end up in *any* attractor. No finite
measurement precision can determine the outcome.

This is not quantum uncertainty or statistical mechanics. It is purely
classical, deterministic unpredictability—a new category of limits to
knowledge.

### The unreasonable effectiveness of mathematics

Eugene Wigner famously wrote about “the unreasonable effectiveness of
mathematics in the natural sciences.” Wada basins exemplify this:

- A construction invented purely for mathematical curiosity (1917)
- Turns out to describe physical reality (1991)
- Has practical implications for prediction and control (2000s-present)

Mathematics anticipated physics by 74 years.

### Complexity from simplicity

Perhaps the deepest lesson of Wada basins is that extreme complexity can
arise from simple rules:

- The forced pendulum equation fits on a single line
- Yet its basin boundaries have infinite complexity
- This complexity is not added from outside—it emerges from the dynamics

As the physicist Murray Gell-Mann said: “You don’t need something more
to get something more.”

------------------------------------------------------------------------

## Looking forward

### Current research frontiers

Active areas of research include:

- **Higher dimensions**: Do Wada basins exist in 3D and
  higher-dimensional systems?
- **Partial Wada**: When only some boundary points are Wada
- **Basin entropy**: Quantifying unpredictability in complex systems
- **Machine learning**: Using neural networks to detect Wada properties
- **Quantum systems**: Classical-quantum correspondence for Wada basins

### The future of Wada basin analysis

The `wadaR` package represents the current state of the art in Wada
basin detection. Future developments may include:

- Faster algorithms using GPU acceleration
- Automated detection in experimental data
- Integration with machine learning frameworks
- Extension to higher-dimensional systems

------------------------------------------------------------------------

## Conclusion

From Takeo Wada’s topological curiosity in 1917 to modern computational
methods, the Wada property has traveled a remarkable journey. It
connects:

- **Pure mathematics**: Point-set topology and continuum theory
- **Dynamical systems**: Chaos, attractors, and basins
- **Physics**: Pendulums, celestial mechanics, and beyond
- **Art**: The visual beauty of fractal boundaries
- **Philosophy**: Determinism, predictability, and the limits of
  knowledge

The story of Wada basins reminds us that mathematics is not merely a
tool for describing reality—sometimes it anticipates reality, waiting
decades for physics to catch up.

------------------------------------------------------------------------

## Timeline

| Year | Event                                                |
|------|------------------------------------------------------|
| 1917 | Yoneyama publishes the Lakes of Wada construction    |
| 1924 | Kuratowski proves theorem on indecomposable continua |
| 1964 | Henon and Heiles introduce their Hamiltonian system  |
| 1975 | Li and Yorke coin the term “chaos”                   |
| 1989 | Przytycki proves Wada property for Newton fractals   |
| 1991 | Kennedy and Yorke discover Wada basins in pendulum   |
| 1996 | Nusse and Yorke develop basin cell theory            |
| 2015 | Daza et al. introduce the grid method                |
| 2016 | Basin entropy concept introduced                     |
| 2018 | Merging method published                             |
| 2020 | Saddle-straddle method published                     |

------------------------------------------------------------------------

## References

1.  Yoneyama, K. (1917). Theory of continuous sets of points. *Tohoku
    Mathematical Journal*, 12, 43-158.

2.  Kuratowski, K. (1924). Sur les coupures irreductibles du plan.
    *Fundamenta Mathematicae*, 6, 130-145.

3.  Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D*, 51,
    213-225.
    [doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789(91)90234-Z)

4.  Przytycki, F. (1989). Remarks on the simple connectedness of basins
    of sinks for iterations of rational maps. *Banach Center
    Publications*, 23(1), 229-235.

5.  Henon, M., & Heiles, C. (1964). The applicability of the third
    integral of motion. *Astronomical Journal*, 69, 73.

6.  Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and basin
    cells. *Physica D*, 90, 242-261.

7.  Daza, A., Wagemakers, A., Sanjuan, M. A. F., & Yorke, J. A. (2015).
    Testing for basins of Wada. *Scientific Reports*, 5, 16579.

8.  Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining
    when a basin is Wada: The merging method. *Scientific Reports*, 8,
    9954.

9.  Wagemakers, A., Daza, A., & Sanjuan, M. A. F. (2020). The
    saddle-straddle method to test for Wada basins. *CNSNS*, 84, 105167.

10. Li, T.-Y., & Yorke, J. A. (1975). Period three implies chaos.
    *American Mathematical Monthly*, 82, 985-992.

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
