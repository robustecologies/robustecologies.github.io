# User-defined dynamical systems

``` r

library(wadaR)
library(ggplot2)
library(knitr)
library(patchwork)
```

## Introduction

While wadaR provides several built-in dynamical systems for studying
Wada basins (the forced damped pendulum, Henon-Heiles system, and Newton
fractals), researchers often need to analyze custom systems from their
own work. This vignette demonstrates how to define arbitrary dynamical
systems that integrate seamlessly with all wadaR methods.

The
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
function creates high-performance, OpenMP-parallelized system objects
that work with:

- [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
  for basin computation
- [`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md)
  for uncertainty quantification [\[1\]](#ref1)
- [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
  for topological Wada detection [\[2\]](#ref2)
- [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
  for the merging algorithm [\[3\]](#ref3)
- [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
  for saddle-straddle analysis [\[4\]](#ref4)
- [`shinywadaR()`](https://robustecologies.github.io/wadaR/reference/shinywadaR.md)
  for interactive exploration

**Key feature**: All user-defined systems are compiled to native C++
code with OpenMP parallelization, providing performance matching the
built-in systems.

In this vignette, we demonstrate the complete wadaR analysis workflow on
three classic nonlinear systems:

1.  **Duffing oscillator**: A paradigmatic example of chaos in driven
    oscillators [\[5\]](#ref5)
2.  **Henon map**: A discrete-time system exhibiting strange attractors
    [\[6\]](#ref6)
3.  **Van der Pol oscillator**: A self-sustained oscillator with limit
    cycle dynamics [\[7\]](#ref7)

  

## The compiled_system architecture

  

### How it works

When you call
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md):

1.  Your C++ dynamics code is embedded into a complete C++ source file
2.  An RK4 integrator (for ODEs) or map iterator is generated
3.  OpenMP parallelization is added for multi-core computation
4.  The code is compiled at runtime using
    [`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)
5.  The compiled function is stored in the system object

This means your custom systems run at full C++ speed with parallel
execution.

  

### C++ dynamics code format

  

#### For ODEs

Write code that computes derivatives and stores them in the `deriv[]`
array:

``` cpp
// Available variables:
// state[i] - current state (read)
// deriv[i] - derivatives to compute (write)
// t        - current time
// dim      - state dimension
// Parameters by name (e.g., damping, omega)

deriv[0] = state[1];                    // dx/dt = v
deriv[1] = -damping * state[1] - ...;   // dv/dt = ...
```

  

#### For discrete maps

Write code that computes the next state in `next_state[]`:

``` cpp
// Available variables:
// state[i]      - current state (read)
// next_state[i] - next iteration (write)
// iter          - current iteration number
// Parameters by name

next_state[0] = 1 - a * state[0] * state[0] + state[1];
next_state[1] = b * state[0];
```

  

### Attractor specifications

wadaR provides helper functions for different attractor types:

| Function | Use case | Required fields |
|:---|:---|:---|
| [`attractor_point()`](https://robustecologies.github.io/wadaR/reference/attractor_point.md) | Fixed points, periodic orbits | `center`, `radius` |
| [`attractor_cycle()`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md) | Limit cycles | `center`, `radius`, optional `period` |
| [`attractor_exit()`](https://robustecologies.github.io/wadaR/reference/attractor_exit.md) | Escape channels | `angle` |
| [`attractor_outcome()`](https://robustecologies.github.io/wadaR/reference/attractor_outcome.md) | Discrete outcomes | `id`, `label` |

Attractor specification functions in wadaR {.table}

------------------------------------------------------------------------

  

## The Duffing oscillator

  

### Mathematical background

The Duffing oscillator is one of the most studied examples in nonlinear
dynamics [\[5\]](#ref5). Named after Georg Duffing, who first analyzed
it in 1918, the equation describes a damped oscillator with a nonlinear
restoring force [\[8\]](#ref8):

\\\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \cos(\omega
t)\\

where:

- \\\delta\\ is the damping coefficient
- \\\alpha\\ controls the linear stiffness (negative for double-well
  potential)
- \\\beta\\ is the nonlinear stiffness parameter
- \\\gamma\\ is the forcing amplitude
- \\\omega\\ is the forcing frequency

Converting to a first-order system with \\v = \dot{x}\\:

\\\begin{aligned} \frac{dx}{dt} &= v \\ \frac{dv}{dt} &= -\delta v -
\alpha x - \beta x^3 + \gamma \cos(\omega t) \end{aligned}\\

The potential energy function is:

\\V(x) = \frac{\alpha}{2}x^2 + \frac{\beta}{4}x^4\\

When \\\alpha \< 0\\ and \\\beta \> 0\\, this creates a **double-well
potential** with stable equilibria at \\x = \pm\sqrt{-\alpha/\beta}\\
and an unstable saddle at \\x = 0\\. The presence of external forcing
can lead to chaotic behavior and fractal basin boundaries
[\[9\]](#ref9).

  

### Defining the system

``` r

duffing <- compiled_system(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -delta * state[1] - alpha * state[0]
                   - beta * state[0] * state[0] * state[0]
                   + gamma_f * cos(omega * t);
    ',
    attractors = list(
        attractor_point(center = c(1, 0), radius = 0.3, label = "Right well"),
        attractor_point(center = c(-1, 0), radius = 0.3, label = "Left well"),
        attractor_point(center = c(0, 0), radius = 0.3, label = "Origin")
    ),
    dim = 2,
    type = "ode",
    params = list(delta = 0.3, alpha = -1, beta = 1, gamma_f = 0.37, omega = 1.2),
    name = "Duffing oscillator",
    verbose = FALSE
)

print(duffing)
#> ⚙ Compiled dynamical system: Duffing oscillator
#> -------------------------------------------------- 
#>   Type:       ODE (C++/OpenMP compiled)
#>   Dimension:  2
#>   Attractors: 3
#>     [1] Right well at (1, 0) (r=0.300)
#>     [2] Left well at (-1, 0) (r=0.300)
#>     [3] Origin at (0, 0) (r=0.300)
#>   Parameters (compiled):
#>     delta = 0.3
#>     alpha = -1
#>     beta = 1
#>     gamma_f = 0.37
#>     omega = 1.2
#> -------------------------------------------------- 
#>   Duffing oscillator (compiled ODE in 2D)
#>   ✔ Ready for parallel computation
```

Note how parameters are passed as a named list. They become compile-time
constants in C++ for maximum performance.

  

#### Alternative: formula-based input via janos

When the companion package
[janos](https://github.com/robustecologies/janos) is installed, the same
Duffing system can be specified through formulas instead of a raw C++
string. wadaR translates the `model_spec` into its basin-kernel C++
template via the exported helper
[`janos::model_spec_rhs_cpp()`](https://robustecologies.github.io/janos/reference/model_spec_rhs_cpp.html);
the resulting `compiled_system` is numerically identical to the raw-C++
route but spares the user from manual array indexing and lets a model
already validated for
[`janos::dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.html)
flow directly into basin analysis.

``` r

duffing_model <- janos::model_spec(
    rhs         = list(
        x ~ y,
        y ~ -delta * y - alpha * x - beta * x * x * x +
            gamma_f * cos(omega * t)
    ),
    state_names = c("x", "y"),
    parms       = list(delta = 0.3, alpha = -1, beta = 1,
                       gamma_f = 0.37, omega = 1.2)
)

duffing_via_janos <- compiled_system(
    model      = duffing_model,
    attractors = list(
        attractor_point(c( 1, 0), 0.3, "Right well"),
        attractor_point(c(-1, 0), 0.3, "Left well"),
        attractor_point(c( 0, 0), 0.3, "Origin")
    ),
    name       = "Duffing oscillator (via janos::model_spec)",
    verbose    = FALSE
)

print(duffing_via_janos)
#> ⚙ Compiled dynamical system: Duffing oscillator (via janos::model_spec)
#> -------------------------------------------------- 
#>   Type:       ODE (C++/OpenMP compiled)
#>   Dimension:  2
#>   Attractors: 3
#>     [1] Right well at (1, 0) (r=0.300)
#>     [2] Left well at (-1, 0) (r=0.300)
#>     [3] Origin at (0, 0) (r=0.300)
#>   Parameters (compiled):
#>     delta = 0.3
#>     alpha = -1
#>     beta = 1
#>     gamma_f = 0.37
#>     omega = 1.2
#> -------------------------------------------------- 
#>   Adapted from janos::model_spec (janos 1.10.1)
#>   ✔ Ready for parallel computation
```

The two routes are interchangeable in downstream basin computations;
either object can be passed to
[`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md),
[`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md),
[`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md)
and the rest of the pipeline.

  

### Computing basins of attraction

``` r

basins_duffing <- compute_basins(
    duffing,
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 300,
    t_max = 100,
    dt = 0.01,
    verbose = TRUE
)
#> ⚙ Computing basins (compiled C++/OpenMP): Duffing oscillator
#>   Grid: 300x300 (90,000 points)
#>   Attractors: 3
#>   Type: ode (2D)
#>   Cores: 19
#>   Integration steps: ~900,000,000
#> ⏱ Completed in 0.28 seconds
```

  

### Visualization

``` r

plot(basins_duffing,
     title = "Duffing oscillator basins of attraction",
     show_boundary = TRUE,
     boundary_color = "white")
```

![Basins of attraction for the Duffing oscillator. Different colors
represent initial conditions that converge to different
attractors.](custom-systems_files/figure-html/duffing-plot-1.png)

Basins of attraction for the Duffing oscillator. Different colors
represent initial conditions that converge to different attractors.

  

### Basin entropy analysis

Basin entropy [\[1\]](#ref1) quantifies the unpredictability of the
final state given an initial condition. Higher entropy indicates more
complex, interleaved basin structures.

``` r

entropy_duffing <- basin_entropy(basins_duffing, box_size = 2)
```

| Metric | Value | Interpretation |
|:---|:---|:---|
| Basin entropy (\\S_b\\) | 0.1384 | Overall unpredictability |
| Boundary basin entropy (\\S\_{bb}\\) | 0.9650 | Uncertainty at boundaries |
| Boundary fraction | 14.35% | Fraction of phase space on boundaries |

Basin entropy metrics for the Duffing oscillator {.table}

``` r

plot(entropy_duffing)
```

![Entropy distribution across the Duffing oscillator phase
space.](custom-systems_files/figure-html/duffing-entropy-plot-1.png)

Entropy distribution across the Duffing oscillator phase space.

  

### Wada property detection

The Wada property [\[2\]](#ref2)[\[3\]](#ref3) states that a single
boundary separates all basins simultaneously. This topological property
implies extreme sensitivity to initial conditions.

  

#### Grid method

The grid method [\[10\]](#ref10) tests whether all basin boundaries
coincide by checking if boundary points neighbor all basin types:

``` r

grid_duffing <- wada_grid_method(basins_duffing, verbose = TRUE)
#> ¡ Basin matrix: 300 x 300, 3 attractors
#> ⚙ Grid method: 300x300 grid, 3 attractors
#> ¡ Classifying boundary boxes (parallel Rcpp, 19 cores)...
#> ⏱ Classification: 0.08 seconds
#> ¡ Refining 16065 boundary boxes with ODE bisection (parallel Rcpp, 19 cores)...
#> ¡ Press Esc to abort computation
#> ⏱ Refinement: 134.64 seconds
#> ¡ Results:
#>   G_2: 9214 boxes (W_2 = 0.4376)
#>   G_3: 11843 boxes (W_3 = 0.5624)
#> ⚠ Basin is PARTIALLY WADA or NOT WADA
```

| Result                     | Value  |
|:---------------------------|:-------|
| Wada detected              | No     |
| W-statistic (\\W\_{N_A}\\) | 0.5624 |
| Total boundary boxes       | 21,057 |

Grid method results for the Duffing oscillator {.table}

``` r

plot(grid_duffing)
```

![Grid method visualization showing boundary
classification.](custom-systems_files/figure-html/duffing-grid-plot-1.png)

Grid method visualization showing boundary classification.

  

#### Merging method

The merging method [\[3\]](#ref3) tests Wada by progressively merging
basins and checking if the boundary remains unchanged:

``` r

merge_duffing <- wada_merging_method(basins_duffing, verbose = TRUE)
#> ¡ Basin matrix: 300 x 300, 3 attractors
#> ⚙ Merging method: 300x300 grid, 3 attractors
#> ¡ Computing boundaries (parallel Rcpp)...
#>   ⚙ Boundaries [=============                           ]  33.3% | ETA: 0s  ⚙ Boundaries [===========================             ]  66.7% | ETA: 0s  ⚙ Boundaries [========================================] 100.0% | ETA: 0s
#>   ⚙ Boundaries [========================================] 100.0% | ETA: 0s
#> ⏱ Boundary extraction: 0.03 seconds
#> ¡ Computing Hausdorff distances (parallel k-d tree)...
#> ⏱ Distance computation: 2.87 seconds
#> ¡ Results:
#>   max_d = 0.7229
#>   min_d = 0.4447
#>   (max_d - min_d) / min_d = 0.6256
#>   max_d / phase_space_size = 0.1278
#> ⚠ Basin is PARTIALLY WADA or NOT WADA
```

| Metric                   | Value |
|:-------------------------|:------|
| Wada detected            | No    |
| Merged boundary fraction | NA    |
| Hausdorff distance (max) | -Inf  |
| Number of pairs tested   | 0     |

Merging method results for the Duffing oscillator {.table}

``` r

plot(merge_duffing)
```

![Merging method visualization showing pairwise boundary
comparisons.](custom-systems_files/figure-html/duffing-merge-plot-1.png)

Merging method visualization showing pairwise boundary comparisons.

------------------------------------------------------------------------

  

## The Henon map

  

### Mathematical background

The Henon map [\[6\]](#ref6) is a discrete-time dynamical system
introduced by Michel Henon in 1976 as a simplified model of the Poincare
section of the Lorenz system. It is defined by:

\\\begin{aligned} x\_{n+1} &= 1 - a x_n^2 + y_n \\ y\_{n+1} &= b x_n
\end{aligned}\\

where:

- \\a\\ controls the nonlinearity (typically \\a = 1.4\\)
- \\b\\ is the dissipation parameter (typically \\b = 0.3\\)

The Jacobian determinant is \\\det(J) = -b\\, making this a dissipative
map when \\\|b\| \< 1\\. For the classical parameters \\(a, b) = (1.4,
0.3)\\, the map exhibits a **strange attractor** with fractal dimension
approximately \\1.26\\ [\[11\]](#ref11).

The fixed points are found by solving \\x = 1 - ax^2 + bx\\ and \\y =
bx\\:

\\x^\* = \frac{(1-b) \pm \sqrt{(1-b)^2 + 4a}}{2a}\\

  

### Defining the system

``` r

henon <- compiled_system(
    cpp_dynamics = '
        next_state[0] = 1.0 - a * state[0] * state[0] + state[1];
        next_state[1] = b * state[0];
    ',
    attractors = list(
        attractor_point(center = c(0.63, 0.19), radius = 0.15,
                       label = "Strange attractor"),
        attractor_point(center = c(-1.13, -0.34), radius = 0.15,
                       label = "Fixed point"),
        attractor_point(center = c(0, 0), radius = 5.0,
                       label = "Escape")
    ),
    dim = 2,
    type = "map",
    params = list(a = 1.4, b = 0.3),
    name = "Henon map",
    verbose = FALSE
)

print(henon)
#> ⚙ Compiled dynamical system: Henon map
#> -------------------------------------------------- 
#>   Type:       Discrete map (C++/OpenMP compiled)
#>   Dimension:  2
#>   Attractors: 3
#>     [1] Strange attractor at (0.63, 0.19) (r=0.150)
#>     [2] Fixed point at (-1.13, -0.34) (r=0.150)
#>     [3] Escape at (0, 0) (r=5.000)
#>   Parameters (compiled):
#>     a = 1.4
#>     b = 0.3
#> -------------------------------------------------- 
#>   Henon map (compiled discrete map in 2D)
#>   ✔ Ready for parallel computation
```

  

### Computing basins

``` r

basins_henon <- compute_basins(
    henon,
    x_range = c(-1.5, 1.5),
    y_range = c(-0.5, 0.5),
    resolution = 400,
    t_max = 500,  # Maximum iterations for maps
    verbose = TRUE
)
#> ⚙ Computing basins (compiled C++/OpenMP): Henon map
#>   Grid: 400x400 (160,000 points)
#>   Attractors: 3
#>   Type: map (2D)
#>   Cores: 19
#>   Max iterations: 500
#> ⏱ Completed in 0.02 seconds
```

  

### Visualization

``` r

plot(basins_henon, title = "Henon map basins of attraction")
```

![Basins of attraction for the Henon map. The fractal structure reflects
the strange attractor's sensitivity to initial
conditions.](custom-systems_files/figure-html/henon-plot-1.png)

Basins of attraction for the Henon map. The fractal structure reflects
the strange attractor’s sensitivity to initial conditions.

  

### Basin entropy analysis

``` r

entropy_henon <- basin_entropy(basins_henon, box_size = 2)
```

| Metric                               | Value  |
|:-------------------------------------|:-------|
| Basin entropy (\\S_b\\)              | 0.0104 |
| Boundary basin entropy (\\S\_{bb}\\) | 0.9471 |
| Boundary fraction                    | 1.10%  |

Basin entropy metrics for the Henon map {.table}

``` r

plot(entropy_henon)
```

![Entropy distribution across the Henon map phase
space.](custom-systems_files/figure-html/henon-entropy-plot-1.png)

Entropy distribution across the Henon map phase space.

  

### Wada property detection

``` r

grid_henon <- wada_grid_method(basins_henon, verbose = TRUE)
#> ¡ Basin matrix: 400 x 400, 3 attractors
#> ⚙ Grid method: 400x400 grid, 3 attractors
#> ¡ Classifying boundary boxes (parallel Rcpp, 19 cores)...
#> ⏱ Classification: 0.01 seconds
#> ¡ Refining 3520 boundary boxes with ODE bisection (parallel Rcpp, 19 cores)...
#> ¡ Press Esc to abort computation
#> ⏱ Refinement: 0.02 seconds
#> ¡ Results:
#>   G_2: 3520 boxes (W_2 = 1.0000)
#>   G_3: 0 boxes (W_3 = 0.0000)
#> ⚠ Basin is PARTIALLY WADA or NOT WADA
```

| Result                     | Value  |
|:---------------------------|:-------|
| Wada detected              | No     |
| W-statistic (\\W\_{N_A}\\) | 0.0000 |
| Total boundary boxes       | 3,520  |

Grid method results for the Henon map {.table}

``` r

plot(grid_henon)
```

![Grid method results for the Henon
map.](custom-systems_files/figure-html/henon-grid-plot-1.png)

Grid method results for the Henon map.

------------------------------------------------------------------------

  

## The Van der Pol oscillator

  

### Mathematical background

The Van der Pol oscillator [\[7\]](#ref7)[\[12\]](#ref12) is a
non-conservative system with nonlinear damping, introduced by Balthasar
van der Pol in 1926 to describe oscillations in vacuum tube circuits.
The equation is:

\\\ddot{x} - \mu(1 - x^2)\dot{x} + x = 0\\

where \\\mu \> 0\\ is the nonlinearity parameter. The key feature is
that:

- When \\\|x\| \< 1\\: negative damping (energy input)
- When \\\|x\| \> 1\\: positive damping (energy dissipation)

This creates a **self-sustained oscillation** with a stable limit cycle.
The period of oscillation increases with \\\mu\\:

\\T \approx (3 - 2\ln 2)\mu + O(\mu^{-1/3}) \quad \text{for } \mu \gg
1\\

When external forcing is added [\[13\]](#ref13):

\\\ddot{x} - \mu(1 - x^2)\dot{x} + x = A\cos(\omega t)\\

the system can exhibit complex behavior including **frequency
entrainment**, **quasiperiodic motion**, and **chaos** [\[14\]](#ref14).

Converting to first-order form with \\v = \dot{x}\\:

\\\begin{aligned} \frac{dx}{dt} &= v \\ \frac{dv}{dt} &= \mu(1 - x^2)v -
x + A\cos(\omega t) \end{aligned}\\

  

### Defining the system

``` r

vdp <- compiled_system(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = mu * (1.0 - x*x) * v - x + A * cos(omega * t);
    ',
    attractors = list(
        attractor_cycle(center = c(2, 0), radius = 0.8,
                       period = 2*pi, label = "Limit cycle 1"),
        attractor_cycle(center = c(-2, 0), radius = 0.8,
                       period = 2*pi, label = "Limit cycle 2"),
        attractor_point(center = c(0, 0), radius = 0.5,
                       label = "Origin")
    ),
    dim = 2,
    type = "ode",
    params = list(mu = 0.8, A = 1.5, omega = 1.0),
    name = "Van der Pol oscillator",
    verbose = FALSE
)

print(vdp)
#> ⚙ Compiled dynamical system: Van der Pol oscillator
#> -------------------------------------------------- 
#>   Type:       ODE (C++/OpenMP compiled)
#>   Dimension:  2
#>   Attractors: 3
#>     [1] Limit cycle 1 at (2, 0) (r=0.800)
#>     [2] Limit cycle 2 at (-2, 0) (r=0.800)
#>     [3] Origin at (0, 0) (r=0.500)
#>   Parameters (compiled):
#>     mu = 0.8
#>     A = 1.5
#>     omega = 1
#> -------------------------------------------------- 
#>   Van der Pol oscillator (compiled ODE in 2D)
#>   ✔ Ready for parallel computation
```

  

### Computing basins

``` r

basins_vdp <- compute_basins(
    vdp,
    x_range = c(-4, 4),
    y_range = c(-4, 4),
    resolution = 300,
    t_max = 200,
    dt = 0.01,
    verbose = TRUE
)
#> ⚙ Computing basins (compiled C++/OpenMP): Van der Pol oscillator
#>   Grid: 300x300 (90,000 points)
#>   Attractors: 3
#>   Type: ode (2D)
#>   Cores: 19
#>   Integration steps: ~1,800,000,000
#> ⏱ Completed in 0.06 seconds
```

  

### Visualization

``` r

plot(basins_vdp,
     title = "Van der Pol oscillator basins of attraction")
```

![Basins of attraction for the forced Van der Pol
oscillator.](custom-systems_files/figure-html/vdp-plot-1.png)

Basins of attraction for the forced Van der Pol oscillator.

  

### Basin entropy analysis

``` r

entropy_vdp <- basin_entropy(basins_vdp, box_size = 2)
```

| Metric                               | Value  |
|:-------------------------------------|:-------|
| Basin entropy (\\S_b\\)              | 0.0291 |
| Boundary basin entropy (\\S\_{bb}\\) | 0.9388 |
| Boundary fraction                    | 3.10%  |

Basin entropy metrics for the Van der Pol oscillator {.table}

``` r

plot(entropy_vdp)
```

![Entropy distribution across the Van der Pol phase
space.](custom-systems_files/figure-html/vdp-entropy-plot-1.png)

Entropy distribution across the Van der Pol phase space.

  

### Wada property detection

``` r

grid_vdp <- wada_grid_method(basins_vdp, verbose = TRUE)
#> ¡ Basin matrix: 300 x 300, 3 attractors
#> ⚙ Grid method: 300x300 grid, 3 attractors
#> ¡ Classifying boundary boxes (parallel Rcpp, 19 cores)...
#> ⏱ Classification: 0.01 seconds
#> ¡ Refining 5108 boundary boxes with ODE bisection (parallel Rcpp, 19 cores)...
#> ¡ Press Esc to abort computation
#> ⏱ Refinement: 47.37 seconds
#> ¡ Results:
#>   G_2: 3423 boxes (W_2 = 0.6396)
#>   G_3: 1929 boxes (W_3 = 0.3604)
#> ⚠ Basin is PARTIALLY WADA or NOT WADA
```

| Result                     | Value  |
|:---------------------------|:-------|
| Wada detected              | No     |
| W-statistic (\\W\_{N_A}\\) | 0.3604 |
| Total boundary boxes       | 5,352  |

Grid method results for the Van der Pol oscillator {.table}

``` r

plot(grid_vdp)
```

![Grid method results for the Van der Pol
oscillator.](custom-systems_files/figure-html/vdp-grid-plot-1.png)

Grid method results for the Van der Pol oscillator.

``` r

merge_vdp <- wada_merging_method(basins_vdp, verbose = TRUE)
#> ¡ Basin matrix: 300 x 300, 3 attractors
#> ⚙ Merging method: 300x300 grid, 3 attractors
#> ¡ Computing boundaries (parallel Rcpp)...
#>   ⚙ Boundaries [=============                           ]  33.3% | ETA: 0s  ⚙ Boundaries [===========================             ]  66.7% | ETA: 0s  ⚙ Boundaries [========================================] 100.0% | ETA: 0s
#>   ⚙ Boundaries [========================================] 100.0% | ETA: 0s
#> ⏱ Boundary extraction: 0.01 seconds
#> ¡ Computing Hausdorff distances (parallel k-d tree)...
#> ⏱ Distance computation: 0.41 seconds
#> ¡ Results:
#>   max_d = 4.0683
#>   min_d = 0.9550
#>   (max_d - min_d) / min_d = 3.2600
#>   max_d / phase_space_size = 0.3596
#> ⚠ Basin is PARTIALLY WADA or NOT WADA
```

``` r

plot(merge_vdp)
```

![Merging method results for the Van der Pol
oscillator.](custom-systems_files/figure-html/vdp-merge-plot-1.png)

Merging method results for the Van der Pol oscillator.

------------------------------------------------------------------------

  

## Automatic attractor detection

If you don’t know where the attractors are located, wadaR provides the
[`detect_attractors()`](https://robustecologies.github.io/wadaR/reference/detect_attractors.md)
function, which uses compiled C++/OpenMP trajectory integration to
automatically identify attractors.

``` r

# Detect attractors in a triple-well system
detected <- detect_attractors(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = -0.5 * v - state[0] * (state[0]*state[0] - 4);
    ',
    params = list(),
    dim = 2,
    x_range = c(-3, 3),
    y_range = c(-3, 3),
    n_trajectories = 300,
    t_max = 100,
    min_attractors = 3,
    max_attractors = 5,
    verbose = TRUE
)
#> ⚙ Detecting attractors (C++/OpenMP)
#>   Trajectories: 300
#>   Integration time: 100 (transient: 50)
#>   Clustering: hclust
#> ⚙ Compiling trajectory integrator...
#> ✔ Compilation successful
#> ⚙ Integrating trajectories...
#> ⏱ Integration completed in 0.05s
#> ¡ 300 valid trajectories for clustering
#> ⚙ Clustering endpoints...
#> ✔ Found 3 attractor clusters
#>     [1] point: center=(-2.000, -0.000), r=0.100, n=154
#>     [2] point: center=(2.000, -0.000), r=0.100, n=59
#>     [3] point: center=(2.000, 0.000), r=0.100, n=87

print(detected)
#> ⚙ Attractor detection results
#> -------------------------------------------------- 
#>   Attractors found: 3
#>   Trajectories analyzed: 300
#>   Clustering method: hclust
#> 
#>   Detected attractors:
#>     [1] point at (-2, 0) (r=0.100)
#>     [2] point at (2, 0) (r=0.100)
#>     [3] point at (2, 0) (r=0.100)
#> --------------------------------------------------
```

``` r

plot(detected)
```

![Automatic attractor detection results showing clustered trajectory
endpoints.](custom-systems_files/figure-html/detect-plot-1.png)

Automatic attractor detection results showing clustered trajectory
endpoints.

You can use the detected attractors directly in a compiled system:

``` r

triple_well <- compiled_system(
    cpp_dynamics = '
        double x = state[0];
        double v = state[1];
        deriv[0] = v;
        deriv[1] = -0.5 * v - state[0] * (state[0]*state[0] - 4);
    ',
    attractors = detected$attractors,
    params = list(),
    name = "Triple-well system",
    verbose = FALSE
)

basins_triple <- compute_basins(
    triple_well,
    x_range = c(-3, 3),
    y_range = c(-3, 3),
    resolution = 200,
    t_max = 100,
    verbose = FALSE
)
```

``` r

plot(basins_triple, title = "Triple-well system (auto-detected attractors)")
```

![Basins of the triple-well system using automatically detected
attractors.](custom-systems_files/figure-html/triple-plot-1.png)

Basins of the triple-well system using automatically detected
attractors.

------------------------------------------------------------------------

  

## Comparative summary

| System | Type | Dimension | Basin entropy | Boundary fraction | Wada (grid) | W-statistic |
|:---|:---|---:|:---|:---|:---|:---|
| Duffing oscillator | ODE | 2 | 0.138 | 14.35% | No | 0.562 |
| Henon map | Map | 2 | 0.010 | 1.10% | No | 0.000 |
| Van der Pol oscillator | ODE | 2 | 0.029 | 3.10% | No | 0.360 |

Comparative analysis of basin properties across the three systems
{.table style="width:100%;"}

------------------------------------------------------------------------

  

## Performance considerations

Since all custom systems compile to C++ with OpenMP, they achieve
performance similar to built-in systems:

| System type           | Implementation   | Relative speed |
|:----------------------|:-----------------|:---------------|
| Built-in pendulum     | C++ (OpenMP)     | 1x (baseline)  |
| Built-in Henon-Heiles | C++ (symplectic) | ~1x            |
| Newton fractals       | C++ (parallel)   | ~0.5x (faster) |
| Custom compiled ODE   | C++ (OpenMP)     | ~1x            |
| Custom compiled map   | C++ (OpenMP)     | ~1x            |

Performance comparison of wadaR implementations {.table}

**Recommendations:**

- Use `resolution = 100-200` for exploration
- Increase to `300-500` for publication-quality figures
- Reduce `t_max` if trajectories converge quickly
- Use larger `dt` for smooth systems (monitor accuracy)

------------------------------------------------------------------------

  

## C++ syntax reference

  

### Working with parameters

Parameters are embedded as `const double` values:

``` cpp
// If params = list(a = 1.5, b = 0.3)
// The C++ code receives:
// const double a = 1.5;
// const double b = 0.3;

deriv[0] = a * state[0] - b * state[1];
```

  

### Available math functions

All C++ `<cmath>` functions are available:

``` cpp
sin(x), cos(x), tan(x)           // Trigonometric
asin(x), acos(x), atan(x)        // Inverse trigonometric
atan2(y, x)                       // Two-argument arctangent
exp(x), log(x), log10(x)         // Exponential and logarithmic
pow(x, y), sqrt(x), cbrt(x)      // Powers and roots
fabs(x), floor(x), ceil(x)       // Absolute value and rounding
sinh(x), cosh(x), tanh(x)        // Hyperbolic functions
```

  

### Conditional logic

Use standard C++ conditionals:

``` cpp
if (state[0] > 0) {
    deriv[0] = -state[0];
} else {
    deriv[0] = state[0] * state[0];
}

// Ternary operator
deriv[1] = (state[1] > 0 ? -1.0 : 1.0) * state[1];
```

  

### Local variables

You can declare local variables for clarity:

``` cpp
double x = state[0];
double v = state[1];
double x2 = x * x;
double x3 = x2 * x;

deriv[0] = v;
deriv[1] = -delta * v - alpha * x - beta * x3 + gamma_f * cos(omega * t);
```

------------------------------------------------------------------------

  

## Bifurcation analysis

Bifurcation analysis reveals how the structure of basins of attraction
changes as system parameters vary. This is crucial for understanding:

- **Basin erosion**: Gradual reduction in basin area [\[15\]](#ref15)
- **Basin metamorphosis**: Qualitative changes in basin shape
- **Crisis events**: Sudden changes in attractor/basin structure
  [\[16\]](#ref16)
- **Wada transitions**: Emergence or loss of the Wada property

  

### Theory

As a parameter \\\mu\\ varies, basins can undergo dramatic changes. The
basin entropy \\S_b(\mu)\\ and Wada statistic \\W(\mu)\\ quantify these
transitions:

\\S_b(\mu) = -\sum\_{i=1}^{N_A} p_i(\mu) \log p_i(\mu)\\

where \\p_i(\mu)\\ is the probability of finding basin \\i\\ at
parameter \\\mu\\.

  

### Duffing oscillator forcing amplitude sweep

We analyze how the Duffing oscillator basins change as the forcing
amplitude \\\gamma\\ varies from 0.15 to 0.45:

``` r

# Parameter sweep for Duffing oscillator
bif_duffing <- bifurcation_basins(
    cpp_dynamics = '
        deriv[0] = state[1];
        deriv[1] = -delta * state[1] - alpha * state[0]
                   - beta * pow(state[0], 3)
                   + gamma_f * cos(omega * t);
    ',
    params = list(delta = 0.3, alpha = -1, beta = 1, omega = 1.2),
    sweep_param = "gamma_f",
    sweep_values = seq(0.15, 0.45, length.out = 12),
    attractors = list(
        attractor_point(c(1, 0), 0.3, "Right well"),
        attractor_point(c(-1, 0), 0.3, "Left well"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-2, 2),
    y_range = c(-2, 2),
    resolution = 100,
    t_max = 80,
    compute_entropy = TRUE,
    compute_wada = FALSE,
    verbose = FALSE
)

print(bif_duffing)
#> Bifurcation analysis of basins of attraction
#> -------------------------------------------------- 
#> Sweep parameter: gamma_f 
#> Parameter range: [0.15, 0.45] 
#> Number of steps: 12 
#> Grid resolution: 100 x 100 
#> 
#> Basin entropy range: [0.0721, 0.5364]
#> Boundary fraction range: [7.78%, 51.88%]
#> 
#> Distinct basins detected: 4, 3
```

  

#### Basin entropy evolution

``` r

plot(bif_duffing, type = "entropy")
```

![Basin entropy as a function of forcing amplitude in the Duffing
oscillator. Higher entropy indicates more complex, interleaved basin
structures.](custom-systems_files/figure-html/bifurcation-entropy-plot-1.png)

Basin entropy as a function of forcing amplitude in the Duffing
oscillator. Higher entropy indicates more complex, interleaved basin
structures.

  

#### Basin snapshots

``` r

plot(bif_duffing, type = "basins", indices = c(1, 4, 8, 12))
```

![Basins of attraction at different forcing amplitudes showing the
transition from regular to complex
structure.](custom-systems_files/figure-html/bifurcation-snapshots-1.png)

Basins of attraction at different forcing amplitudes showing the
transition from regular to complex structure.

  

#### Boundary fraction evolution

``` r

plot(bif_duffing, type = "boundary")
```

![Boundary fraction evolution showing increased mixing at higher
forcing.](custom-systems_files/figure-html/bifurcation-boundary-1.png)

Boundary fraction evolution showing increased mixing at higher forcing.

| Parameter (\\\gamma\\) | Basin entropy | Boundary fraction | Distinct basins |
|:-----------------------|:--------------|:------------------|----------------:|
| 0.150                  | 0.0865        | 9.24%             |               4 |
| 0.177                  | 0.0793        | 8.52%             |               4 |
| 0.205                  | 0.0721        | 7.78%             |               4 |
| 0.232                  | 0.1175        | 12.44%            |               4 |
| 0.259                  | 0.2071        | 22.58%            |               4 |
| 0.286                  | 0.2460        | 26.89%            |               4 |
| 0.314                  | 0.4519        | 45.96%            |               4 |
| 0.341                  | 0.4505        | 45.60%            |               3 |
| 0.368                  | 0.4514        | 46.08%            |               3 |
| 0.395                  | 0.4744        | 46.96%            |               3 |
| 0.423                  | 0.5364        | 51.88%            |               3 |
| 0.450                  | 0.5245        | 50.48%            |               3 |

Bifurcation analysis results for the Duffing oscillator {.table}

  

### Extracting basins at specific parameters

The
[`extract_basin()`](https://robustecologies.github.io/wadaR/reference/extract_basin.md)
function retrieves a full basin result at any parameter value for
further analysis:

``` r

# Extract basin at gamma = 0.35
basin_035 <- extract_basin(bif_duffing, param_value = 0.35)
print(basin_035)
#> Wada basins of attraction
#> ---------------------------------------- 
#> Resolution: 100 x 100 
#> Number of attractors: 3 
#> X range: [-2.000, 2.000] 
#> Y range: [-2.000, 2.000]

# Compute detailed entropy
entropy_035 <- basin_entropy(basin_035, box_size = 3)
cat(sprintf("Basin entropy at gamma=0.35: %.4f\n", entropy_035$S_b))
#> Basin entropy at gamma=0.35: 0.6473
```

------------------------------------------------------------------------

  

## 3D visualization

For higher-dimensional systems, wadaR provides 3D visualization
capabilities using plotly for interactive exploration.

  

### The Lorenz system

The Lorenz system [\[17\]](#ref17) is a classic 3D chaotic system:

\\\frac{dx}{dt} = \sigma(y - x)\\ \\\frac{dy}{dt} = x(\rho - z) - y\\
\\\frac{dz}{dt} = xy - \beta z\\

With standard parameters \\\sigma = 10\\, \\\rho = 28\\, \\\beta =
8/3\\, the system exhibits the famous butterfly-shaped strange attractor
with two wings.

``` r

# Compute 3D basins for Lorenz system
basins_lorenz <- compute_basins_3d(
    cpp_dynamics = '
        deriv[0] = sigma * (state[1] - state[0]);
        deriv[1] = state[0] * (rho - state[2]) - state[1];
        deriv[2] = state[0] * state[1] - beta_l * state[2];
    ',
    params = list(sigma = 10, rho = 28, beta_l = 8.0/3.0),
    dim = 3,
    attractors = list(
        attractor_point(c(sqrt(72), sqrt(72), 27), 5, "Right wing"),
        attractor_point(c(-sqrt(72), -sqrt(72), 27), 5, "Left wing"),
        attractor_exit(0, "Escape")
    ),
    x_range = c(-25, 25),
    y_range = c(-35, 35),
    z_range = c(0, 50),
    resolution = 40,
    t_max = 30,
    verbose = FALSE
)

print(basins_lorenz)
#> 3D basin of attraction result
#> ---------------------------------------- 
#> Grid dimensions: 40 x 40 x 40 
#> Total points: 64,000 
#> X range: [-25, 25] 
#> Y range: [-35, 35] 
#> Z range: [0, 50] 
#> Number of attractors: 3 
#> 
#> Basin distribution:
#>   Basin 0: 35.1%
#>   Basin 1: 32.5%
#>   Basin 2: 32.5%
```

  

#### Interactive 3D visualization

``` r

if (requireNamespace("plotly", quietly = TRUE)) {
    plot_3d_basins(basins_lorenz, type = "scatter",
                   subsample = 0.3, opacity = 0.7)
} else {
    cat("Install 'plotly' for interactive 3D visualization")
}
```

Interactive 3D scatter plot of Lorenz system basins. Points colored by
their basin of attraction.

  

#### 2D slices through 3D basins

We can extract 2D cross-sections at specific z-values to analyze the
basin structure:

``` r

# Extract slices at different z levels
slice_z10 <- slice_3d_basins(basins_lorenz, z_slice = 10)
slice_z27 <- slice_3d_basins(basins_lorenz, z_slice = 27)
slice_z40 <- slice_3d_basins(basins_lorenz, z_slice = 40)

# Create plots
p_z10 <- plot(slice_z10, show_boundary = FALSE) +
    ggplot2::ggtitle("z = 10")
p_z27 <- plot(slice_z27, show_boundary = FALSE) +
    ggplot2::ggtitle("z = 27 (attractor height)")
p_z40 <- plot(slice_z40, show_boundary = FALSE) +
    ggplot2::ggtitle("z = 40")

if (requireNamespace("patchwork", quietly = TRUE)) {
    p_z10 + p_z27 + p_z40
} else {
    print(p_z27)
}
```

![2D slices through the Lorenz system basins at different z-values,
showing the basin structure at various
heights.](custom-systems_files/figure-html/lorenz-slices-1.png)

2D slices through the Lorenz system basins at different z-values,
showing the basin structure at various heights.

  

#### Basin entropy of slices

``` r

# Compute entropy for different slices
entropy_z10 <- basin_entropy(slice_z10, box_size = 2)
entropy_z27 <- basin_entropy(slice_z27, box_size = 2)
entropy_z40 <- basin_entropy(slice_z40, box_size = 2)

cat(sprintf("Entropy at z=10: %.4f\n", entropy_z10$S_b))
#> Entropy at z=10: 0.1116
cat(sprintf("Entropy at z=27: %.4f\n", entropy_z27$S_b))
#> Entropy at z=27: 0.1845
cat(sprintf("Entropy at z=40: %.4f\n", entropy_z40$S_b))
#> Entropy at z=40: 0.2079
```

| z-level | Basin entropy | Boundary boxes | Interpretation |
|---:|:---|---:|:---|
| 10 | 0.1116 | 44 | Below attractor, complex mixing |
| 27 | 0.1845 | 70 | At attractor height, maximum complexity |
| 40 | 0.2079 | 80 | Above attractor, simpler structure |

Basin entropy at different z-slices through the Lorenz system {.table}

------------------------------------------------------------------------

  

## Summary

The
[`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
architecture provides a high-performance framework for analyzing Wada
basins in arbitrary dynamical systems:

1.  **High performance**: C++/OpenMP compilation matches built-in system
    speed
2.  **Unified interface**: All wadaR methods work seamlessly with
    compiled systems
3.  **Two system types**: ODEs and discrete maps
4.  **Flexible attractors**: Points, cycles, exits, and discrete
    outcomes
5.  **Automatic detection**:
    [`detect_attractors()`](https://robustecologies.github.io/wadaR/reference/detect_attractors.md)
    finds attractors automatically
6.  **Esc interrupt**: Press Esc to cancel long computations
7.  **Bifurcation analysis**: Parameter sweeps reveal basin transitions
8.  **3D visualization**: Interactive plotly visualizations for
    higher-dimensional systems

For advanced usage, the C++ dynamics code format allows full use of C++
syntax including conditionals, local variables, and all standard math
functions.

------------------------------------------------------------------------

  

## References

**\[1\]** Daza, A., Wagemakers, A., Georgeot, B., Guery-Odelin, D., &
Sanjuan, M. A. F. (2016). Basin entropy: A new tool to analyze
uncertainty in dynamical systems. *Scientific Reports*, 6, 31416. DOI:
[10.1038/srep31416](https://doi.org/10.1038/srep31416)

**\[2\]** Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica
D: Nonlinear Phenomena*, 51(1-3), 213-225. DOI:
[10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789(91)90234-Z)

**\[3\]** Nusse, H. E., & Yorke, J. A. (1996). Wada basin boundaries and
basin cells. *Physica D: Nonlinear Phenomena*, 90(3), 242-261. DOI:
[10.1016/0167-2789(95)00249-9](https://doi.org/10.1016/0167-2789(95)00249-9)

**\[4\]** Sweet, D., Ott, E., & Yorke, J. A. (1999). Topology in chaotic
scattering. *Nature*, 399(6734), 315-316. DOI:
[10.1038/20573](https://doi.org/10.1038/20573)

**\[5\]** Guckenheimer, J., & Holmes, P. (1983). *Nonlinear
oscillations, dynamical systems, and bifurcations of vector fields*.
Springer-Verlag. ISBN: 978-0-387-90819-9. DOI:
[10.1007/978-1-4612-1140-2](https://doi.org/10.1007/978-1-4612-1140-2)

**\[6\]** Henon, M. (1976). A two-dimensional mapping with a strange
attractor. *Communications in Mathematical Physics*, 50(1), 69-77. DOI:
[10.1007/BF01608556](https://doi.org/10.1007/BF01608556)

**\[7\]** Van der Pol, B. (1926). On relaxation-oscillations. *The
London, Edinburgh, and Dublin Philosophical Magazine and Journal of
Science*, 2(11), 978-992. DOI:
[10.1080/14786442608564127](https://doi.org/10.1080/14786442608564127)

**\[8\]** Duffing, G. (1918). *Erzwungene Schwingungen bei
veranderlicher Eigenfrequenz und ihre technische Bedeutung*. Vieweg,
Braunschweig. (Forced oscillations with variable natural frequency and
their technical significance)

**\[9\]** Ueda, Y. (1980). Steady motions exhibited by Duffing’s
equation: A picture book of regular and chaotic motions. In: Holmes,
P.J. (eds) *New Approaches to Nonlinear Problems in Dynamics*. SIAM,
Philadelphia, 311-322.

**\[10\]** Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018).
Ascertaining when a basin is Wada: the merging method. *Scientific
Reports*, 8, 9954. DOI:
[10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

**\[11\]** Grassberger, P., & Procaccia, I. (1983). Measuring the
strangeness of strange attractors. *Physica D: Nonlinear Phenomena*,
9(1-2), 189-208. DOI:
[10.1016/0167-2789(83)90298-1](https://doi.org/10.1016/0167-2789(83)90298-1)

**\[12\]** Strogatz, S. H. (2015). *Nonlinear dynamics and chaos: With
applications to physics, biology, chemistry, and engineering* (2nd ed.).
Westview Press. ISBN: 978-0-8133-4910-7

**\[13\]** Cartwright, M. L., & Littlewood, J. E. (1945). On non-linear
differential equations of the second order: I. *Journal of the London
Mathematical Society*, s1-20(3), 180-189. DOI:
[10.1112/jlms/s1-20.3.180](https://doi.org/10.1112/jlms/s1-20.3.180)

**\[14\]** Parlitz, U., & Lauterborn, W. (1987). Period-doubling
cascades and devil’s staircases of the driven Van der Pol oscillator.
*Physical Review A*, 36(3), 1428-1434. DOI:
[10.1103/PhysRevA.36.1428](https://doi.org/10.1103/PhysRevA.36.1428)

**\[15\]** Aguirre, J., Viana, R. L., & Sanjuan, M. A. F. (2009).
Fractal structures in nonlinear dynamics. *Reviews of Modern Physics*,
81(1), 333-386. DOI:
[10.1103/RevModPhys.81.333](https://doi.org/10.1103/RevModPhys.81.333)

**\[16\]** Grebogi, C., Ott, E., & Yorke, J. A. (1983). Crises, sudden
changes in chaotic attractors, and transient chaos. *Physica D:
Nonlinear Phenomena*, 7(1-3), 181-200. DOI:
[10.1016/0167-2789(83)90126-4](https://doi.org/10.1016/0167-2789(83)90126-4)

**\[17\]** Lorenz, E. N. (1963). Deterministic nonperiodic flow.
*Journal of the Atmospheric Sciences*, 20(2), 130-141. DOI:
[10.1175/1520-0469(1963)020\<0130:DNF\>2.0.CO;2](https://doi.org/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2)
