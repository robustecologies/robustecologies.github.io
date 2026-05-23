# Estimate largest Lyapunov exponent

Computes the largest Lyapunov exponent using the Wolf or Rosenstein
algorithm, applicable to any turbulence simulation.

## Usage

``` r
turbulence_lyapunov(
  x,
  dt = NULL,
  embed_dim = 10L,
  tau = 10L,
  method = c("wolf", "rosenstein"),
  epsilon = 0.1,
  evolve_steps = 10L,
  max_scale = 1,
  theiler_window = 50L,
  max_steps = 500L
)
```

## Arguments

- x:

  Time series (numeric vector) or simulation object.

- dt:

  Time step of the simulation. If `x` is a simulation object, `dt` is
  extracted automatically. The raw Lyapunov exponent (per step) is
  converted to continuous time units via \\\lambda = \lambda\_{step} /
  dt\\. If `NULL` and cannot be inferred, defaults to 1 (no conversion).

- embed_dim:

  Embedding dimension for reconstruction.

- tau:

  Time delay for embedding.

- method:

  Algorithm: "wolf" (default) or "rosenstein".

- epsilon:

  Initial neighborhood radius (Wolf method).

- evolve_steps:

  Steps between neighbor replacements (Wolf).

- max_scale:

  Maximum scale before replacement (Wolf).

- theiler_window:

  Temporal separation for neighbors (Rosenstein).

- max_steps:

  Maximum divergence steps to track (Rosenstein).

## Value

Object of class "turbulence_lyapunov" containing the exponent estimate
and diagnostic information.

## Details

The largest Lyapunov exponent \\\lambda_1\\ quantifies the average rate
of exponential divergence (or convergence) of nearby trajectories in
phase space: \$\$\|\delta(t)\| \sim \|\delta(0)\| e^{\lambda_1 t}\$\$
where \\\delta(t)\\ represents the separation between initially close
trajectories at time \\t\\. A positive exponent (\\\lambda_1 \> 0\\)
indicates deterministic chaos with sensitive dependence on initial
conditions; trajectories diverge exponentially, making long-term
prediction impossible despite the underlying determinism. Negative
exponents indicate asymptotically stable behavior where perturbations
decay, while zero exponents characterize marginally stable systems or
limit cycles.

**Wolf algorithm (1985):** The original method for Lyapunov estimation
from experimental data. The algorithm tracks a single fiducial
trajectory and its nearest neighbor in the reconstructed phase space.
When the separation exceeds a threshold (`max_scale`), the current
neighbor is replaced with a new one lying within `epsilon` of the
fiducial point, maintaining similar orientation. The exponent is
computed as: \$\$\lambda_1 = \frac{1}{t_M} \sum\_{k=1}^{M} \log_2
\frac{L'(t_k)}{L(t\_{k-1})}\$\$ where \\L(t_k)\\ and \\L'(t_k)\\ are
distances before and after replacement. The Wolf method is robust but
computationally intensive and requires careful parameter tuning. The
`evolve_steps` parameter controls how many time steps pass between
potential neighbor replacements.

**Rosenstein algorithm (1993):** A simpler, faster method particularly
suitable for small datasets. For each point on the trajectory, the
algorithm finds the nearest neighbor (with temporal separation exceeding
`theiler_window` to avoid correlated neighbors) and tracks the
divergence over time. The mean logarithmic divergence is computed across
all reference points: \$\$S(k) = \frac{1}{N} \sum\_{j=1}^{N} \ln
d_j(k)\$\$ where \\d_j(k)\\ is the distance between the \\j\\-th pair of
neighbors after \\k\\ time steps. The Lyapunov exponent is estimated
from the slope of the linear region of \\S(k)\\ versus \\k\\. The
`max_steps` parameter controls the maximum divergence time tracked. This
method produces a divergence curve that can be visualized with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

**Method selection:** Use `"rosenstein"` for exploratory analysis and
visualization of divergence behavior; it is faster and produces
interpretable diagnostics. Use `"wolf"` for comparison with historical
literature values, as it remains the standard reference algorithm. Both
methods require adequate embedding dimension (`embed_dim`) determined
via
[`turbulence_cao`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md)
and appropriate time delay (`tau`).

**Reference Lyapunov exponents** for systems at standard parameters:

- Lorenz (\\\sigma\\=10, \\\rho\\=28, \\\beta\\=8/3): \\\lambda_1
  \approx 0.906\\

- Rössler (a=0.2, b=0.2, c=5.7): \\\lambda_1 \approx 0.07\\

- Lorenz-84 (a=0.25, b=4, F=8, G=1): \\\lambda_1 \approx 0.07\\

- vonkarman (\\\mu\\=0.3): \\\lambda_1 \approx 0.02-0.1\\

- charney_devore: typically positive when multistable

- stommel/hasselmann: may be negative (stable) or near zero

## References

Wolf, A., Swift, J.B., Swinney, H.L. & Vastano, J.A. (1985). Determining
Lyapunov exponents from a time series. *Physica D*, 16(3), 285-317.
[doi:10.1016/0167-2789(85)90011-9](https://doi.org/10.1016/0167-2789%2885%2990011-9)

Rosenstein, M.T., Collins, J.J. & De Luca, C.J. (1993). A practical
method for calculating largest Lyapunov exponents from small data sets.
*Physica D*, 65(1-2), 117-134.
[doi:10.1016/0167-2789(93)90009-P](https://doi.org/10.1016/0167-2789%2893%2990009-P)

## See also

[`print.turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_lyapunov.md),
[`summary.turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_lyapunov.md),
[`plot.turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_lyapunov.md),
[`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md),
[`turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md),
[`turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md),
[`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Lorenz system (reference lambda_1 = 0.906)
sim_lorenz <- lorenz_sim(sigma = 10, rho = 28, beta = 8/3,
                         n_steps = 200000, seed = 42)

# Wolf algorithm
lyap_wolf <- turbulence_lyapunov(sim_lorenz, embed_dim = 10, tau = 10,
                                  method = "wolf")
print(lyap_wolf)

# Rosenstein algorithm with divergence curve
lyap_ros <- turbulence_lyapunov(sim_lorenz, embed_dim = 10, tau = 10,
                                 method = "rosenstein", max_steps = 300)
print(lyap_ros)
plot(lyap_ros)

# Von Karman stochastic attractor
sim_vk <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)
lyap_vk <- turbulence_lyapunov(sim_vk, embed_dim = 10, tau = 10)
print(lyap_vk)

# Rossler system (reference lambda_1 = 0.07)
sim_rossler <- rossler_sim(a = 0.2, b = 0.2, c = 5.7,
                           n_steps = 500000, seed = 42)
lyap_rossler <- turbulence_lyapunov(sim_rossler, embed_dim = 10, tau = 15,
                                     method = "rosenstein")
print(lyap_rossler)
plot(lyap_rossler)
} # }
```
