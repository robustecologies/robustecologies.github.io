# Compute EVT-based instantaneous dynamical metrics

Compute the instantaneous dimension d(t) and inverse persistence
theta(t) for each point in a phase-space trajectory using extreme value
theory (EVT). These metrics characterize the local geometry and dynamics
of the attractor.

## Usage

``` r
evt_metrics(
  x,
  q = 0.98,
  metric = c("euclidean", "manhattan", "maximum"),
  theta_method = c("sueveges", "ferro", "northrop"),
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- x:

  Numeric matrix. Phase-space trajectory with dimensions (n_points,
  embedding_dimension). Typically the output of
  [`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md).
  Must have at least 100 rows.

- q:

  Numeric. Quantile for threshold selection, between 0.9 and 0.999
  (default 0.98). Higher values focus on more extreme recurrences,
  giving more local estimates but with higher variance.

- metric:

  Character. Distance metric for computing recurrences:

  "euclidean"

  :   (Default) Standard L2 norm. Most commonly used.

  "manhattan"

  :   L1 norm (sum of absolute differences).

  "maximum"

  :   L-infinity norm (Chebyshev distance).

- theta_method:

  Character. Method for extremal index estimation:

  "sueveges"

  :   (Default) Sueveges estimator based on inter- exceedance times.
      Robust for dependent sequences.

  "ferro"

  :   Ferro-Segers estimator. Alternative approach.

  "northrop"

  :   Northrop estimator with bias correction.

- n_cores:

  Integer. Number of CPU cores for parallel computation (default 1). The
  computation is O(n^2) in trajectory length, so parallelization
  provides significant speedup for large datasets.

- verbose:

  Logical. Print progress information (default FALSE).

## Value

A list of class `"evt_metrics"` containing:

- d:

  Numeric vector of length n_points. Instantaneous dimension for each
  state. Values typically range from 1 to the embedding dimension, with
  NA for states where estimation failed.

- theta:

  Numeric vector of length n_points. Inverse persistence (extremal
  index) for each state. Values range from 0 to 1, with NA for states
  where estimation failed.

- mean_d:

  Numeric. Mean instantaneous dimension across all states.

- mean_theta:

  Numeric. Mean inverse persistence across all states.

- sd_d:

  Numeric. Standard deviation of instantaneous dimension.

- sd_theta:

  Numeric. Standard deviation of inverse persistence.

- q:

  Numeric. Quantile used for threshold selection.

- n:

  Integer. Number of trajectory points.

- p:

  Integer. Embedding dimension.

- metric:

  Character. Distance metric used.

- theta_method:

  Character. Extremal index estimation method.

S3 methods available: [`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

The approach is based on the observation that recurrence statistics to a
reference state zeta\* follow extreme value distributions. Specifically,
for distances falling in the tail of the distribution (controlled by
quantile q), the exceedances above a threshold follow a generalized
Pareto distribution:

\$\$F(X, \zeta^\*) \simeq \exp\left\[-\theta(\zeta^\*) \cdot d(\zeta^\*)
\cdot X(\zeta^\*)\right\]\$\$

where d is the local dimension and theta is the extremal index (inverse
persistence).

**Physical interpretation:**

- **Local dimension d(t)**: Measures the effective number of degrees of
  freedom active around state zeta(t). Low d indicates the trajectory is
  confined to a low-dimensional subspace. High d indicates all
  dimensions are dynamically active. Variations in d across the
  attractor reveal heterogeneous geometry.

- **Inverse persistence theta(t)**: Measures the probability that the
  trajectory escapes the neighborhood of zeta(t) after a recurrence.
  theta = 1 indicates memoryless (Poisson-like) recurrences. theta \< 1
  indicates clustering of recurrences, i.e., the trajectory tends to
  stay near the reference point. This reflects local stability.

**Algorithmic details:** For each reference state zeta\*, the algorithm:

1.  Computes distances from zeta\* to all other states

2.  Transforms distances to g = -log(distance)

3.  Identifies the q-quantile threshold

4.  Extracts exceedances above the threshold

5.  Fits an exponential (GPD with xi=0) to estimate d =
    1/mean(exceedances)

6.  Estimates theta from inter-exceedance times using the Sueveges
    method

**Computational complexity:** O(n^2) where n is the trajectory length.
The C++ implementation with OpenMP parallelization helps, but runtime
still grows quickly with n; consider subsampling for very long
trajectories.

## References

Lucarini V, Faranda D, Wouters J (2012). Universal behaviour of extreme
value statistics for selected observables of dynamical systems. J. Stat.
Phys. 147:63-73.
[doi:10.1007/s10955-012-0468-z](https://doi.org/10.1007/s10955-012-0468-z)

Faranda D, Messori G, Yiou P (2017). Dynamical proxies of North Atlantic
predictability and extremes. Scientific Reports 7:41278.
[doi:10.1038/srep41278](https://doi.org/10.1038/srep41278)

Sueveges M (2007). Likelihood estimation of the extremal index. Extremes
10:41-55.
[doi:10.1007/s10687-007-0034-2](https://doi.org/10.1007/s10687-007-0034-2)

Alberti T et al. (2023). Chameleon attractors in turbulent flows. Chaos,
Solitons and Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

## See also

[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for the scale-resolved counterpart computed across MEMD modes;
[`print.evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/print.evt_metrics.md),
[`summary.evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/summary.evt_metrics.md),
[`plot.evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/plot.evt_metrics.md)
for inspection of the returned object;
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for the full embedding plus EVT pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create Lorenz-like embedded data
set.seed(42)
t <- seq(0, 100, by = 0.01)
x <- sin(t) + 0.5*sin(2.3*t) + 0.1*rnorm(length(t))
embedded <- takens_embed(x, m = 3, tau = 20)

# Compute EVT metrics
metrics <- evt_metrics(embedded, q = 0.98, verbose = TRUE)

# Inspect results
print(metrics)
summary(metrics)

# Visualize metrics
plot(metrics)
plot(metrics, type = "histogram")
plot(metrics, type = "scatter")
} # }
```
