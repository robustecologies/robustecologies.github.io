# Takens time-delay embedding

Reconstruct phase space from a univariate time series using Takens'
delay embedding theorem. This is a fundamental technique for analyzing
dynamical systems from scalar observations.

## Usage

``` r
takens_embed(
  x,
  m = 3,
  tau = 1,
  method_tau = c("acf_half", "acf_1e", "mutual_info"),
  verbose = FALSE
)
```

## Arguments

- x:

  Numeric vector. Univariate time series to embed. Must have at least 10
  observations and no missing or infinite values.

- m:

  Integer. Embedding dimension, the number of delay coordinates. If NULL
  (default), estimated automatically using the false nearest neighbors
  (FNN) method. The minimum dimension required to unfold the attractor
  without self-intersections. Typical values: 2-10.

- tau:

  Integer. Time delay in samples (units of the sampling interval). If
  NULL (default), estimated from the autocorrelation function. Should be
  large enough for decorrelation but small enough to preserve dynamical
  information. Typical values: 1-100 depending on sampling rate.

- method_tau:

  Character. Method for automatic tau estimation:

  "acf_half"

  :   (Default) First lag where ACF drops below 0.5. Balances
      decorrelation and information preservation.

  "acf_1e"

  :   First lag where ACF drops below 1/e. More conservative, may lead
      to better unfolding.

  "mutual_info"

  :   First local minimum of mutual information between x(t) and
      x(t+lag). Information-theoretic approach, computationally more
      expensive but often more reliable.

- verbose:

  Logical. If `TRUE`, print informational messages reporting the values
  chosen by automatic estimation of `tau` and `m` (only printed when the
  corresponding argument is `NULL`). Default `FALSE`.

## Value

A matrix of class `"takens_embedding"` with dimensions (n_points, m),
where n_points = length(x) - (m-1)\*tau. Row i contains the delay vector
\\(x_i, x\_{i+\tau}, \ldots, x\_{i+(m-1)\tau})\\.

The returned object has attributes:

- m:

  Embedding dimension used

- tau:

  Time delay used

- n_original:

  Length of original time series

S3 methods available:
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), which creates
an interactive 3D visualization of the reconstructed attractor.

## Details

The time series is transformed into an m-dimensional state space via the
delay coordinate map:

\$\$\Theta(t) \rightarrow \Theta\_{m,\tau}(t) = \[\Theta(t),
\Theta(t-\tau), \ldots, \Theta(t-(m-1)\tau)\]^\top\$\$

Under suitable conditions (Whitney embedding theorem), this
reconstruction is diffeomorphic to the original attractor, preserving
its geometric and topological properties.

**Theoretical background:**

Takens' theorem (1981) guarantees that for a deterministic dynamical
system, the delay coordinate map provides a faithful reconstruction of
the attractor if the embedding dimension m \> 2\*d, where d is the
attractor's box-counting dimension. In practice, smaller m often
suffices.

**Parameter selection:**

- **Time delay tau**: Too small leads to strongly correlated coordinates
  (the attractor collapses to the diagonal). Too large leads to loss of
  dynamical connection between coordinates. The ACF-based methods find a
  balance automatically.

- **Embedding dimension m**: Too small leads to self- intersecting
  (folded) attractors. Too large adds noise dimensions. The FNN method
  identifies the minimum dimension where false neighbors (caused by
  projection) are eliminated.

**Computational notes:** The core embedding is implemented in C++ for
efficiency. FNN computation uses subsampling for large time series to
maintain reasonable runtime.

## References

Takens F (1981). Detecting strange attractors in turbulence. In: Lecture
Notes in Mathematics, 898. Springer.
[doi:10.1007/BFb0091924](https://doi.org/10.1007/BFb0091924)

Kennel MB, Brown R, Abarbanel HDI (1992). Determining embedding
dimension for phase-space reconstruction using a geometrical
construction. Physical Review A 45:3403.
[doi:10.1103/PhysRevA.45.3403](https://doi.org/10.1103/PhysRevA.45.3403)

Fraser AM, Swinney HL (1986). Independent coordinates for strange
attractors from mutual information. Physical Review A 33:1134.
[doi:10.1103/PhysRevA.33.1134](https://doi.org/10.1103/PhysRevA.33.1134)

Alberti T et al. (2023). Chameleon attractors in turbulent flows. Chaos,
Solitons and Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate Lorenz-like data
set.seed(42)
t <- seq(0, 50, by = 0.01)
x <- sin(t) + 0.5 * sin(2.3 * t) + 0.1 * rnorm(length(t))

# Embed with specified parameters
embedded <- takens_embed(x, m = 3, tau = 20)

# Inspect the embedding
print(embedded)
summary(embedded)

# Visualize the reconstructed attractor
plot(embedded)

# Automatic parameter estimation
embedded_auto <- takens_embed(x, m = NULL, tau = NULL)
print(embedded_auto)
} # }
```
