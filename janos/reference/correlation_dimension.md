# Correlation dimension

Estimates the Grassberger-Procaccia correlation dimension \\D_2\\ of the
strange set visited by a `dyn_sim` trajectory. The correlation integral
\\C(\varepsilon)\\ counts pairs of attractor points within distance
\\\varepsilon\\; on a fractal set it scales as \\C(\varepsilon) \sim
\varepsilon^{D_2}\\ over a scaling region bounded below by noise and
above by the attractor diameter. The slope of \\\log C\\ versus \\\log
\varepsilon\\ in that window is \\D_2\\.

## Usage

``` r
correlation_dimension(
  x,
  vars = NULL,
  n_points = 2000L,
  n_eps = 40L,
  eps_range = NULL,
  theiler = 10L,
  fit_range = c(0.2, 0.7),
  seed = 42
)
```

## Arguments

- x:

  A `dyn_sim` object.

- vars:

  Optional character vector of state variables to embed. By default all
  state variables are used.

- n_points:

  Number of attractor points to sample (default: 2000).

- n_eps:

  Number of epsilon values to scan (default: 40).

- eps_range:

  Optional `c(eps_min, eps_max)` fixing the scan window. Defaults to
  `c(q_0.01, q_0.5)` of the pairwise distance distribution.

- theiler:

  Theiler window in number of samples (default: 10).

- fit_range:

  Fraction of the log-epsilon range to use for the linear fit, as
  `c(lo, hi)` in \\\[0, 1\]\\. Default: `c(0.2, 0.7)`.

- seed:

  Random seed for sampling (default: 42).

## Value

An S3 object of class `correlation_dimension` with components `D2`
(scalar estimate), `eps`, `C` (numeric vectors of the correlation
integral), `fit` (the linear fit), `n_points`, `theiler`.

## Details

Grassberger-Procaccia correlation dimension of an attractor

The implementation uses the Theiler window correction to skip pairs of
samples closer in time than `theiler` steps, which removes spurious
small distances due to temporal contiguity. Distances are computed on
the full state; pass `vars` to restrict the embedding dimension.

## References

Grassberger, P., & Procaccia, I. (1983). Characterization of strange
attractors. *Physical Review Letters*, 50(5), 346-349.
[doi:10.1103/PhysRevLett.50.346](https://doi.org/10.1103/PhysRevLett.50.346)

Theiler, J. (1986). Spurious dimension from correlation algorithms
applied to limited time-series data. *Physical Review A*, 34(3),
2427-2432.
[doi:10.1103/PhysRevA.34.2427](https://doi.org/10.1103/PhysRevA.34.2427)

## See also

[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md),
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 100, discard_transient = 10,
               solver = solver_rk4(dt = 0.01))
cd <- correlation_dimension(run)
print(cd); plot(cd)
} # }
```
