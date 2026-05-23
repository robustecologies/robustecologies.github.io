# Bifurcation diagram

For each value of a chosen parameter the model is simulated to the
attractor, the transient discarded, and a thin slice of the attractor
(all iterates for a map, or local maxima of the observable for a flow)
is recorded. The resulting scatter of (parameter, observable) is the
classical bifurcation diagram that reveals period-doubling cascades,
periodic windows, band mergings and the chaotic sea. Optionally the
leading Lyapunov exponent is computed at each parameter value and
displayed on a second panel beneath the attractor scatter.

Unlike
[`bifurcation_sweep`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
which traces equilibrium branches by pseudo-arclength continuation, this
routine samples the genuine attractor; it complements that tool in
regimes where the long-time dynamics is not an equilibrium.

## Usage

``` r
bifurcation_diagram(
  model,
  par_name,
  par_range,
  observable = NULL,
  n_par = 300L,
  t_max = NULL,
  discard_transient = NULL,
  n_keep = 150L,
  solver = NULL,
  compute_lyapunov = TRUE,
  lyap_kwargs = list(),
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec`.

- par_name:

  Name of the parameter to scan.

- par_range:

  Numeric `c(lo, hi)`.

- observable:

  Name of the state variable to plot. Default: first state.

- n_par:

  Number of parameter values (default 300).

- t_max:

  Simulation horizon at each parameter value (default 400 for flows,
  1500 iterations for maps).

- discard_transient:

  Transient to discard at each parameter value. Default: 60 percent of
  `t_max`.

- n_keep:

  Maximum iterates retained per parameter value after the transient
  (default 150).

- solver:

  Solver spec for flows. Default
  [`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md).

- compute_lyapunov:

  Logical; if TRUE (default) the leading Lyapunov exponent is evaluated
  at every parameter value and stored for rendering on a second panel.

- lyap_kwargs:

  Optional list of extra arguments forwarded to
  [`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
  when `compute_lyapunov = TRUE`. When not supplied, family-aware
  defaults are chosen so that the exponent converges over the scan
  horizon: `t_max = 5000, renorm_interval = 1` for maps, and
  `t_max = 120, dt = 0.02` for flows.

- parallel:

  Logical; dispatch the independent work units across workers. On Unix
  this uses
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (fork); on Windows a PSOCK cluster is created transparently. Default:
  TRUE.

- n_cores:

  Integer or NULL; number of workers. When NULL, defaults to
  `max(1, parallel::detectCores() - 1)`. Silently clamped to 1 inside
  `R CMD check` (`_R_CHECK_LIMIT_CORES_`). The work is chunked so that
  pressing Esc returns a partial result rather than discarding the run;
  the returned S3 object carries `$interrupted` and
  `$metadata$n_completed` when this happens.

- verbose:

  Progress bar (serial mode only) and a one-line announcement of the
  parallel configuration.

## Value

An S3 object of class `bifurcation_diagram` with components `data` (data
frame of `par` and `y`), `lyapunov` (data frame of `par` and `lambda1`,
or NULL), `par_name`, `observable`, `family`, `par_range`, `n_par`,
`t_max`, `discard_transient`, `n_cores`.

## Details

Bifurcation diagram: attractor-level bifurcation over a parameter

The scan is embarrassingly parallel. By default `n_cores` is set to
`max(1, parallel::detectCores() - 1)`; on Unix the workers are forked
through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html), on
Windows they run in a transparent
[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html).
Work is chunked so that pressing Esc returns a partial result and marks
`$interrupted = TRUE` rather than discarding the scan. A progress bar is
shown only in serial mode to avoid interleaved worker output.

When `compute_lyapunov = TRUE` the leading exponent \\\lambda_1\\ is
evaluated via
[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
at each parameter value using the same family dispatch (ODE or discrete
map). Additional arguments for the exponent computation may be passed
through the `lyap_kwargs` list, for instance
`list(t_max = 80, dt = 0.01, renorm_interval = 1)`. The result is stored
in the `lyapunov` slot and automatically rendered by the plot method as
a lower panel stacked via patchwork.

## References

Feigenbaum, M. J. (1978). Quantitative universality for a class of
nonlinear transformations. *Journal of Statistical Physics*, 19(1),
25-52. [doi:10.1007/BF01020332](https://doi.org/10.1007/BF01020332)

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos*, 2nd ed. Westview
Press. ISBN 978-0813349107.

## See also

[`bifurcation_sweep`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
[`continuation`](https://robustecologies.github.io/janos/reference/continuation.md),
[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md),
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
logistic <- system_spec(
    map = list(x ~ r * x * (1 - x)),
    state_names = "x",
    parms = list(r = 3.5),
    init = c(x = 0.2)
)
bd <- bifurcation_diagram(logistic, par_name = "r",
                          par_range = c(2.8, 4.0), n_par = 500,
                          compute_lyapunov = TRUE)
plot(bd)
} # }
```
