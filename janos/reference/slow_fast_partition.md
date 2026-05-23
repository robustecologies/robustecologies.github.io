# Slow-fast partition of an ODE system

Detect or accept a partition of the state vector into a "fast" block and
a "slow" block. Three discovery methods are implemented:

- `eigen`:

  Linearise the right-hand side at a reference state, compute the
  Jacobian eigenvalues, and partition by spectral gap. Components whose
  participation in the largest-magnitude eigenvector exceeds a threshold
  are assigned to the fast block.

- `graph`:

  Inspect the Jacobian sparsity pattern; weakly coupled blocks are
  identified as candidate partitions. Useful when the model has an
  obvious structural separation, for instance a fast biochemistry
  coupled to a slow regulatory layer, or a fast reaction coupled to a
  slow transport process.

- `user`:

  Accept a user-supplied `fast_idx` or `slow_idx` vector of state
  indices.

## Usage

``` r
slow_fast_partition(
  model,
  method = c("eigen", "graph", "user"),
  state = NULL,
  parms = NULL,
  tau_ratio_cutoff = 100,
  fast_idx = NULL,
  slow_idx = NULL
)
```

## Arguments

- model:

  A
  [`janos::system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)
  or any list carrying an `rhs` closure together with a `state` initial
  condition and a `parms` field. For IMEX-RK splits, the model may
  instead carry `rhs_explicit` and `rhs_implicit`; the partition is
  computed against their sum.

- method:

  One of `"eigen"`, `"graph"`, `"user"`.

- state:

  Reference state at which to linearise (default: model\$state).

- parms:

  Parameter values (default: model\$parms).

- tau_ratio_cutoff:

  Minimum spectral-gap ratio required for the partition to be declared
  non-trivial (default 100).

- fast_idx, slow_idx:

  For `method = "user"`, explicit index vectors.

## Value

An object of class `slow_fast_partition` with fields `fast_idx`,
`slow_idx`, `tau_fast`, `tau_slow`, `coupling_norm`, `eigvals`,
`method`.

## References

Fenichel, N. (1979). Geometric singular perturbation theory for ordinary
differential equations. *J. Differential Equations* 31(1), 53-98.
[doi:10.1016/0022-0396(79)90152-9](https://doi.org/10.1016/0022-0396%2879%2990152-9)
. Sastry, S. and Desoer, C. A. (1981). The robustness of controllability
and observability of linear time-varying systems. *IEEE Trans. Automatic
Control* 27(4), 933-939.

## See also

[`stiffness_ratio`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md),
[`gsp_reduce`](https://robustecologies.github.io/janos/reference/gsp_reduce.md),
[`solver_imex_ark`](https://robustecologies.github.io/janos/reference/solver_imex_ark.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Three-state slow-fast system with a clear spectral gap.
rhs <- function(t, y, p) c(
    -100 * (y[1] - sin(y[3])),
    -100 * (y[2] - cos(y[3])),
    p$eps * (1 - y[3])
)
model <- list(rhs = rhs,
              state = c(x1 = 0, x2 = 1, s = 0),
              parms = list(eps = 1e-3))
part <- slow_fast_partition(model, method = "eigen")
print(part)
plot(part)
} # }
```
