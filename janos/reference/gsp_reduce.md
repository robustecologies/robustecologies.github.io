# Geometric singular perturbation reduction

Reduce a slow-fast ODE system to its slow-manifold dynamics, either at
zeroth order (Tikhonov) or with the first-order Fenichel correction.
Returns a closure for the reduced right-hand side over the slow
variables.

## Usage

``` r
gsp_reduce(model, partition, order = 0)
```

## Arguments

- model:

  A janos::system_spec or list with `rhs`, `state`, `parms`.

- partition:

  A `slow_fast_partition` object identifying fast and slow blocks.

- order:

  Reduction order: 0 (Tikhonov) or 1 (Fenichel) (default 0).

## Value

A list of class `gsp_reduction` with fields:

- rhs:

  Closure `function(t, y_slow, parms)` returning the time derivative of
  the slow variables on the slaved manifold.

- slow_names:

  Character vector of slow-variable names.

- fast_names:

  Character vector of fast-variable names (slaved).

- order:

  Integer reduction order in `{0L, 1L}`.

- parent_model:

  The model passed in, retained for downstream use.

## Details

At order 0 (Tikhonov 1952) the fast variables are slaved to the slow
through the algebraic relation f_fast(t, y_fast, y_slow) = 0; the
Newton-solved slaved values are substituted back into f_slow. At order 1
the Fenichel 1979 correction adds a term proportional to the slow-time
derivative of the slaved fast state.

The algebraic solve at order 0 runs a finite-difference Newton iteration
on the fast residual, warm-started from the current fast state. Higher-
order corrections and global manifold continuation are not provided; a
user-supplied closure may substitute the slaving map if a richer
construction is needed.

## References

Fenichel, N. (1979). Geometric singular perturbation theory for ordinary
differential equations. *J. Differential Equations* 31(1), 53-98.
[doi:10.1016/0022-0396(79)90152-9](https://doi.org/10.1016/0022-0396%2879%2990152-9)
.

Tikhonov, A. N. (1952). Systems of differential equations containing
small parameters in the derivatives. *Matematicheskii Sbornik* 31(73),
575-586.

## See also

[`slow_fast_partition`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md),
[`stiffness_ratio`](https://robustecologies.github.io/janos/reference/stiffness_ratio.md).

## Examples

``` r
if (FALSE) { # \dontrun{
rhs <- function(t, y, p) c(
    -100 * (y[1] - sin(y[3])),
    -100 * (y[2] - cos(y[3])),
    p$eps * (1 - y[3])
)
model <- list(rhs = rhs,
              state = c(x1 = 0, x2 = 1, s = 0),
              parms = list(eps = 1e-3))
part <- slow_fast_partition(model, method = "eigen")
red  <- gsp_reduce(model, part, order = 0L)
print(red)
} # }
```
