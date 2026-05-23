# Stiffness ratio of an ODE system at a reference state

Returns lambda_max / lambda_min where lambdas are the absolute values of
the real parts of the Jacobian eigenvalues at a reference state. Also
reports the spectral abscissa (largest real part) and the
field-of-values radius (approximation by Hermitian average).

## Usage

``` r
stiffness_ratio(model, state = NULL, parms = NULL)
```

## Arguments

- model:

  A janos::system_spec or a list with at least `rhs` and `state`.

- state:

  Reference state.

- parms:

  Parameter values.

## Value

A list with fields `ratio`, `spectral_abscissa`, `fov_radius`,
`eigvals`.

## References

Hairer, E. and Wanner, G. (1996). *Solving Ordinary Differential
Equations II*. Springer, Section IV.1.

## See also

[`slow_fast_partition`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md).

## Examples

``` r
if (FALSE) { # \dontrun{
rhs <- function(t, y, p) c(-1e4 * y[1] + y[2], -y[2])
model <- list(rhs = rhs, state = c(y1 = 1, y2 = 1), parms = list())
s <- stiffness_ratio(model)
s$ratio
} # }
```
