# Drift potential and quasi-potential reconstruction

For a one-dimensional system with drift \\f(x)\\, reconstructs the
deterministic potential \\V(x) = -\int\_{x_0}^x f(s)\\ds\\ via
trapezoidal quadrature. When `epsilon` is provided, also computes the
quasi-potential \\\Phi(x) = -(\varepsilon/2)\ln u(x)\\ from the
numerically computed Fokker-Planck stationary density.

## Usage

``` r
fp_potential(model, xlim, n_grid = 200L, epsilon = NULL, verbose = TRUE)
```

## Arguments

- model:

  A `system_spec` with exactly one state variable.

- xlim:

  Numeric vector of length 2 giving the spatial domain.

- n_grid:

  Integer number of grid points (default 200).

- epsilon:

  Optional positive numeric noise intensity. If provided, the
  quasi-potential is computed from the FP stationary density.

- verbose:

  Logical; print progress messages (default TRUE).

## Value

An S3 object of class `fp_potential` containing:

- x:

  Numeric vector of grid points.

- V:

  Numeric vector of drift potential values.

- Phi:

  Numeric vector of quasi-potential values (if `epsilon` provided), or
  `NULL`.

- wells:

  Data frame of potential wells (local minima) with columns x, V, and
  idx.

- barriers:

  Data frame of potential barriers (local maxima) with columns x, V, and
  idx.

- barrier_heights:

  Numeric vector of barrier heights above each well (one per
  well-barrier pair), or `NULL`.

- epsilon:

  The noise intensity used, or `NULL`.

- xlim:

  Domain limits.

- model:

  The input `system_spec`.

- parms:

  Parameter values used.

## Details

In one dimension, the potential \\V\\ is the negative antiderivative of
the drift, so that \\f = -V'\\. Local minima of \\V\\ correspond to
stable equilibria (wells), local maxima to unstable equilibria
(barriers). The barrier height \\\Delta V\\ between two wells governs
the Kramers escape rate as \\r \sim \exp(-2\Delta V/\varepsilon)\\.

The quasi-potential \\\Phi(x) = -(\varepsilon/2)\ln u(x)\\ is the
Freidlin-Wentzell large-deviation rate function. For gradient systems in
1D, it coincides with \\V(x)\\ up to an additive constant. For systems
with absorbing boundaries or large \\\varepsilon\\, \\\Phi\\ deviates
from \\V\\, capturing the genuine stochastic landscape seen by the
system.

## References

Freidlin, M. I. & Wentzell, A. D. (2012). *Random perturbations of
dynamical systems* (3rd ed.). Springer. ISBN: 978-3-642-25846-6. DOI:
[doi:10.1007/978-3-642-25847-3](https://doi.org/10.1007/978-3-642-25847-3)
.

Nolting, B. C. & Abbott, K. C. (2016). Balls, cups, and
quasi-potentials: quantifying stability in stochastic systems. *Journal
of the Royal Society Interface*, 13(120), 20150772. DOI:
[doi:10.1098/rsif.2015.0772](https://doi.org/10.1098/rsif.2015.0772) .

## See also

[`fp_stationary`](https://robustecologies.github.io/janos/reference/fp_stationary.md),
[`fp_kramers_rate`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Double-well potential from Allee effect
allee <- system_spec(
    rhs = list(x ~ r * x * (x / A - 1) * (1 - x / K)),
    state_names = "x",
    parms = list(r = 1, A = 0.3, K = 1),
    init = c(x = 0.5)
)

pot <- fp_potential(allee, xlim = c(-0.2, 1.3), epsilon = 0.05)

print(pot)
summary(pot)
plot(pot)
plot(pot, type = "both")
} # }
```
