# Compute the stationary density of the Fokker-Planck equation

For a one-dimensional drift field \\f(x)\\ and noise intensity
\\\varepsilon\\, the Fokker-Planck equation \$\$\partial_t u =
\varepsilon\\\partial\_{xx} u - \partial_x(f(x)\\u)\$\$ has a unique
positive stationary solution \\u^{f,\varepsilon}(x)\\ under reflecting
boundary conditions. This function computes that stationary density
numerically by finding the null vector of the discretized Fokker-Planck
operator, constructed using the Chang-Cooper scheme for non-negative,
probability-conserving discretizations.

## Usage

``` r
fp_stationary(
  model,
  xlim,
  n_grid = 200L,
  epsilon = 0.1,
  gradient = TRUE,
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with exactly one state variable (1D). The model may be
  an ODE (drift only) or SDE (drift + diffusion); only the drift
  component is used.

- xlim:

  Numeric vector of length 2 giving the spatial domain
  `c(x_min, x_max)`.

- n_grid:

  Integer number of grid points (default 200).

- epsilon:

  Positive numeric noise intensity (default 0.1).

- gradient:

  Logical; if `TRUE` (default), also compute the drift potential and
  Boltzmann distribution for comparison.

- normalize:

  Logical; if `TRUE` (default), normalize the density so that it
  integrates to 1 over the domain.

- verbose:

  Logical; if `TRUE` (default), print progress messages.

## Value

An S3 object of class `fp_stationary` containing:

- x:

  Numeric vector of grid points.

- density:

  Numeric vector of stationary density values.

- potential:

  Numeric vector of drift potential \\V(x)\\ (if `gradient = TRUE`), or
  `NULL`.

- boltzmann:

  Numeric vector of Boltzmann density (if `gradient = TRUE`), or `NULL`.

- l2_error:

  Relative \\L^2\\ error between numerical and Boltzmann densities (if
  `gradient = TRUE`), or `NULL`.

- epsilon:

  The noise intensity used.

- spectral_gap:

  Magnitude of the second eigenvalue \\\|\lambda_2\|\\ of the FP
  operator.

- eigenvalues:

  The 6 eigenvalues closest to zero.

- xlim:

  Domain limits used.

- n_grid:

  Number of grid points.

- model:

  The input `model_spec`.

- parms:

  Parameter values used.

## Details

The discretization follows Chang & Cooper (1970). At each
half-grid-point \\x\_{i+1/2}\\, the probability flux is written as a
weighted average of neighboring densities with weight \\\delta = 1/w -
1/(e^w - 1)\\ where \\w = f \Delta x / \varepsilon\\ is the local Peclet
number. This upwind-biased scheme guarantees positive densities and
exact probability conservation, unlike central-difference
discretizations which produce spurious oscillations when the Peclet
number exceeds 2.

The stationary density is computed by forward recursion from the
zero-flux condition \\J\_{i+1/2} = 0\\, which is numerically stable and
exact for the discretized system. The spectral gap \\\lambda_2\\,
defined as the magnitude of the second-largest eigenvalue of the
operator, is obtained separately via full eigendecomposition.

In one dimension every drift field is a gradient system: \\f(x) =
-V'(x)\\ where \\V(x) = -\int\_{x_0}^x f(s)\\ds\\. When
`gradient = TRUE`, the function also computes the exact Boltzmann
distribution \\u_B(x) \propto \exp(-2V(x)/\varepsilon)\\ and its \\L^2\\
relative error against the numerically computed density, providing a
built-in accuracy diagnostic.

## References

Chang, J. S. & Cooper, G. (1970). A practical difference scheme for
Fokker-Planck equations. *Journal of Computational Physics*, 6(1), 1-16.
DOI:
[doi:10.1016/0021-9991(70)90001-X](https://doi.org/10.1016/0021-9991%2870%2990001-X)
.

Risken, H. (1996). *The Fokker-Planck equation: methods of solution and
applications* (2nd ed.). Springer. ISBN: 978-3-540-61530-9.

Zeeman, E. C. (1988). Stability of dynamical systems. *Nonlinearity*,
1(1), 115-155. DOI:
[doi:10.1088/0951-7715/1/1/005](https://doi.org/10.1088/0951-7715/1/1/005)
.

## See also

[`fp_potential`](https://robustecologies.github.io/janos/reference/fp_potential.md),
[`fp_kramers_rate`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Ornstein-Uhlenbeck process: exact Gaussian stationary density
ou <- model_spec(
    rhs = list(x ~ -a * x),
    diffusion = list(x ~ sigma),
    state_names = "x",
    parms = list(a = 1.0, sigma = 1.0),
    init = c(x = 0)
)

fp <- fp_stationary(ou, xlim = c(-4, 4), epsilon = 0.5)

# Inspect results
print(fp)
summary(fp)

# Visualize
plot(fp)
plot(fp, type = "potential")
plot(fp, type = "both")
} # }
```
