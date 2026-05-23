# Stationary density of the 2D Fokker-Planck equation

For a two-dimensional drift field \\(f_1(x,y), f_2(x,y))\\ and noise
intensity \\\varepsilon\\, the stationary Fokker-Planck equation is
\$\$-\partial_x(f_1 p) - \partial_y(f_2 p) +
\frac{\varepsilon}{2}\partial\_{xx}(g_1^2 p) +
\frac{\varepsilon}{2}\partial\_{yy}(g_2^2 p) = 0\$\$ where \\g_1(x,y)\\
and \\g_2(x,y)\\ are the diffusion coefficients. This function computes
the stationary density numerically using the Chang-Cooper scheme on a 2D
tensor-product grid, solving the linear system \\Lp = 0\\ with the
constraint \\\sum p = 1\\ via power iteration.

## Usage

``` r
fp_stationary_2d(
  model,
  xlim,
  ylim,
  n_grid = c(100L, 100L),
  epsilon = 0.1,
  max_iter = 10000L,
  tol = 1e-08,
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec` with exactly two state variables. The model may be an
  ODE or SDE; only the drift and (if present) diffusion components are
  used.

- xlim:

  Numeric vector of length 2 giving the domain for the first state
  variable.

- ylim:

  Numeric vector of length 2 giving the domain for the second state
  variable.

- n_grid:

  Integer vector of length 1 or 2 giving the number of grid points in
  each dimension (default `c(100L, 100L)`).

- epsilon:

  Positive numeric noise intensity (default 0.1).

- max_iter:

  Maximum number of power iterations (default 10000).

- tol:

  Convergence tolerance on the relative change in the density vector
  (default 1e-8).

- normalize:

  Logical; if `TRUE` (default), normalize the density so that it
  integrates to 1 over the domain.

- verbose:

  Logical; if `TRUE` (default), print progress messages.

## Value

An S3 object of class `fp_stationary_2d` containing:

- x_grid:

  Numeric vector of grid points for the first state.

- y_grid:

  Numeric vector of grid points for the second state.

- density:

  Matrix (Nx x Ny) of stationary density values.

- state_names:

  Character vector of length 2 with state variable names.

- epsilon:

  The noise intensity used.

- xlim:

  Domain limits for the first state.

- ylim:

  Domain limits for the second state.

- converged:

  Logical indicating whether power iteration converged.

- n_iter:

  Number of iterations performed.

- model:

  The input `system_spec`.

- parms:

  Parameter values used.

## Details

The discretization extends Chang & Cooper (1970) to two dimensions via
operator splitting. The full operator \\L = L_x + L_y\\ is decomposed
into x-sweeps and y-sweeps, each using the 1D Chang-Cooper upwind-biased
weighting at half-grid points. The stationary density is obtained by
finding the null vector of \\L\\ through inverse power iteration:
repeatedly solving \\(L - \sigma I)p^{(k+1)} = p^{(k)}\\ with a small
shift \\\sigma\\ to regularize the singular system, then normalizing.

For SDE models, the diffusion coefficients \\g_1(x,y)\\ and \\g_2(x,y)\\
are read from the model's diffusion formulas. For ODE models, unit
diffusion (additive noise) is assumed.

## References

Chang, J. S. & Cooper, G. (1970). A practical difference scheme for
Fokker-Planck equations. *Journal of Computational Physics*, 6(1), 1-16.
DOI:
[doi:10.1016/0021-9991(70)90001-X](https://doi.org/10.1016/0021-9991%2870%2990001-X)
.

Risken, H. (1996). *The Fokker-Planck equation: methods of solution and
applications* (2nd ed.). Springer. ISBN: 978-3-540-61530-9.

## See also

[`fp_stationary`](https://robustecologies.github.io/janos/reference/fp_stationary.md),
[`fp_kramers_rate`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md),
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 2D Ornstein-Uhlenbeck process
ou2d <- system_spec(
    rhs = list(x ~ -a * x, y ~ -b * y),
    diffusion = list(x ~ sigma, y ~ sigma),
    state_names = c("x", "y"),
    parms = list(a = 1, b = 2, sigma = 1),
    init = c(x = 0, y = 0)
)

fp2d <- fp_stationary_2d(ou2d, xlim = c(-3, 3), ylim = c(-3, 3),
                           n_grid = c(80, 80), epsilon = 0.5)

# Inspect results
print(fp2d)
summary(fp2d)

# Visualize
plot(fp2d)
plot(fp2d, type = "surface")
} # }
```
