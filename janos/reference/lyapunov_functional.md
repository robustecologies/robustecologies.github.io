# Energy Lyapunov functional for a reaction-diffusion PDE

Builds a Lyapunov functional \\V\[u\] = \int\_\Omega \left( \tfrac{D}{2}
\|\nabla u\|^2 + F(u) \right) dx\\ for a 1D reaction-diffusion PDE
specified via
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
with `pde = list(u ~ D * d2x(u) + f(u, ...))`. The functional is a
Lyapunov function along the flow when the reaction \\f\\ derives from a
potential \\F\\ (i.e. \\f = -\partial F / \partial u\\), a property that
is tested symbolically.

## Usage

``` r
lyapunov_functional(model, n_sample_ic = 5L, verbose = TRUE, ...)
```

## Arguments

- model:

  A `system_spec` object with `is_pde = TRUE` (1D PDE; 2D is not yet
  supported at the functional level).

- n_sample_ic:

  Integer, number of random initial conditions at which the time
  derivative `dV/dt` is numerically evaluated to corroborate the
  algebraic argument. Defaults to 5.

- verbose:

  Logical. Narrate the construction. Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function` with `type = "functional"`.
The `params` field stores the symbolic potential `F_expr`, the diffusion
coefficient representation, the grid spacing `dx`, and a vector
`dV_dt_samples` of measured time derivatives at sample initial
conditions.

## Details

The gradient test compares \\\partial f_i / \partial u_j\\ across the
reaction term for mixed partial symmetry; if symmetric the reaction is a
gradient field and \\F\\ is recovered by integrating the formula \\-\int
f_1 \\ du_1\\. For scalar PDE (\\n = 1\\), any continuous \\f(u)\\
derives from a potential trivially. For systems, the test is strict.

When the PDE has state-dependent diffusion or non-gradient reaction, the
constructor returns `feasible = FALSE` with a clear diagnostic. The
functional is computed on the discretised state vector using central
differences for \\\nabla u\\ and the trapezoidal rule for the spatial
integral.

## References

Henry, D. (1981). *Geometric Theory of Semilinear Parabolic Equations*.
Springer. [doi:10.1007/BFb0089647](https://doi.org/10.1007/BFb0089647)

Temam, R. (1997). *Infinite-Dimensional Dynamical Systems in Mechanics
and Physics* (2nd ed.). Springer.
[doi:10.1007/978-1-4612-0645-3](https://doi.org/10.1007/978-1-4612-0645-3)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fisher <- system_spec(
    pde = list(u ~ D * d2x(u) + r * u * (1 - u)),
    state_names = "u",
    parms = list(D = 0.01, r = 1),
    spatial = list(domain = c(0, 1), n_grid = 51,
        bc = list(u = list(type = "neumann", left = 0, right = 0))),
    init = function(x) 0.5 + 0.1 * sin(pi * x)
)
lf <- lyapunov_functional(fisher)
print(lf); summary(lf)
} # }
```
