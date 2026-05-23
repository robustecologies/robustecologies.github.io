# Plot method for a solver_spec object

Visualises the absolute-stability region of a `solver_spec` object in
the complex plane. For explicit Runge-Kutta methods (RK4, RK45) the
bounded region where \\\|R(z)\| \le 1\\ is shaded; for the implicit
Rosenbrock Rodas3 the A-stable left half-plane is shaded with the
imaginary axis as the boundary. Stochastic, event-driven and
method-of-lines solvers do not carry a conventional stability region;
for these, a schematic text card is rendered instead, naming the
algorithm and its event-loop or CFL constraint.

## Usage

``` r
# S3 method for class 'solver_spec'
plot(x, n_grid = 401L, title = NULL, ...)
```

## Arguments

- x:

  A `solver_spec` object.

- n_grid:

  Integer, grid resolution for the stability contour. Defaults to 401.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## Details

For an explicit Runge-Kutta method applied to \\\dot y = \lambda y\\
with complex \\\lambda\\ and step \\h\\, the numerical solution is
\\y\_{n+1} = R(z) y_n\\ with \\z = h\lambda\\. The absolute-stability
region is \\S = \\z \in \mathbb{C} : \|R(z)\| \le 1\\\\. The RK4
stability polynomial is \\R_4(z) = 1 + z + z^2/2 + z^3/6 + z^4/24\\ and
RK4(5) is plotted via its order-5 component for the main region.

## References

Hairer, E., Norsett, S. P., & Wanner, G. (1993). *Solving Ordinary
Differential Equations I* (2nd ed.). Springer.
[doi:10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)

Hairer, E., & Wanner, G. (1996). *Solving Ordinary Differential
Equations II* (2nd ed.). Springer.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)

## See also

[`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md)
and siblings . constructors;
[`print.solver_spec()`](https://robustecologies.github.io/janos/reference/print.solver_spec.md)
. compact header;
[`summary.solver_spec()`](https://robustecologies.github.io/janos/reference/summary.solver_spec.md)
. method class and order.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(solver_rk45())
plot(solver_rosenbrock())
} # }
```
