# Plot method for an equilibrium_point object

Renders the Jacobian spectrum at an equilibrium as a scatter of
eigenvalues in the complex plane. Points are coloured by the sign of the
real part, the imaginary axis separates stable (left) from unstable
(right) eigenvalues for continuous systems, and the unit circle is
overlaid when `type = "discrete"` to mark the stability boundary of a
discrete-time map.

## Usage

``` r
# S3 method for class 'equilibrium_point'
plot(x, type = c("continuous", "discrete"), title = NULL, ...)
```

## Arguments

- x:

  An `equilibrium_point` object.

- type:

  Character, either `"continuous"` (default) or `"discrete"`. Continuous
  systems use the imaginary axis as the stability boundary; discrete
  systems use the unit circle.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## Details

For a continuous-time ODE \\\dot x = f(x)\\, the equilibrium \\x^\*\\ is
asymptotically stable if every eigenvalue of \\Df(x^\*)\\ has negative
real part. The corresponding plot uses the left half-plane as the stable
region. For a discrete-time map \\x\_{n+1} = F(x_n)\\, asymptotic
stability requires every eigenvalue of \\DF(x^\*)\\ to lie strictly
inside the unit circle, and the plot overlays the unit circle to mark
the boundary. The stability classification from the classifier (stable
node, stable focus, saddle, centre, etc.) appears in the subtitle.

## References

Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos* (2nd ed.).
Westview Press. ISBN: 978-0813349107.

Kuznetsov, Y. A. (2004). *Elements of Applied Bifurcation Theory* (3rd
ed.). Springer.
[doi:10.1007/978-1-4757-3978-7](https://doi.org/10.1007/978-1-4757-3978-7)

## See also

[`equilibrium()`](https://robustecologies.github.io/janos/reference/equilibrium.md)
— constructor;
[`print.equilibrium_point()`](https://robustecologies.github.io/janos/reference/print.equilibrium_point.md)
— compact header;
[`summary.equilibrium_point()`](https://robustecologies.github.io/janos/reference/summary.equilibrium_point.md)
— eigenvalue table with condition number;
[`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
— pseudo-arclength continuation along a parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(
    rhs = list(x ~ x * (1 - x - 0.3 * y),
               y ~ y * (0.8 - 0.2 * x - y)),
    state_names = c("x", "y"), parms = list(),
    init = c(x = 0.5, y = 0.5)
)
eq <- equilibrium(m)
plot(eq)
} # }
```
