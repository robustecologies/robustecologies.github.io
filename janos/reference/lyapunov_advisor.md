# Advise on feasible Lyapunov techniques for a model

Inspects a
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
object, classifies its dynamical family and enumerates the
Lyapunov-function construction techniques that are theoretically
applicable. The advisor does not build a Lyapunov function; it returns a
ranked plan that
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
later follows, and optionally narrates the reasoning step by step when
`verbose = TRUE`.

## Usage

``` r
lyapunov_advisor(
  model,
  verbose = TRUE,
  max_poly_degree = 6L,
  symmetry_tol = 0.3,
  ...
)
```

## Arguments

- model:

  A `system_spec` object.

- verbose:

  Logical. When `TRUE`, progress messages use coloured symbols (internal
  `get_colored_symbol()` helper) to report the inspection, the family
  detected, which theory applies, and which methods have been rejected
  and why. Defaults to `TRUE`.

- max_poly_degree:

  Integer. Upper bound on the polynomial degree at which the
  sum-of-squares method is declared feasible. Defaults to 6.

- symmetry_tol:

  Numeric tolerance for the interaction-matrix symmetry test used by the
  MacArthur method. Defaults to 0.3.

- ...:

  Additional arguments, currently unused.

## Value

An S3 object of class `lyapunov_advisor` with fields `family`
(character), `subtype` (character), `dim` (integer), `equilibrium`
(numeric or `NULL`), `feasible_methods` (character vector, ranked),
`rejected_methods` (named list of rejection reasons), `theory_notes`
(character vector), `warnings` (character vector), `details`
(method-specific diagnostics) and `model` (a reference back to the
input).

## Details

Classification proceeds in layers. The top-level family is read from the
`is_*` flags of the `system_spec`: ODE, map, SDE, DDE, jump-diffusion,
PDMP, CTMC (with or without spatial RDME extension), PDE. Within the ODE
family, internal structural detectors (linearity, gLV, polynomial
degree, symmetry, gradient field, additive noise) inspect the parsed
right-hand side and classify the subtype as linear, generalized
Lotka-Volterra, polynomial (with degree), gradient-field or general
nonlinear.

The advisor then associates each subtype with the theorems that apply,
following the decision cascade of
[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
for ODE subtypes and branching to family-specific constructions for the
rest: the discrete Lyapunov equation \\A^\top P A - P = -Q\\ for maps
(Khalil, 2002); the stochastic Lyapunov equation \\A^\top P + PA = -Q\\
with Ito correction \\\mathcal{L}V = \nabla V \cdot f +
\tfrac{1}{2}\mathrm{tr}(g g^\top \nabla^2 V)\\ for SDEs (Khasminskii,
2012); Lyapunov-Krasovskii functionals solved by LMI for DDEs (Gu,
Kharitonov and Chen, 2003); common or mode-dependent quadratic LMIs for
PDMPs (Costa, Fragoso and Marques, 2005); the Foster-Lyapunov drift
criterion lifted from the fluid-limit ODE for CTMCs (Meyn and Tweedie,
2009); and the energy functional \\V\[u\] = \int (D/2)\\\nabla u\\^2 +
F(u)\\dx\\ for reaction-diffusion PDEs with gradient reaction (Henry,
1981; Temam, 1997).

When no construction applies (chaotic attractors without an equilibrium,
pure RDME processes, genuinely non-gradient PDEs, etc.) the advisor
returns an empty `feasible_methods` vector and a populated `warnings`
field explaining why.
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
then returns `feasible = FALSE` rather than fabricating a certificate.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

Khasminskii, R. (2012). *Stochastic Stability of Differential Equations*
(2nd ed.). Springer.
[doi:10.1007/978-3-642-23280-0](https://doi.org/10.1007/978-3-642-23280-0)

Gu, K., Kharitonov, V. L., & Chen, J. (2003). *Stability of Time-Delay
Systems*. Birkhauser.
[doi:10.1007/978-1-4612-0039-0](https://doi.org/10.1007/978-1-4612-0039-0)

Costa, O. L. V., Fragoso, M. D., & Marques, R. P. (2005). *Discrete-Time
Markov Jump Linear Systems*. Springer.
[doi:10.1007/b138575](https://doi.org/10.1007/b138575)

Meyn, S. P., & Tweedie, R. L. (2009). *Markov Chains and Stochastic
Stability* (2nd ed.). Cambridge University Press.
[doi:10.1017/CBO9780511626630](https://doi.org/10.1017/CBO9780511626630)

Henry, D. (1981). *Geometric Theory of Semilinear Parabolic Equations*.
Springer. [doi:10.1007/BFb0089647](https://doi.org/10.1007/BFb0089647)

Temam, R. (1997). *Infinite-Dimensional Dynamical Systems in Mechanics
and Physics* (2nd ed.). Springer.
[doi:10.1007/978-1-4612-0645-3](https://doi.org/10.1007/978-1-4612-0645-3)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md),
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`print.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/print.lyapunov_advisor.md),
[`summary.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_advisor.md),
[`plot.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
glv <- system_spec(
    rhs = list(x ~ x * (1 - x - 0.3 * y),
               y ~ y * (0.8 - 0.2 * x - y)),
    state_names = c("x", "y"),
    parms = list(),
    init = c(x = 0.5, y = 0.5)
)
adv <- lyapunov_advisor(glv, verbose = TRUE)
print(adv)
summary(adv)
plot(adv)
} # }
```
