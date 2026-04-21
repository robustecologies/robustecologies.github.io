# Unified Lyapunov analysis for a model_spec

Analyses a
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
object and, when theory allows, constructs a Lyapunov function (or
functional) with an algebraic or local-algebraic certificate. The entry
point works uniformly for continuous ODEs, discrete maps, stochastic
differential equations, delay differential equations, piecewise
deterministic Markov processes, continuous-time Markov chain reaction
networks and reaction-diffusion partial differential equations; the
family dispatcher is driven by
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md).

## Usage

``` r
analyse_lyapunov(
  model,
  method = c("auto", "quadratic", "goh", "macarthur", "gilpin", "sos", "rbf", "massera",
    "cpa", "discrete", "stochastic", "krasovskii", "pdmp", "foster", "functional"),
  x_star = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- model:

  A `model_spec` object.

- method:

  Character, one of `"auto"` (default) or a specific family method:
  `"quadratic"`, `"goh"`, `"macarthur"`, `"gilpin"`, `"sos"`, `"rbf"`,
  `"massera"`, `"cpa"` (all ODE), `"discrete"`, `"stochastic"`,
  `"krasovskii"`, `"pdmp"`, `"foster"`, `"functional"`.

- x_star:

  Optional numeric vector; the equilibrium or linearisation point.
  Defaults depend on the family.

- verbose:

  Logical. When `TRUE`, the advisor and the family constructor narrate
  each step of the reasoning using coloured symbols (internal
  `get_colored_symbol()` helper). Defaults to `TRUE`.

- ...:

  Additional arguments forwarded to the family constructor.

## Value

An S3 object of class `lyapunov_report` with fields:

- `advisor`:

  The `lyapunov_advisor` object returned by
  [`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md).

- `family`:

  The top-level family (ODE, map, SDE, DDE, PDMP, CTMC, PDE).

- `method`:

  The method actually used.

- `lyapunov`:

  The inner `lyapunov_function` object, or `NULL` when no construction
  was possible.

- `feasible`:

  Logical.

- `certificate_type`:

  One of `"algebraic"`, `"local_algebraic"`, `"numerical"`, `"none"`.

- `reason`:

  Human-readable rationale.

- `elapsed`:

  Seconds.

## Details

The dispatcher follows the plan emitted by
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md).
When `method = "auto"`, the top entry of `advisor$feasible_methods` is
selected; if multiple ODE methods are feasible, the cascade follows the
heuristic of
[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md):
quadratic for linear, Goh then MacArthur then quadratic for gLV, SOS for
polynomial within dimension and degree caps, and RBF, Massera or
quadratic otherwise. For non-ODE families the single family constructor
is called directly. When the advisor reports no feasible method,
`analyse_lyapunov` returns `feasible = FALSE` with the rationale, rather
than fabricating a certificate.

The `certificate_type` field summarises the epistemic status of the
result. `"algebraic"` means a closed-form proof (a solved Lyapunov
equation, a feasible LMI, a gradient-structure theorem);
`"local_algebraic"` means an algebraic proof that holds only in a
neighbourhood (Lyapunov's indirect method, Khasminskii local stability
in probability); `"numerical"` means the certificate is a numerical spot
check on a grid or sample.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

Khasminskii, R. (2012). *Stochastic Stability of Differential Equations*
(2nd ed.). Springer.
[doi:10.1007/978-3-642-23280-0](https://doi.org/10.1007/978-3-642-23280-0)

Gu, K., Kharitonov, V. L., & Chen, J. (2003). *Stability of Time-Delay
Systems*. Birkhauser.
[doi:10.1007/978-1-4612-0039-0](https://doi.org/10.1007/978-1-4612-0039-0)

Meyn, S. P., & Tweedie, R. L. (2009). *Markov Chains and Stochastic
Stability* (2nd ed.). Cambridge University Press.
[doi:10.1017/CBO9780511626630](https://doi.org/10.1017/CBO9780511626630)

Henry, D. (1981). *Geometric Theory of Semilinear Parabolic Equations*.
Springer. [doi:10.1007/BFb0089647](https://doi.org/10.1007/BFb0089647)

## See also

[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md),
[`lyapunov_discrete()`](https://robustecologies.github.io/janos/reference/lyapunov_discrete.md),
[`lyapunov_stochastic()`](https://robustecologies.github.io/janos/reference/lyapunov_stochastic.md),
[`lyapunov_krasovskii()`](https://robustecologies.github.io/janos/reference/lyapunov_krasovskii.md),
[`lyapunov_pdmp()`](https://robustecologies.github.io/janos/reference/lyapunov_pdmp.md),
[`lyapunov_foster()`](https://robustecologies.github.io/janos/reference/lyapunov_foster.md),
[`lyapunov_functional()`](https://robustecologies.github.io/janos/reference/lyapunov_functional.md),
[`print.lyapunov_report()`](https://robustecologies.github.io/janos/reference/print.lyapunov_report.md),
[`summary.lyapunov_report()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_report.md),
[`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(
    rhs = list(x ~ x * (1 - x - 0.3 * y),
               y ~ y * (0.8 - 0.2 * x - y)),
    state_names = c("x", "y"),
    parms = list(),
    init = c(x = 0.5, y = 0.5)
)
rep <- analyse_lyapunov(m, verbose = TRUE)
print(rep)
summary(rep)
plot(rep)
} # }
```
