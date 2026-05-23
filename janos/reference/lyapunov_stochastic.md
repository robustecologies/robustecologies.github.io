# Stochastic Lyapunov function for an SDE

Constructs a quadratic stochastic Lyapunov function \\V(x) = (x -
x^\*)^\top P (x - x^\*)\\ for an Ito stochastic differential equation
\\dX = f(X) dt + g(X) dW\\ specified via
[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
with the `diffusion` argument. The drift \\f\\ is linearised at \\x^\*\\
to \\A = Df(x^\*)\\ and the continuous Lyapunov equation \\A^\top P + P
A = -Q\\ is solved. The infinitesimal generator satisfies
\\\mathcal{L}V(x) = \nabla V \cdot f + \tfrac{1}{2} \mathrm{tr}(g g^\top
\nabla^2 V)\\; exponential stability in mean square is certified when
the drift is Hurwitz and the diffusion is additive, and in probability
by Khasminskii's theorem for state-dependent diffusion.

## Usage

``` r
lyapunov_stochastic(
  model,
  x_star = NULL,
  Q = NULL,
  tol = sqrt(.Machine$double.eps),
  verbose = TRUE,
  ...
)
```

## Arguments

- model:

  A `system_spec` object with `is_sde = TRUE`.

- x_star:

  Numeric vector, the drift equilibrium. If `NULL`,
  [`equilibrium()`](https://robustecologies.github.io/janos/reference/equilibrium.md)
  is used with the model's defaults.

- Q:

  Symmetric positive-definite weight matrix. Defaults to the identity.

- tol:

  Hurwitz tolerance for the linearised drift. Defaults to
  `sqrt(.Machine$double.eps)`.

- verbose:

  Logical. When `TRUE`, narrates the construction with coloured symbols.
  Defaults to `TRUE`.

- ...:

  Unused.

## Value

An S3 object of class `lyapunov_function`. The `type` field is
`"stochastic"`. The `params` field additionally records the evaluated
diffusion matrix \\G(x^\*)\\ and the trace correction `tr_GtPG` that
appears in the generator.

## Details

When the diffusion is independent of the state (detected by the internal
`detect_additive_noise` helper), \\\mathcal{L}V(x) = -(x - x^\*)^\top Q
(x - x^\*) + \mathrm{tr}(G^\top P G)\\. The second term is a constant
that shifts the stationary distribution but does not prevent asymptotic
stability; the certificate label is `"algebraic (additive linear SDE)"`
in that case (Khasminskii, 2012, Chapter 5).

When the diffusion depends on the state, the advisor reports that the
certificate is `"local_algebraic (Khasminskii)"`: \\V\\ still certifies
local exponential stability in probability provided the drift
linearisation is Hurwitz and the diffusion is Lipschitz-bounded
(Khasminskii, 2012, Theorem 5.5). The decrease condition is verified
numerically at random points in a neighbourhood of \\x^\*\\.

## References

Khasminskii, R. (2012). *Stochastic Stability of Differential Equations*
(2nd ed.). Springer.
[doi:10.1007/978-3-642-23280-0](https://doi.org/10.1007/978-3-642-23280-0)

Mao, X. (2007). *Stochastic Differential Equations and Applications*
(2nd ed.). Woodhead Publishing. ISBN: 978-1-904275-34-3.

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
ou <- system_spec(
    rhs       = list(x ~ -theta * x),
    diffusion = list(x ~ sigma),
    state_names = "x",
    parms = list(theta = 1, sigma = 0.3),
    init = c(x = 0)
)
lf <- lyapunov_stochastic(ou)
print(lf); summary(lf)
} # }
```
