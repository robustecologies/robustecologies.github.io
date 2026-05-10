# Construct a Lyapunov function for an ODE dynamical system

Given a continuous-time dynamical system \\\dot x = f(x)\\ and
(optionally) an equilibrium point \\x^\*\\, constructs a Lyapunov
function \\V(x)\\ certifying asymptotic stability. Eight computational
strategies are available, selected automatically from the structural
type detected on the input or chosen explicitly via the `method`
argument. The input is uniformly a
[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
object; two ergonomic shortcuts are kept for numerical work that
bypasses symbolic specification.

## Usage

``` r
lyapunov_function(
  model,
  x_star = NULL,
  method = c("auto", "quadratic", "goh", "macarthur", "gilpin", "sos", "rbf", "massera",
    "cpa"),
  degree = 4L,
  domain = NULL,
  Q = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- model:

  A
  [`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
  object, a square numeric matrix \\A\\ interpreted as the linear system
  \\\dot x = Ax\\, or a named list with elements `r` and `alpha`
  interpreted as a generalized Lotka-Volterra system \\\dot x_i =
  x_i(r_i + \sum_j \alpha\_{ij} x_j)\\. Both shortcuts construct the
  same internal envelope that a `model_spec` would produce after
  structural detection, so the downstream algebra is identical.

- x_star:

  Numeric vector, the equilibrium point. If `NULL`, it is computed
  automatically: the origin for linear systems, the interior fixed point
  \\x^\* = -\alpha^{-1} r\\ for gLV systems, and a Newton solve via
  [`equilibrium()`](https://robustecologies.github.io/janos/reference/equilibrium.md)
  for general `model_spec` inputs whose RHS is parsed.

- method:

  Character, one of `"auto"` (default), `"quadratic"`, `"goh"`,
  `"macarthur"`, `"gilpin"`, `"sos"`, `"rbf"`, `"massera"` or `"cpa"`.
  See Details for the dispatch logic of `"auto"` and the theorem
  realised by each method.

- degree:

  Integer, polynomial degree for the SOS method. Defaults to 4. Ignored
  by other methods.

- domain:

  Numeric matrix with 2 rows and \\n\\ columns giving the lower and
  upper bounds of the search domain for RBF, Massera and CPA methods.
  Ignored by other methods.

- Q:

  Symmetric positive definite weight matrix for the quadratic method,
  defining the right-hand side of the Lyapunov equation \\A^\top P + PA
  = -Q\\. Defaults to the identity.

- verbose:

  Logical; when `TRUE`, the dispatcher and the method-specific solver
  narrate their progress via coloured status symbols.

- ...:

  Additional arguments forwarded to the method-specific solver (e.g.
  `n_collocation` for RBF, `T_horizon` and `n_eval` for Massera and
  Gilpin).

## Value

An S3 object of class `lyapunov_function` with fields:

- `V`:

  Function evaluating \\V(x)\\ at a point.

- `grad_V`:

  Function evaluating \\\nabla V(x)\\.

- `V_dot`:

  Function evaluating \\\dot V(x) = \nabla V \cdot f(x)\\.

- `type`:

  Character naming the method used.

- `params`:

  Method-specific parameters, for instance `P` and `condition_P` for
  quadratic, `c_weights` for Goh, the SDP dual for SOS, collocation
  weights for RBF.

- `x_star`:

  The equilibrium point.

- `system`:

  The internal envelope (`$f`, `$jacobian`, `$type`, `$n`, `$params`,
  `$state_names`) produced from the input model. Not an S3 object;
  exposed for inspection only.

- `doa`:

  Domain-of-attraction estimate, with fields `type` (`"global"`,
  `"local"`, `"positive_orthant"`), `level`, `volume` and an explanatory
  `note`.

- `certificate`:

  Certification details: `V_posdef`, `V_dot_negdef`, `method_details`.

- `residual`:

  Numerical residual of the defining condition.

- `elapsed`:

  Computation time in seconds.

## Details

The `"auto"` dispatch logic selects the method based on the structural
type detected on the input and the dimension \\n\\ of the state space:

- linear: `"quadratic"`.

- gLV: `"goh"` first (the ecologically canonical Lyapunov function);
  `"macarthur"` second (near-symmetric interactions); `"quadratic"`
  (linearisation) as the last resort.

- polynomial with \\n \le 6\\: `"sos"`.

- general with \\n \le 3\\: `"rbf"`.

- general with \\n \le 5\\: `"massera"`.

- otherwise: `"quadratic"` (local certificate by Lyapunov's indirect
  method).

**Quadratic** (`method = "quadratic"`). Solves the continuous Lyapunov
equation \\A^\top P + PA = -Q\\ at the Jacobian \\A = Df(x^\*)\\ using
the Bartels-Stewart algorithm (\\O(n^3)\\). Returns \\V(x) = (x -
x^\*)^\top P (x - x^\*)\\. For nonlinear systems this is a local
certificate by Lyapunov's indirect method.

**Goh** (`method = "goh"`). For generalized Lotka-Volterra systems,
constructs the logarithmic Lyapunov function \$\$V(x) = \sum_i c_i
\bigl(x_i - x^\*\_i - x^\*\_i \ln(x_i / x^\*\_i)\bigr)\$\$ due to Goh
(1977). The weights \\c_i = \mathrm{diag}(D)\\ come from the
Volterra-Lyapunov stability certificate (the diagonal witness matrix
\\D\\ solving \\D\alpha + \alpha^\top D \prec 0\\). Certifies global
asymptotic stability in the positive orthant whenever VL-stability
holds.

**MacArthur** (`method = "macarthur"`). For gLV with near-symmetric
interactions, the quadratic form \\V(x) = -2 r^\top x - x^\top \alpha
x\\ is a Lyapunov function whenever \\\alpha + \alpha^\top\\ is negative
definite (MacArthur, 1969). Reports the skew-symmetric fraction of
\\\alpha\\ as a proxy for applicability.

**Gilpin** (`method = "gilpin"`). Constructs the line-integral Lyapunov
function \\V(x) = \int\_{x^\*}^x f(\xi) \cdot d\xi\\ along straight-line
paths in state space for gLV systems (Gilpin, 1974). Certifies \\\dot V
\le 0\\ by construction.

**SOS** (`method = "sos"`). For polynomial systems of dimension \\n \le
6\\, searches for a polynomial Lyapunov function of degree `degree` via
sum-of-squares programming (Parrilo, 2000). The SDP is solved through
CVXR; feasibility gives an algebraic certificate of asymptotic stability
on the whole polynomial level set.

**RBF** (`method = "rbf"`). Approximates a Lyapunov function as a radial
basis function expansion \\V(x) = \sum_k w_k \phi(\\x - \xi_k\\)\\ and
enforces the decrease condition at collocation points via quadratic
programming (Giesl, 2007).

**Massera** (`method = "massera"`). Numerically evaluates the Massera
converse construction \$\$V(x) = \int_0^T \\\varphi(t,x) - x^\*\\^2 \\
dt\$\$ by integrating trajectories forward in time (Massera, 1949).
Operational whenever the flow is locally asymptotically stable.

**CPA** (`method = "cpa"`). Constructs a continuous piecewise affine
Lyapunov function on a Delaunay triangulation of the search domain via
linear programming (Hafstein, 2004). Practical for \\n = 2\\ or \\3\\.

## References

Lyapunov, A. M. (1892). *The general problem of the stability of
motion*. Translation in: Int. J. Control, 55(3), 531-773 (1992).
[doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253)

Goh, B. S. (1977). Global stability in many-species systems. *The
American Naturalist*, 111(977), 135-143.
[doi:10.1086/283144](https://doi.org/10.1086/283144)

MacArthur, R. H. (1969). Species packing, and what competition
minimizes. *Proc. Natl. Acad. Sci. USA*, 64(4), 1369-1371.
[doi:10.1073/pnas.64.4.1369](https://doi.org/10.1073/pnas.64.4.1369)

Gilpin, M. E. (1974). A Liapunov function for competition communities.
*Journal of Theoretical Biology*, 44(1), 35-48.
[doi:10.1016/S0022-5193(74)80028-7](https://doi.org/10.1016/S0022-5193%2874%2980028-7)

Parrilo, P. A. (2000). *Structured Semidefinite Programs and
Semialgebraic Geometry Methods in Robustness and Optimization*. PhD
thesis, California Institute of Technology.

Giesl, P. (2007). *Construction of Global Lyapunov Functions Using
Radial Basis Functions*. Lecture Notes in Mathematics 1904, Springer.
[doi:10.1007/978-3-540-69909-5](https://doi.org/10.1007/978-3-540-69909-5)

Massera, J. L. (1949). On Liapounoff's conditions of stability. *Annals
of Mathematics*, 50(3), 705-721.
[doi:10.2307/1969558](https://doi.org/10.2307/1969558)

Hafstein, S. F. (2004). A constructive converse Lyapunov theorem on
exponential stability. *Discrete and Continuous Dynamical Systems*,
10(3), 657-678.
[doi:10.3934/dcds.2004.10.657](https://doi.org/10.3934/dcds.2004.10.657)

## See also

[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
— the canonical model constructor;
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
— higher-level family-aware dispatcher;
[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md)
— ranks methods by structural applicability;
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
— S3 trio for the returned object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear system via the matrix shortcut
A <- matrix(c(-1, 0.5, 0, -2), 2, 2)
lf <- lyapunov_function(A)
print(lf)
summary(lf)
plot(lf)

# Generalized Lotka-Volterra via the list shortcut (Goh method)
r <- c(1, 0.5)
alpha <- matrix(c(-1, -0.3, -0.2, -1), 2, 2)
lf_glv <- lyapunov_function(list(r = r, alpha = alpha))
print(lf_glv)
plot(lf_glv)

# Full model_spec ingress (the canonical form)
m <- model_spec(
    rhs = list(x ~ x * (1 - x - 0.3 * y),
               y ~ y * (0.8 - 0.2 * x - y)),
    state_names = c("x", "y"),
    parms = list(),
    init = c(x = 0.5, y = 0.5)
)
lf_m <- lyapunov_function(m)
plot(lf_m, type = "phase_portrait")
} # }
```
