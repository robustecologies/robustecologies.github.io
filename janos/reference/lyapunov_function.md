# Construct a Lyapunov function for a dynamical system

Given a continuous-time dynamical system \\\dot x = f(x)\\ and an
equilibrium point \\x^\*\\, constructs a Lyapunov function \\V(x)\\
certifying asymptotic stability. Multiple computational strategies are
available, selected automatically based on the system type or specified
explicitly via the `method` argument.

## Usage

``` r
lyapunov_function(
  system,
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

- system:

  A `dynamical_system` object (from
  [`dynamical_system`](https://robustecologies.github.io/janos/reference/dynamical_system.md)),
  a square numeric matrix \\A\\ (interpreted as the linear system \\\dot
  x = Ax\\), or a named list with elements `r` and `alpha` (interpreted
  as a generalized Lotka-Volterra system).

- x_star:

  Numeric vector, the equilibrium point. If `NULL`, computed
  automatically for linear systems (the origin) and gLV systems (\\x^\*
  = -\alpha^{-1} r\\).

- method:

  Character, one of `"auto"` (default), `"quadratic"`, `"goh"`, `"sos"`,
  `"rbf"`, `"massera"` or `"cpa"`. See Details for the dispatch logic of
  `"auto"`.

- degree:

  Integer, polynomial degree for the SOS method. Defaults to 4.

- domain:

  Numeric matrix with 2 rows and \\n\\ columns, giving the lower and
  upper bounds of the search domain for RBF, Massera and CPA methods.

- Q:

  Symmetric positive definite weight matrix for the quadratic method.
  Defaults to the identity.

- verbose:

  Logical; print progress information.

- ...:

  Additional arguments forwarded to the method-specific solver.

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

  Method-specific parameters (e.g. `P` for quadratic, `c_weights` for
  Goh).

- `x_star`:

  The equilibrium point.

- `system`:

  The `dynamical_system` object.

- `doa`:

  Domain of attraction estimate.

- `certificate`:

  Certification details: `V_posdef`, `V_dot_negdef`, `method_details`.

- `residual`:

  Numerical residual of the defining condition.

- `elapsed`:

  Computation time in seconds.

## Details

The `"auto"` dispatch logic selects the method based on the system type
and dimension:

- `"linear"`: always `"quadratic"`.

- `"glv"`: tries `"goh"` first (the ecologically canonical Lyapunov
  function); falls back to `"quadratic"` (linearisation) if VL-stability
  fails.

- `"polynomial"` with \\n \le 6\\: `"sos"`.

- `"general"` with \\n \le 3\\: `"rbf"`.

- `"general"` with \\n \le 5\\: `"massera"`.

- Otherwise: `"quadratic"` (linearisation).

**Quadratic** (`method = "quadratic"`). Solves the continuous Lyapunov
equation \\A^\top P + PA = -Q\\ at the Jacobian \\A = Df(x^\*)\\ using
the Bartels-Stewart algorithm (\\O(n^3)\\) for \\n \> 10\\ or the
Kronecker identity for smaller systems. Returns \\V(x) = (x - x^\*)^\top
P (x - x^\*)\\. For nonlinear systems this is a local certificate by
Lyapunov's indirect method.

**Goh** (`method = "goh"`). For generalized Lotka-Volterra systems,
constructs the logarithmic Lyapunov function \\V(x) = \sum_i c_i (x_i -
x^\*\_i - x^\*\_i \ln(x_i / x^\*\_i))\\ due to Goh (1977). The weights
\\c_i\\ are obtained from the Volterra-Lyapunov stability certificate
(the diagonal witness matrix \\D\\ from the internal VL-stability
solver). When VL-stability holds, this function certifies global
asymptotic stability in the positive orthant.

**SOS** (`method = "sos"`). For polynomial systems, searches for a
polynomial Lyapunov function of degree `degree` via sum-of-squares
programming (Parrilo, 2000). Requires CVXR.

**RBF** (`method = "rbf"`). Approximates a Lyapunov function as a radial
basis function expansion and enforces the decrease condition at
collocation points via quadratic programming.

**Massera** (`method = "massera"`). Numerically evaluates the Massera
converse construction \\V(x) = \int_0^T \\\varphi(t,x) - x^\*\\^2\\ dt\\
by integrating trajectories forward.

**CPA** (`method = "cpa"`). Constructs a continuous piecewise affine
Lyapunov function on a Delaunay triangulation via linear programming.
Practical for \\n = 2\\ or \\3\\.

## References

Lyapunov, A. M. (1892). *The general problem of the stability of
motion*. Translation in: Int. J. Control, 55(3), 531-773 (1992).
[doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253)

Goh, B. S. (1977). Global stability in many-species systems. *The
American Naturalist*, 111(977), 135-143.
[doi:10.1086/283144](https://doi.org/10.1086/283144)

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

[`dynamical_system()`](https://robustecologies.github.io/janos/reference/dynamical_system.md),
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear system
A <- matrix(c(-1, 0.5, 0, -2), 2, 2)
lf <- lyapunov_function(A)
print(lf)
summary(lf)
plot(lf)

# Generalized Lotka-Volterra (Goh method)
r <- c(1, 0.5)
alpha <- matrix(c(-1, -0.3, -0.2, -1), 2, 2)
lf_glv <- lyapunov_function(list(r = r, alpha = alpha))
print(lf_glv)
plot(lf_glv)

# Forced quadratic on gLV
lf_quad <- lyapunov_function(list(r = r, alpha = alpha), method = "quadratic")
} # }
```
