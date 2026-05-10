# Pseudo-arclength bifurcation continuation along a parameter

Traces a branch of equilibria as a single parameter varies, using
pseudo-arclength continuation with bifurcation detection. The method
follows the solution curve through fold (saddle-node) bifurcations where
natural parameterization fails, and detects both fold and Hopf
bifurcation points along the branch.

## Usage

``` r
continuation(
  model,
  par_name,
  par_range,
  init_eq = NULL,
  ds = 0.01,
  max_points = 1000L,
  newton_tol = 1e-08,
  newton_max_iter = 20L,
  bifurcation_tol = 1e-06,
  detect_bifurcations = TRUE,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with formula-based RHS.

- par_name:

  Character string naming the parameter to continue in.

- par_range:

  Numeric vector of length 2 giving the parameter range
  `c(par_start, par_end)`.

- init_eq:

  Optional numeric vector or `equilibrium_point` object giving the
  starting equilibrium. If `NULL`, Newton's method is used to find an
  equilibrium at `par_range[1]`.

- ds:

  Arclength step size (default: 0.01). Smaller values give finer
  resolution but longer computation. The algorithm does not adapt step
  size automatically.

- max_points:

  Maximum number of continuation points (default: 1000).

- newton_tol:

  Convergence tolerance for the Newton corrector (default: 1e-8).

- newton_max_iter:

  Maximum Newton iterations per corrector step (default: 20).

- bifurcation_tol:

  Tolerance for bisection-based bifurcation location (default: 1e-6).

- detect_bifurcations:

  Logical; whether to detect fold and Hopf bifurcations (default: TRUE).

- verbose:

  Logical; print progress messages (default: TRUE).

## Value

An S3 object of class `continuation_result` containing:

- branch:

  Data frame with columns for arclength s, the continuation parameter,
  each state variable at equilibrium, the real parts of the eigenvalues,
  and a logical stability flag

- bifurcations:

  Data frame of detected bifurcation points with type, parameter value,
  state values, and critical eigenvalues. Empty if none detected.

- model:

  The model_spec used

- par_name:

  Name of the continuation parameter

- par_range:

  The requested parameter range

- ds:

  Step size used

## Details

The algorithm implements Keller's pseudo-arclength continuation. At each
point on the branch, the extended system is

\$\$F(y, p) = f(y, p) = 0\$\$ \$\$N(y, p) = (y - y_0)^T \dot{y}\_0 +
(p - p_0) \dot{p}\_0 - ds = 0\$\$

where \\(y_0, p_0)\\ is the current point, \\(\dot{y}\_0, \dot{p}\_0)\\
is the tangent direction, and \\ds\\ is the arclength step. The bordered
system is solved using the standard decomposition: solve \\J u = -f\\
and \\J v = -f_p\\, then assemble the correction from the arclength
constraint.

The tangent predictor at each step is computed from the null-space
direction of the augmented Jacobian \\\[J \| f_p\]\\. After prediction,
Newton correction on the augmented system brings the iterate back to the
branch.

Bifurcation detection monitors two indicators along the branch. A fold
(saddle-node) bifurcation is detected when \\dp/ds\\ changes sign,
equivalently when the determinant of J crosses zero. A Hopf bifurcation
is detected when a complex conjugate pair of eigenvalues crosses the
imaginary axis, identified by a sign change in the real part of the
rightmost complex eigenvalue pair. When a bifurcation is detected
between consecutive points, bisection on the arclength locates it within
`bifurcation_tol`.

The parameter derivative \\f_p = \partial f / \partial p\\ is computed
symbolically using the same differentiation engine that generates the
Jacobian, applied with respect to the continuation parameter rather than
a state variable.

## References

Keller, H. B. (1977). Numerical solution of bifurcation and nonlinear
eigenvalue problems. In *Applications of bifurcation theory*, pp.
359-384. Academic Press.

Seydel, R. (2010). *Practical bifurcation and stability analysis*.
Springer. ISBN: 978-1-4419-1739-3.

## See also

[`equilibrium`](https://robustecologies.github.io/janos/reference/equilibrium.md),
[`bifurcation_sweep`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md),
[`bifurcation_diagram`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Saddle-node bifurcation in normal form
sn <- model_spec(
    rhs = list(x ~ r + x^2),
    state_names = "x",
    parms = list(r = -1),
    init = c(x = -1)
)
result <- continuation(sn, par_name = "r",
                       par_range = c(-2, 0.5), ds = 0.01)
print(result)
summary(result)
plot(result)
plot(result, type = "eigenvalue")
} # }
```
