# Equilibrium of an ODE system

Locates a steady-state solution of f(y, p) = 0 using Newton's method
with the symbolically compiled Jacobian from the model specification.
The equilibrium is classified by the eigenvalues of the Jacobian
evaluated at the fixed point.

## Usage

``` r
equilibrium(
  model,
  init_guess = NULL,
  parms = NULL,
  tol = 1e-10,
  max_iter = 50L
)
```

## Arguments

- model:

  A `system_spec` with formula-based RHS. The symbolic Jacobian is
  compiled automatically on first use.

- init_guess:

  Optional numeric vector of initial guess for the equilibrium. If
  `NULL`, uses the model's default initial conditions.

- parms:

  Optional named list of parameter values. If `NULL`, uses the model's
  default parameters.

- tol:

  Convergence tolerance for the infinity norm of f (default: 1e-10).

- max_iter:

  Maximum number of Newton iterations (default: 50).

## Value

An S3 object of class `equilibrium_point` containing:

- y:

  Named numeric vector of the equilibrium state

- eigenvalues:

  Complex vector of eigenvalues of J at the equilibrium

- stable:

  Logical; TRUE if all eigenvalues have Re(lambda) \< 0

- classification:

  Character string describing the equilibrium type

- jacobian:

  The Jacobian matrix evaluated at the equilibrium

- residual:

  Infinity norm of f at the equilibrium

- n_iter:

  Number of Newton iterations used

- model:

  The system_spec (with compiled Jacobian)

- parms:

  Parameter values used

## Details

Newton iteration proceeds as y(k+1) = y(k) - inv(J(y(k))) \* f(y(k))
where J is the Jacobian matrix df/dy. Convergence is quadratic near the
solution. The iteration terminates when the infinity norm of f drops
below `tol` or `max_iter` iterations are reached. When convergence
fails, the function raises an error rather than returning an unreliable
result.

Stability classification uses the eigenvalues of J at the equilibrium. A
fixed point is stable when all eigenvalues have negative real parts. The
classification distinguishes stable nodes (all real, negative), unstable
nodes (all real, at least one positive), saddle points (real eigenvalues
with mixed signs), stable foci (complex pair with negative real part),
unstable foci (complex pair with positive real part), and centers
(purely imaginary eigenvalues, within numerical tolerance).

## References

Seydel, R. (2010). *Practical bifurcation and stability analysis*.
Springer. ISBN: 978-1-4419-1739-3.

## See also

[`continuation`](https://robustecologies.github.io/janos/reference/continuation.md),
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Lotka-Volterra equilibrium
lv <- system_spec(
    rhs = list(x ~ a * x - b * x * y, y ~ d * x * y - c * y),
    state_names = c("x", "y"),
    parms = list(a = 1.0, b = 0.1, d = 0.075, c = 1.5),
    init = c(x = 20, y = 10)
)
eq <- equilibrium(lv)
print(eq)
summary(eq)
} # }
```
