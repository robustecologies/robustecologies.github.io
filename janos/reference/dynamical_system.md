# Specify a dynamical system for Lyapunov analysis

Constructs an S3 object of class `dynamical_system` that wraps a
continuous-time ODE \\\dot x = f(x)\\ in a form suitable for Lyapunov
function construction. Four system types are supported: `"linear"`
(\\f(x) = Ax\\), `"glv"` (generalized Lotka-Volterra), `"polynomial"`
and `"general"`.

## Usage

``` r
dynamical_system(
  f = NULL,
  jacobian = NULL,
  type = c("general", "linear", "glv", "polynomial"),
  r = NULL,
  alpha = NULL,
  A = NULL,
  coefficients = NULL,
  n = NULL,
  state_names = NULL
)
```

## Arguments

- f:

  A function taking a numeric vector `x` of length \\n\\ and returning
  \\\dot x\\ as a numeric vector of the same length. Required for
  `type = "general"`; auto-generated for `"linear"` and `"glv"`.

- jacobian:

  Optional function taking `x` and returning the \\n \times n\\ Jacobian
  matrix \\\partial f / \partial x\\. If `NULL`, a numerical Jacobian
  via numDeriv is used when needed.

- type:

  Character, one of `"general"` (default), `"linear"`, `"glv"` or
  `"polynomial"`.

- r:

  Numeric vector of intrinsic growth rates (length \\S\\). Required when
  `type = "glv"`.

- alpha:

  Square interaction matrix (\\S \times S\\). Required when
  `type = "glv"`.

- A:

  Square system matrix. Required when `type = "linear"`.

- coefficients:

  For `type = "polynomial"`: a list describing the polynomial vector
  field. Each element corresponds to one component \\f_i\\ and is itself
  a list with entries `indices` (integer matrix, each row a multi-index)
  and `values` (numeric vector of coefficients).

- n:

  State dimension. Inferred from other arguments when possible.

- state_names:

  Optional character vector of names for the state variables.

## Value

An S3 object of class `dynamical_system` with fields `$f`, `$jacobian`,
`$type`, `$n`, `$params`, `$state_names`.

## Details

For `type = "glv"` the vector field is \\f_i(x) = x_i (r_i + \sum_j
\alpha\_{ij} x_j)\\ and the analytic Jacobian is \\J\_{ij}(x) =
\delta\_{ij}(r_i + \sum_k \alpha\_{ik} x_k) + x_i \alpha\_{ij}\\, or in
matrix form \\J(x) = \mathrm{diag}(r + \alpha x) + \mathrm{diag}(x)
\alpha\\. Both are generated automatically so the user need only supply
`r` and `alpha`.

For `type = "linear"` the Jacobian is constant and equal to \\A\\.

## References

Hofbauer, J., & Sigmund, K. (1998). *Evolutionary Games and Population
Dynamics*. Cambridge University Press.
[doi:10.1017/CBO9781139173179](https://doi.org/10.1017/CBO9781139173179)

Takeuchi, Y. (1996). *Global Dynamical Properties of Lotka-Volterra
Systems*. World Scientific. ISBN: 978-981-02-2471-4.

## See also

[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md),
[`print.dynamical_system()`](https://robustecologies.github.io/janos/reference/print.dynamical_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear system
A <- matrix(c(-1, 0.5, 0, -2), 2, 2)
sys_lin <- dynamical_system(A = A, type = "linear")
print(sys_lin)

# Generalized Lotka-Volterra
r <- c(1, 0.5)
alpha <- matrix(c(-1, -0.3, -0.2, -1), 2, 2)
sys_glv <- dynamical_system(r = r, alpha = alpha, type = "glv")
sys_glv$f(c(0.5, 0.5))

# General nonlinear (Van der Pol, reversed for stability)
sys_vdp <- dynamical_system(
  f = function(x) c(-x[2], x[1] - (1 - x[1]^2) * x[2]),
  n = 2L, type = "general"
)
} # }
```
