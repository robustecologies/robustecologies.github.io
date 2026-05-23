# Lyapunov spectrum via QR renormalisation

Computes the full Lyapunov spectrum \\\\\lambda_1 \ge \lambda_2 \ge
\dots \ge \lambda_n\\\\ of a continuous-time or discrete-time dynamical
system by propagating an orthonormal basis of tangent vectors alongside
the reference trajectory and periodically restoring its orthogonality
via QR decomposition. The method follows Benettin, Galgani, Giorgilli
and Strelcyn (1980) for flows and Shimada and Nagashima (1979) for maps.

## Usage

``` r
lyapunov_spectrum(
  model,
  t_max = 200,
  dt = 0.01,
  renorm_interval = NULL,
  discard_transient = NULL,
  init = NULL,
  parms = NULL,
  n_exponents = NULL,
  fd_eps = 1e-06,
  verbose = TRUE
)
```

## Arguments

- model:

  A `system_spec` describing an ODE or a discrete map.

- t_max:

  Total integration time (flows) or number of iterations (maps).

- dt:

  Sub-step for the RK4 integration of the combined (state, tangent)
  system. Ignored for maps.

- renorm_interval:

  Renormalisation period \\\tau\\ in time units (flows) or iterations
  (maps). Default: 1.0 for flows, 1 for maps.

- discard_transient:

  Time (flows) or iterations (maps) discarded before the accumulation of
  logarithms begins. Default: 10 percent of `t_max`.

- init:

  Optional initial condition overriding `model$init`.

- parms:

  Optional parameter list overriding `model$parms`.

- n_exponents:

  Number of leading exponents to compute. Default: full spectrum
  (`length(state_names)`).

- fd_eps:

  Relative step size for the finite-difference Jacobian.

- verbose:

  Logical, print progress.

## Value

An S3 object of class `lyapunov_spectrum` with components `exponents`
(numeric vector, sorted descending), `convergence` (matrix, running
estimates of each exponent at every renormalisation step),
`kaplan_yorke` (scalar), `sum` (sum of exponents, equal to the
time-average divergence of the flow), `family` ("ode" or "map"),
`t_max`, `renorm_interval`, `model`.

## Details

Numerical Lyapunov spectrum of an ODE or discrete-map model

Let \\\dot y = f(y)\\ and \\J(y) = \partial f / \partial y\\. The
variational equation \\\dot Q = J(y(t)) Q(t)\\ is integrated for an
orthonormal initial basis \\Q(0) = I\\. After every renormalisation
interval of length \\\tau\\ the matrix \\Q(t_k)\\ is factorised as \\Q
R\\ with \\Q\\ orthogonal and \\R\\ upper triangular with positive
diagonal; the growth of the \\i\\-th axis is accumulated as \\\log
R\_{ii}\\. After \\N\\ intervals the exponents are recovered as
\\\lambda_i = (1/N\tau) \sum_k \log R\_{ii}^{(k)}\\. For a discrete map
\\x\_{n+1} = F(x_n)\\ the same scheme runs with \\J(x_n)\\ and \\\tau\\
counted in iterations.

The Jacobian is evaluated by centred finite differences on the
system_spec callable RHS, so any formula-based ODE or map works without
additional compilation. The reference trajectory is integrated by a
fourth-order Runge-Kutta scheme for flows and by direct iteration for
maps; both co-evolve the tangent basis in lockstep.

The Kaplan-Yorke (Lyapunov) dimension is returned as a byproduct:
\\D\_{KY} = k + (\sum\_{i=1}^{k} \lambda_i) / \|\lambda\_{k+1}\|\\, with
\\k\\ the largest index such that the partial sum is non-negative.

## References

Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J.-M. (1980).
Lyapunov characteristic exponents for smooth dynamical systems and for
Hamiltonian systems; a method for computing all of them. *Meccanica*,
15(1), 9-30.
[doi:10.1007/BF02128236](https://doi.org/10.1007/BF02128236)

Shimada, I., & Nagashima, T. (1979). A numerical approach to the ergodic
problem of dissipative dynamical systems. *Progress of Theoretical
Physics*, 61(6), 1605-1616.
[doi:10.1143/PTP.61.1605](https://doi.org/10.1143/PTP.61.1605)

Kaplan, J. L., & Yorke, J. A. (1979). Chaotic behavior of
multidimensional difference equations. In *Functional Differential
Equations and Approximation of Fixed Points*, Lecture Notes in
Mathematics 730, Springer, pp. 204-227.
[doi:10.1007/BFb0064319](https://doi.org/10.1007/BFb0064319)

## See also

[`correlation_dimension`](https://robustecologies.github.io/janos/reference/correlation_dimension.md),
[`poincare_section`](https://robustecologies.github.io/janos/reference/poincare_section.md),
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md),
[`bifurcation_diagram`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md),
[`spectral_analysis`](https://robustecologies.github.io/janos/reference/spectral_analysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lorenz <- system_spec(
    rhs = list(
        x ~ sigma * (y - x),
        y ~ x * (rho - z) - y,
        z ~ x * y - beta * z
    ),
    state_names = c("x", "y", "z"),
    parms = list(sigma = 10, rho = 28, beta = 8 / 3),
    init = c(x = 1, y = 1, z = 1)
)
ls <- lyapunov_spectrum(lorenz, t_max = 200, dt = 0.01)
print(ls); summary(ls); plot(ls)
} # }
```
