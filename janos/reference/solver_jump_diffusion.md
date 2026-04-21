# Jump-diffusion SDE solver

Creates a solver specification for Ito jump-diffusion processes of the
form dX = f(X)dt + g(X)dW + sum_k h_k(X, J_k) dN_k(t), where dN_k are
independent Poisson processes with state-dependent intensity. Per time
step, the drift and diffusion are handled by either Euler-Maruyama or
Milstein, then for each jump channel the number of arrivals is drawn
from Poisson(lambda_k \* dt) and jump sizes are applied.

## Usage

``` r
solver_jump_diffusion(
  dt = 0.01,
  method = "euler_maruyama",
  seed = 42L,
  subsample = 1L,
  dg_eps = 1e-06
)
```

## Arguments

- dt:

  Time step (default: 0.01).

- method:

  Diffusion integration scheme: "euler_maruyama" (default) or
  "milstein".

- seed:

  Random seed (default: 42).

- subsample:

  Store every nth step (default: 1).

- dg_eps:

  Relative perturbation for Milstein finite difference (default: 1e-6).
  Ignored for Euler-Maruyama.

## Value

A list of class `solver_spec` with `method = "jump_diffusion"`.

## References

Merton, R. C. (1976). Option pricing when underlying stock returns are
discontinuous. *Journal of Financial Economics*, 3(1-2), 125-144.
[doi:10.1016/0304-405X(76)90022-2](https://doi.org/10.1016/0304-405X%2876%2990022-2)

Cont, R., & Tankov, P. (2004). *Financial Modelling with Jump
Processes*. Chapman and Hall/CRC. ISBN: 978-1584884132.

## See also

[`solver_euler_maruyama`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
[`solver_milstein`](https://robustecologies.github.io/janos/reference/solver_milstein.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Merton jump-diffusion model
merton <- model_spec(
    rhs = list(x ~ mu * x),
    diffusion = list(x ~ sigma * x),
    jumps = list(
        x ~ list(intensity = ~ lambda,
                 size_distribution = "normal",
                 size_mean = ~ mu_J * x, size_sd = ~ sigma_J * x)
    ),
    state_names = "x",
    parms = list(mu = 0.05, sigma = 0.2, lambda = 1.0,
                 mu_J = -0.1, sigma_J = 0.05),
    init = c(x = 1.0)
)
result <- dyn_sim(merton, t_max = 5,
                  solver = solver_jump_diffusion(dt = 0.001),
                  discard_transient = 0)
print(result)
summary(result)
plot(result)
} # }
```
