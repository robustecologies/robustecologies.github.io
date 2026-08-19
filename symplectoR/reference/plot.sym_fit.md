# Plot a symplectic optimization fit

Diagnostic views of a
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
run. The default `"trace"` view shows the objective against iteration on
a logarithmic scale, the single most informative summary of solver
behaviour; `"dashboard"` assembles the informative views of the
particular fit into one figure.

## Usage

``` r
# S3 method for class 'sym_fit'
plot(
  x,
  type = c("trace", "path", "phase", "energy", "density", "gradient", "increments",
    "samples", "observed", "residuals", "dashboard"),
  ...
)
```

## Arguments

- x:

  A `sym_fit` object.

- type:

  View: `"trace"` (objective and gradient norm against iteration),
  `"path"` (iterate path over the objective contours, dimension 2 only),
  `"phase"` (position against a finite-difference velocity proxy, first
  coordinate), `"energy"` (the Lyapunov energy recorded by the leapfrog
  kernels; monotone decay indicates a faithful integration), `"density"`
  (final probability density marginals, quantum Hamiltonian descent
  only), `"gradient"` (gradient-norm decay), `"increments"` (per-step
  objective decrement), `"samples"` (measured objective values, quantum
  Hamiltonian descent only), `"observed"` (observations against the
  fitted predictions, for fits of objectives built by
  [`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md)),
  `"residuals"` (the corresponding residuals), or `"dashboard"` (a
  four-panel composite selecting the views available for this fit).

- ...:

  Ignored.

## Value

The ggplot object, invisibly; for `type = "dashboard"` the assembled
grid graphical object.

## See also

[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[`print.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/print.sym_fit.md),
[`summary.sym_fit()`](https://robustecologies.github.io/symplectoR/reference/summary.sym_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- sym_optim(sym_objective(sym_benchmark("quadratic", d = 2)),
                 x0 = c(-1, 2), method = "rgd", keep_path = "full")
plot(fit, type = "trace")
plot(fit, type = "path")
plot(fit, type = "dashboard")

## Fits of a sym_inverse() objective additionally carry the data views
main <- subset(nicholson_blowfly, population == "main")
N <- main$count; k <- 7L; idx <- seq.int(k + 1L, length(N) - 1L)
obj <- sym_inverse(function(th) log(exp(th[1]) * N[idx - k] *
                                    exp(-N[idx - k] / exp(th[2])) +
                                    N[idx] * exp(-exp(th[3])) + 1),
                   data = log(N[idx + 1L] + 1), obs_times = main$day[idx + 1L],
                   theta_bounds = list(lo = c(logP = -3, logN0 = 4, logdelta = -5),
                                       hi = c(logP = 4, logN0 = 10, logdelta = 1)))
bf <- sym_optim(obj, method = "qhd", seed = 1)
plot(bf, type = "observed")
plot(bf, type = "residuals")
plot(bf, type = "dashboard")
} # }
```
