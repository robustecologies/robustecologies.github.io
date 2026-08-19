# The non-smooth global optimization benchmark suite

Names and ground truth of the 12 non-smooth, non-convex, box-constrained
test functions of Leng et al. (2025) shipped through
[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md).
Each carries its published domain, global minimizer and minimum value;
none is smooth, and several are non-differentiable at or near the
minimizer, which is the regime where quantum Hamiltonian descent carries
its global convergence guarantee and subgradient methods can fail.

## Usage

``` r
sym_benchmark_suite()
```

## Value

A data frame with one row per suite function: `name`, `dim`, `f_star`,
and the box bounds pasted as text.

## Details

The suite members, all callable as `sym_benchmark("<name>")`: `wf`,
`crownedcross`, `bukin06`, `keane`, `schwefel` (one-dimensional),
`ackley`, `xinsheyang04`, `carromtable`, `rana` (two-dimensional except
as noted), `dropwave`, `layeb04` (three-dimensional), `damavandi`.
Domains are the published ones, several asymmetric (`bukin06` on
`[-15, -5] x [-3, 3]`, `crownedcross` on `[-10, 15]^2`, `ackley` on
`[-15, 30]^2`, `damavandi` on `[0, 14]^2`).

## References

Leng, J., Zheng, Y., Jia, Z., Fan, L., Zhao, C., Peng, Y., & Wu, X.
(2025). Quantum Hamiltonian descent for non-smooth optimization. *arXiv
preprint* arXiv:2503.15878.
[doi:10.48550/arXiv.2503.15878](https://doi.org/10.48550/arXiv.2503.15878)

## See also

[`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md),
[`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sym_benchmark_suite()
bm <- sym_benchmark("xinsheyang04")
obj <- sym_objective(bm)
fit <- sym_optim(obj, method = "qhd", seed = 1,
                 control = sym_control("qhd", N_grid = 256, K = 2000))
print(fit)
plot(fit, type = "density")
plot(fit, type = "dashboard")
} # }
```
