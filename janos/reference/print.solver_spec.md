# Print method for a solver_spec object

Formats a `solver_spec` object as a compact one-screen header naming the
numerical method (Dormand-Prince RK4(5), fixed-step RK4, Rosenbrock
Rodas3, Gillespie direct / NRM / modified NRM, adaptive tau-leap, hybrid
SSA/CLE, Euler-Maruyama, Milstein, jump-diffusion, Lewis-Shedler PDMP,
method-of-lines, RDME) and its key control parameters (step size,
tolerances, seeds, maximum reactions).

## Usage

``` r
# S3 method for class 'solver_spec'
print(x, ...)
```

## Arguments

- x:

  A `solver_spec` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## See also

[`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md),
[`solver_rk4()`](https://robustecologies.github.io/janos/reference/solver_rk4.md),
[`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md),
[`solver_map()`](https://robustecologies.github.io/janos/reference/solver_map.md),
[`solver_dde()`](https://robustecologies.github.io/janos/reference/solver_dde.md),
[`solver_euler_maruyama()`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md),
[`solver_milstein()`](https://robustecologies.github.io/janos/reference/solver_milstein.md),
[`solver_jump_diffusion()`](https://robustecologies.github.io/janos/reference/solver_jump_diffusion.md),
[`solver_ssa_direct()`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`solver_tau_leap()`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md),
[`solver_hybrid()`](https://robustecologies.github.io/janos/reference/solver_hybrid.md),
[`solver_pdmp()`](https://robustecologies.github.io/janos/reference/solver_pdmp.md),
[`solver_rdme()`](https://robustecologies.github.io/janos/reference/solver_rdme.md),
[`solver_mol()`](https://robustecologies.github.io/janos/reference/solver_mol.md),
[`solver_mol2d()`](https://robustecologies.github.io/janos/reference/solver_mol2d.md)
. constructors;
[`summary.solver_spec()`](https://robustecologies.github.io/janos/reference/summary.solver_spec.md)
. method class and order;
[`plot.solver_spec()`](https://robustecologies.github.io/janos/reference/plot.solver_spec.md)
. absolute-stability region.

## Examples

``` r
if (FALSE) { # \dontrun{
print(solver_rk45())
print(solver_ssa_direct())
} # }
```
