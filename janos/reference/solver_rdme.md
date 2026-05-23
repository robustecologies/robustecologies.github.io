# RDME solver (reaction-diffusion master equation)

Creates a solver specification for the reaction-diffusion master
equation (RDME) using Gillespie direct simulation over a 1D grid of
voxels. The algorithm tracks integer molecule counts per voxel and
handles both intra-voxel reactions (same propensities as the well-mixed
SSA) and inter-voxel diffusion hops at rate D_k / dx^2 per molecule.

## Usage

``` r
solver_rdme(seed = 42, output_dt = 1, max_events = 1e+08)
```

## Arguments

- seed:

  Random seed for reproducible simulations (default: 42).

- output_dt:

  Time interval between stored output snapshots (default: 1.0).

- max_events:

  Maximum number of events (reactions + hops) before stopping (default:
  1e8).

## Value

A list of class `solver_spec` with `method = "rdme"`.

## Details

The RDME discretizes a spatial domain into voxels, each containing an
independent copy of the reaction system coupled by nearest-neighbor
diffusion. For species k in voxel i with count X(k,i), the diffusion hop
rate to an adjacent voxel is D(k) \* X(k,i) / dx^2, where dx is the
voxel width and D(k) is the macroscopic diffusion coefficient. The
Gillespie direct algorithm selects the next event (reaction or hop)
across all voxels, with propensity recomputation restricted to affected
voxels and their neighbors.

## References

Elf, J. and Ehrenberg, M. (2004). Spontaneous separation of bi-stable
biochemical systems into spatial domains of opposite phases. *Systems
Biology*, 1(2), 230-236.
[doi:10.1049/sb:20045021](https://doi.org/10.1049/sb%3A20045021)

## See also

[`solver_ssa_direct`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md),
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bd <- system_spec(
    stoichiometry = list(birth = c(N = 1L), death = c(N = -1L)),
    propensities = list(birth ~ lambda * N, death ~ mu * N),
    state_names = "N",
    parms = list(lambda = 1.0, mu = 0.5, D_N = 0.1),
    init = c(N = 100),
    spatial = list(
        domain = c(0, 10), n_voxels = 50,
        diffusion_rates = list(N = ~ D_N),
        bc = list(type = "reflecting")
    )
)
result <- dyn_sim(bd, t_max = 20, solver = solver_rdme(),
                  discard_transient = 0)
print(result)
summary(result)
plot(result)
} # }
```
