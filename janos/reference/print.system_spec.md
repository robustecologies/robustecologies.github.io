# Print method for a system_spec object

Formats a `system_spec` object as a compact one-screen synopsis: the
detected model type (ODE, DDE, SDE, jump-diffusion, CTMC, SSA, PDMP,
1D/2D PDE, 1D RDME, graph RDME), the state and parameter names, the
default initial condition, the positivity-clamp set, the compiler
backend, and any environmental or demographic noise attached to the
model.

## Usage

``` r
# S3 method for class 'system_spec'
print(x, ...)
```

## Arguments

- x:

  A `system_spec` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## See also

[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
. constructor;
[`summary.system_spec()`](https://robustecologies.github.io/janos/reference/summary.system_spec.md)
. tabulated states, parameters and equations;
[`plot.system_spec()`](https://robustecologies.github.io/janos/reference/plot.system_spec.md)
. schematic of the model (vector field, reaction graph or spatial grid);
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
. downstream simulator.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
print(m)
} # }
```
