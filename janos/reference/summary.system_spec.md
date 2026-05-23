# Summary method for a system_spec object

Produces a structured summary of a `system_spec` object: detected model
type, tabulated state names and default initial conditions, tabulated
parameter names and their default values, a flattened textual listing of
the governing equations (RHS formulas for continuous models, map
formulas for discrete maps, reactions and propensities for CTMC and RDME
models, drift and diffusion for SDE, delay structure for DDE, per-regime
drift and transitions for PDMP, spatial operators and boundary
conditions for PDE), the compiled backend status, and the
positivity-clamp specification.

## Usage

``` r
# S3 method for class 'system_spec'
summary(object, ...)
```

## Arguments

- object:

  A `system_spec` object.

- ...:

  Unused.

## Value

A list of class `summary.system_spec`, invisibly.

## Details

The summary is a faithful extract of the parsed model fields; it does
not re-parse or compile the model. For reaction networks the
stoichiometry matrix is reported by its species-by-reaction dimensions;
for spatial models the grid resolution and boundary conditions are
reported on one line per state variable.

## See also

[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
. constructor;
[`print.system_spec()`](https://robustecologies.github.io/janos/reference/print.system_spec.md)
. compact header;
[`plot.system_spec()`](https://robustecologies.github.io/janos/reference/plot.system_spec.md)
. schematic of the model.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
summary(m)
} # }
```
