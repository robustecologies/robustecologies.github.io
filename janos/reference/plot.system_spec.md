# Plot method for a system_spec object

Renders a schematic representation of a `system_spec` object, dispatched
on the detected model type. Two-dimensional ODE, SDE and DDE models
render the vector field as a grid of arrows at the default parameters;
higher-dimensional continuous models render an equation card listing the
governing equations. CTMC, SSA and RDME models render the
species-reaction bipartite graph (species as filled circles, reactions
as open squares, stoichiometry edges as arrows). PDE models render the
spatial grid and boundary conditions as a horizontal bar. Discrete maps
in 1D render a cobweb of the iterate.

## Usage

``` r
# S3 method for class 'system_spec'
plot(x, type = c("auto", "equations", "grid"), title = NULL, ...)
```

## Arguments

- x:

  A `system_spec` object.

- type:

  Character, the plot type. `"auto"` (default) dispatches on the model
  type. `"equations"` forces the text-card rendering of the governing
  equations. `"grid"` forces the vector-field grid for 2D continuous
  models.

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## Details

The vector-field grid is a 24 x 24 grid of arrows over a square region
centred on the default initial condition; for higher-dimensional models
the equation card is always the fallback because a vector field cannot
be drawn in more than two dimensions. The reaction-graph renderer places
reactions in a single column and species around them, with arrow heads
marking the net stoichiometry sign. The PDE spatial bar annotates the
grid resolution and the boundary condition type per state variable.

## See also

[`system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.md)
. constructor;
[`print.system_spec()`](https://robustecologies.github.io/janos/reference/print.system_spec.md)
. compact header;
[`summary.system_spec()`](https://robustecologies.github.io/janos/reference/summary.system_spec.md)
. tabulated states, parameters and equations.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ x * (1 - x - 0.3 * y),
                           y ~ y * (0.8 - 0.2 * x - y)),
                state_names = c("x", "y"), parms = list(),
                init = c(x = 0.5, y = 0.5))
plot(m)
} # }
```
