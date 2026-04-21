# Plot method for map_portrait objects

Renders a map portrait as a 2D displacement-field and fixed-point
diagram (ggplot2) or as a cobweb diagram for 1D maps.

## Usage

``` r
# S3 method for class 'map_portrait'
plot(
  x,
  type = "default",
  title = NULL,
  arrow_scale = 1,
  arrow_type = c("open", "closed"),
  feasible = FALSE,
  ...
)
```

## Arguments

- x:

  A `map_portrait` object

- type:

  Character string selecting which components to render. One of
  `"default"` (all computed components), `"field"` (displacement field
  only), `"isoclines"`, `"fixed_points"`, `"manifolds"`, or `"scatter"`.

- title:

  Optional plot title. If `NULL`, uses the model name.

- arrow_scale:

  Numeric scaling factor for displacement arrows (default 1).

- arrow_type:

  Arrow head style: `"open"` (default) or `"closed"`.

- feasible:

  Controls shading of infeasible state-space regions. `FALSE` (default)
  disables shading. `TRUE` shades regions where any state variable is
  negative. A named list (e.g. `list(x = c(0, 1))`) shades regions
  outside the specified per-state bounds; states not listed are
  unconstrained.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.

## Details

Fixed points are encoded by both shape and color following the same
convention as
[`plot.phase_portrait`](https://robustecologies.github.io/janos/reference/plot.phase_portrait.md),
except that "center" is replaced by "non-hyperbolic" (open green
circle):

- Stable node: filled blue circle

- Unstable node: filled red circle

- Stable focus: filled blue diamond

- Unstable focus: filled red diamond

- Non-hyperbolic: open green circle

- Saddle: filled orange triangle

- Saddle-focus: filled orange inverted triangle

Manifold branches are colored blue (stable) or red (unstable).

## Examples

``` r
if (FALSE) { # \dontrun{
logistic <- model_spec(
    map = list(x ~ r * x * (1 - x)),
    state_names = "x",
    parms = list(r = 3.2),
    init = c(x = 0.1)
)
mp <- map_portrait(logistic, xlim = c(0, 1))
plot(mp)
} # }
```
