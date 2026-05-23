# Plot method for a lyapunov_report object

Visualises a `lyapunov_report`. Two report-level plot types are
available: `"certificate_stack"` renders a three-panel summary of the
algebraic certificate (positivity of V, negativity of dV/dt, residual of
the Lyapunov equation) and `"timeline"` shows the wall-clock sequence of
advisor, classifier, constructor and verifier. Any other value of `type`
is forwarded to
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
on the inner Lyapunov function. When the report is infeasible the
advisor plan is rendered instead.

## Usage

``` r
# S3 method for class 'lyapunov_report'
plot(x, type = "level_sets", ...)
```

## Arguments

- x:

  A `lyapunov_report` object.

- type:

  Character. Report-level types are `"certificate_stack"` and
  `"timeline"`. Any other value (default `"level_sets"`) is forwarded to
  [`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md).
  Family-specific types such as `"cobweb"`, `"generator_field"`,
  `"ensemble_decay"`, `"noise_ellipse"`, `"lmi_spectrum"`,
  `"delay_margin"`, `"regime_lmi"`, `"switched_trajectory"`,
  `"drift_grid"`, `"fluid_vs_ctmc"`, `"energy_decay"`,
  `"gradient_field_check"` and `"profile"` are recognised by the inner
  Lyapunov function via family dispatch.

- ...:

  Additional arguments forwarded to the delegated plot method.

## Value

A ggplot or plotly object.

## Details

Dispatch rule: if `type` is a report-level type (`"certificate_stack"`,
`"timeline"`), the corresponding panel is drawn from the report itself.
Otherwise, if `x$feasible` is `TRUE` and the inner Lyapunov function
exists, the request is forwarded to
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
with its family-aware dispatcher. When the report is infeasible the
advisor plan is drawn so the reader can see at a glance which methods
were considered and why none of them worked.

## References

Goh, B. S. (1977). Global stability in many-species systems. *The
American Naturalist*, 111(977), 135-143.
[doi:10.1086/283144](https://doi.org/10.1086/283144)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`print.lyapunov_report()`](https://robustecologies.github.io/janos/reference/print.lyapunov_report.md),
[`summary.lyapunov_report()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_report.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ x * (1 - x - 0.3 * y),
                           y ~ y * (0.8 - 0.2 * x - y)),
                state_names = c("x","y"), parms = list(),
                init = c(x = 0.5, y = 0.5))
plot(analyse_lyapunov(m, verbose = FALSE))
} # }
```
