# Plot method for a lyapunov_advisor object

Renders one of three visual summaries of the advisor state:
`"family_tree"` (default) shows the decision tree with family, subtype,
feasible and rejected methods; `"detector_scores"` shows a horizontal
bar plot with the score of every structural detector that ran; `"radar"`
shows a radial confidence polygon across the seven supported families.

## Usage

``` r
# S3 method for class 'lyapunov_advisor'
plot(x, type = c("family_tree", "detector_scores", "radar"), ...)
```

## Arguments

- x:

  A `lyapunov_advisor` object.

- type:

  Character, one of `"family_tree"` (default), `"detector_scores"` or
  `"radar"`.

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object.

## Details

`"family_tree"` is the canonical view: rounded boxes connected by line
segments, with feasible methods filled in RElab orange and rejected
methods greyed out with the rejection reason inline. `"detector_scores"`
reads the internal detector results stored in `x$details` (linearity
probe, gLV recovery, polynomial degree, symmetry, additive noise) and
plots each score against the acceptance threshold 0.5, with the
numerical evidence annotated to the right of each bar. `"radar"` is
useful when the advisor's family classification is the question of
interest: the unit-radius vertex marks the detected family, and
alternative families sit at a floor radius unless a detector explicitly
supports them.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

## See also

[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`print.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/print.lyapunov_advisor.md),
[`summary.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_advisor.md),
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(rhs = list(x ~ x * (1 - x - 0.3 * y),
                           y ~ y * (0.8 - 0.2 * x - y)),
                state_names = c("x","y"), parms = list(),
                init = c(x = 0.5, y = 0.5))
plot(lyapunov_advisor(m, verbose = FALSE))
plot(lyapunov_advisor(m, verbose = FALSE), type = "detector_scores")
plot(lyapunov_advisor(m, verbose = FALSE), type = "radar")
} # }
```
