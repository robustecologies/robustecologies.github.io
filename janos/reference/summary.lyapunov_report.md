# Summary method for a lyapunov_report object

Combines the advisor summary and the inner Lyapunov-function summary
into a single report, with the top-level family, subtype and reasoning
followed by the numerical diagnostics of the constructed Lyapunov
function (when available).

## Usage

``` r
# S3 method for class 'lyapunov_report'
summary(object, ...)
```

## Arguments

- object:

  A `lyapunov_report` object.

- ...:

  Passed to
  [`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md)
  when applicable.

## Value

A list of class `summary.lyapunov_report`, invisibly.

## Details

This method is the preferred way to inspect a full analysis: it surfaces
the advisor's decision tree, the rejected methods with their rationales,
and the spot-check results from the inner Lyapunov function in one
place.

## References

Khasminskii, R. (2012). *Stochastic Stability of Differential Equations*
(2nd ed.). Springer.
[doi:10.1007/978-3-642-23280-0](https://doi.org/10.1007/978-3-642-23280-0)

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`print.lyapunov_report()`](https://robustecologies.github.io/janos/reference/print.lyapunov_report.md),
[`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md),
[`summary.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_advisor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
summary(analyse_lyapunov(m, verbose = FALSE))
} # }
```
