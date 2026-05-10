# Print method for a lyapunov_report object

Formats a `lyapunov_report` as a compact panel with the detected family,
the method used, the certificate type and the inner Lyapunov function
(when present). When the report is infeasible the rationale takes the
place of the inner function.

## Usage

``` r
# S3 method for class 'lyapunov_report'
print(x, ...)
```

## Arguments

- x:

  A `lyapunov_report` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## Details

`certificate_type` follows the epistemic classification used by
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md).
`"algebraic"` means that the proof is exact; `"local_algebraic"` means
that the proof holds only in a neighbourhood; `"numerical"` means the
evidence is a grid or sample check; `"none"` means no construction was
possible.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

## See also

[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md),
[`summary.lyapunov_report()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_report.md),
[`plot.lyapunov_report()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md),
[`print.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/print.lyapunov_advisor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
print(analyse_lyapunov(m, verbose = FALSE))
} # }
```
