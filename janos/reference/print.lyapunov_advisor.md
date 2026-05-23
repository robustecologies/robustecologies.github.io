# Print method for a lyapunov_advisor object

Formats a `lyapunov_advisor` object as a compact report listing the
family, subtype, dimension, ranked feasible methods, rejected methods
with rationales, and key theory notes.

## Usage

``` r
# S3 method for class 'lyapunov_advisor'
print(x, ...)
```

## Arguments

- x:

  A `lyapunov_advisor` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## Details

The feasible list reflects a cascade from most informative (global) to
most generic (local linearisation). Rejections indicate why a theorem
does not apply; a missing rejection entry for a method means the advisor
did not even consider that method relevant to the detected family.

## References

Khalil, H. K. (2002). *Nonlinear Systems* (3rd ed.). Prentice Hall.
ISBN: 978-0-13-067389-3.

## See also

[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`summary.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_advisor.md),
[`plot.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md),
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
print(lyapunov_advisor(m, verbose = FALSE))
} # }
```
