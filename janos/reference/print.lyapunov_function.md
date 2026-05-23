# Print method for a lyapunov_function object

Formats a `lyapunov_function` S3 object as a compact single-screen
report: the construction method, the envelope structural type and
dimension, the equilibrium point, the certification status for \\V \succ
0\\ and \\\dot V \prec 0\\, the numerical residual, the
domain-of-attraction category, any method-specific note and the elapsed
wall-clock time.

## Usage

``` r
# S3 method for class 'lyapunov_function'
print(x, ...)
```

## Arguments

- x:

  A `lyapunov_function` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## Details

The coloured glyphs indicate whether the positive-definiteness of \\V\\
and the negative-definiteness of \\\dot V\\ have been certified. A
missing \\V\\ (indicated by a red cross) means the construction failed;
the `certificate$method_details` field explains why and is surfaced in
the `Note:` line.

## References

Lyapunov, A. M. (1892). *The general problem of the stability of
motion*. Translation in: Int. J. Control, 55(3), 531-773 (1992).
[doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253)

## See also

[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
. constructor;
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md)
. numerical spot check;
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
. level sets and surfaces;
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
. higher-level family-aware dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
A <- diag(c(-1, -2))
print(lyapunov_function(A))
} # }
```
