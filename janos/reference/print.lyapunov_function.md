# Print method for a lyapunov_function object

Formats a `lyapunov_function` S3 object as a short report with the
construction method, equilibrium point, certification status, residual
and domain of attraction estimate.

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
the `certificate$method_details` field explains why.

## References

Lyapunov, A. M. (1892). *The general problem of the stability of
motion*. Translation in: Int. J. Control, 55(3), 531-773 (1992).
[doi:10.1080/00207179208934253](https://doi.org/10.1080/00207179208934253)

## See also

[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md),
[`summary.lyapunov_function()`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md),
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md),
[`dynamical_system()`](https://robustecologies.github.io/janos/reference/dynamical_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
A <- diag(c(-1, -2))
print(lyapunov_function(A))
} # }
```
