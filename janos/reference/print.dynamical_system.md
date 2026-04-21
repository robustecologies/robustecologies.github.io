# Print method for a dynamical_system object

Formats a `dynamical_system` S3 object as a short report with the system
type, state dimension, variable names and any structured parameters.

## Usage

``` r
# S3 method for class 'dynamical_system'
print(x, ...)
```

## Arguments

- x:

  A `dynamical_system` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## Details

For `"glv"` systems the growth rate vector and the diagonal of the
interaction matrix are shown inline. For `"linear"` systems the spectral
abscissa is reported.

## References

Hofbauer, J., & Sigmund, K. (1998). *Evolutionary Games and Population
Dynamics*. Cambridge University Press.
[doi:10.1017/CBO9781139173179](https://doi.org/10.1017/CBO9781139173179)

## See also

[`dynamical_system()`](https://robustecologies.github.io/janos/reference/dynamical_system.md),
[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
r <- c(1, 0.5); alpha <- matrix(c(-1, -0.3, -0.2, -1), 2, 2)
print(dynamical_system(r = r, alpha = alpha, type = "glv"))
} # }
```
