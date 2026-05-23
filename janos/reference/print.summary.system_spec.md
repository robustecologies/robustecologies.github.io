# Print method for a summary.system_spec object

Formats a `summary.system_spec` object as a tabular report of states,
parameters and governing equations.

## Usage

``` r
# S3 method for class 'summary.system_spec'
print(x, ...)
```

## Arguments

- x:

  A `summary.system_spec` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.system_spec()`](https://robustecologies.github.io/janos/reference/summary.system_spec.md)
. summary constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
print(summary(m))
} # }
```
