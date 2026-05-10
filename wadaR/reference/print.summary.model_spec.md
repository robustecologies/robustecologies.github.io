# Print method for a summary.model_spec object

Formats a `summary.model_spec` object as a tabular report of states,
parameters and governing equations.

## Usage

``` r
# S3 method for class 'summary.model_spec'
print(x, ...)
```

## Arguments

- x:

  A `summary.model_spec` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.model_spec()`](https://robustecologies.github.io/janos/reference/summary.model_spec.md)
— summary constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
print(summary(m))
} # }
```
