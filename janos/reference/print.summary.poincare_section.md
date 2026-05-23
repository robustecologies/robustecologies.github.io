# Print method for a summary.poincare_section object

Formats a `summary.poincare_section` object as a tabular report of
crossing statistics. Called automatically when a
`summary.poincare_section` object is printed.

## Usage

``` r
# S3 method for class 'summary.poincare_section'
print(x, ...)
```

## Arguments

- x:

  A `summary.poincare_section` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.poincare_section()`](https://robustecologies.github.io/janos/reference/summary.poincare_section.md)
. summary constructor;
[`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
. upstream constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100)
print(summary(poincare_section(run, var = "z", value = 0.1)))
} # }
```
