# Print method for a poincare_section object

Prints the section plane, the crossing direction and the number of
crossings detected on the attractor. Compact, one-screen header.

## Usage

``` r
# S3 method for class 'poincare_section'
print(x, ...)
```

## Arguments

- x:

  A `poincare_section` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
. constructor;
[`summary.poincare_section()`](https://robustecologies.github.io/janos/reference/summary.poincare_section.md)
. crossing statistics;
[`plot.poincare_section()`](https://robustecologies.github.io/janos/reference/plot.poincare_section.md)
. scatter of crossings on the section.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(rossler, t_max = 500, discard_transient = 100,
               solver = solver_rk45())
print(poincare_section(run, var = "z", value = 0.1))
} # }
```
