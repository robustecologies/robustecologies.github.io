# Summary method for a lyapunov_advisor object

Produces a verbose summary of a `lyapunov_advisor` object with the full
reasoning, the structural detection details (linearity probe, gLV
recovery residuals, symmetry ratio, polynomial degree, additive noise
test) and the theory applied for each feasible method.

## Usage

``` r
# S3 method for class 'lyapunov_advisor'
summary(object, ...)
```

## Arguments

- object:

  A `lyapunov_advisor` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

A list of class `summary.lyapunov_advisor` with `feasible`, `rejected`,
`theory`, `details`, invisibly.

## Details

This method is intended for interactive inspection when deciding whether
to trust the advisor or override it. Printing a summary shows every
rationale the dispatcher would otherwise consult.

## References

Goh, B. S. (1977). Global stability in many-species systems. *The
American Naturalist*, 111(977), 135-143.
[doi:10.1086/283144](https://doi.org/10.1086/283144)

## See also

[`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md),
[`print.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/print.lyapunov_advisor.md),
[`plot.lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
m <- system_spec(rhs = list(x ~ -x), state_names = "x",
                parms = list(), init = c(x = 1))
summary(lyapunov_advisor(m, verbose = FALSE))
} # }
```
