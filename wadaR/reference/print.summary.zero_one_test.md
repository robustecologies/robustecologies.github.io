# Print method for a summary.zero_one_test object

Formats the 0-1 test summary as a short report with the median \\K\\,
its empirical 95% CI and IQR, the verdict distribution across random
frequencies, and the effective sample size after subsampling.

## Usage

``` r
# S3 method for class 'summary.zero_one_test'
print(x, ...)
```

## Arguments

- x:

  A `summary.zero_one_test` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.zero_one_test()`](https://robustecologies.github.io/janos/reference/summary.zero_one_test.md)
— summary constructor;
[`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
— upstream constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
print(summary(zero_one_test(run)))
} # }
```
