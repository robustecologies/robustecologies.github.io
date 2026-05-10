# Print method for a summary.correlation_dimension object

Formats the correlation-dimension summary as a short report: point
estimate, confidence interval, \\R^2\\, scaling-window sample size and
the scaling range in log units.

## Usage

``` r
# S3 method for class 'summary.correlation_dimension'
print(x, ...)
```

## Arguments

- x:

  A `summary.correlation_dimension` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.correlation_dimension()`](https://robustecologies.github.io/janos/reference/summary.correlation_dimension.md)
— summary constructor;
[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
— upstream constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
print(summary(correlation_dimension(run)))
} # }
```
