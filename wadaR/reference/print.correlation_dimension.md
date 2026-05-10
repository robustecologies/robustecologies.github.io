# Print method for a correlation_dimension object

Prints the Grassberger-Procaccia \\D_2\\ estimate, the sample size used
for the correlation sum, the Theiler window for temporal decorrelation,
and the scaling window \\\[\varepsilon\_{\min}, \varepsilon\_{\max}\]\\
over which the power-law fit was performed.

## Usage

``` r
# S3 method for class 'correlation_dimension'
print(x, ...)
```

## Arguments

- x:

  A `correlation_dimension` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
— constructor;
[`summary.correlation_dimension()`](https://robustecologies.github.io/janos/reference/summary.correlation_dimension.md)
— fit diagnostics;
[`plot.correlation_dimension()`](https://robustecologies.github.io/janos/reference/plot.correlation_dimension.md)
— log-log scaling plot.

## Examples

``` r
if (FALSE) { # \dontrun{
run <- dyn_sim(lorenz, t_max = 200, discard_transient = 50)
print(correlation_dimension(run))
} # }
```
