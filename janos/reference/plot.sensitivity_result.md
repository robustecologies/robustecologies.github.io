# Plot method for sensitivity_result objects

Plot method for sensitivity_result objects

## Usage

``` r
# S3 method for class 'sensitivity_result'
plot(x, type = "gradient", ...)
```

## Arguments

- x:

  A `sensitivity_result` object

- type:

  Character; `"gradient"` (default) shows a bar chart of absolute
  gradient values, `"trajectory"` shows the adjoint variable lambda(t)
  over time.

- ...:

  Additional arguments (ignored)

## Value

A ggplot2 object.
