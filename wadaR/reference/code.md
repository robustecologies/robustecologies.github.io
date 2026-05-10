# Extract generated code from an object

Generic function to extract generated source code from objects that
produce code programmatically.

S3 method: apply `code()` to objects of class `model_spec`.

## Usage

``` r
code(x, ...)

# S3 method for class 'model_spec'
code(x, ...)
```

## Arguments

- x:

  See Details.

- ...:

  See Details.

## Value

A character string of formatted source code.

See Details.

## Details

Implementation of `code`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `code.model_spec`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## Methods (by class)

- `code(model_spec)`: Extract generated Model function code from a
  model_spec.

## See also

`code.model_spec`,
[`compile.model_spec`](https://robustecologies.github.io/lucifer/reference/compile.model_spec.md),
[`plot.model_spec`](https://robustecologies.github.io/lucifer/reference/plot.model_spec.md),
[`print.model_spec`](https://robustecologies.github.io/lucifer/reference/print.model_spec.md),
[`summary.model_spec`](https://robustecologies.github.io/lucifer/reference/summary.model_spec.md).

[`model_spec`](https://robustecologies.github.io/lucifer/reference/model_spec.md),
[`compile.model_spec`](https://robustecologies.github.io/lucifer/reference/compile.model_spec.md),
[`plot.model_spec`](https://robustecologies.github.io/lucifer/reference/plot.model_spec.md),
[`print.model_spec`](https://robustecologies.github.io/lucifer/reference/print.model_spec.md),
[`summary.model_spec`](https://robustecologies.github.io/lucifer/reference/summary.model_spec.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving code
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving code.model_spec
} # }
```
