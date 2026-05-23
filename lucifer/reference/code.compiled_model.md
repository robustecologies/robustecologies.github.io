# Display or save generated C++ code from a compiled model

Prints the generated C++ source code with optional file output. When
`file` is provided, writes the code to disk and returns the path
invisibly; otherwise prints to the console.

## Usage

``` r
# S3 method for class 'compiled_model'
code(x, section = c("all", "eval", "grad", "data"), file = NULL, ...)
```

## Arguments

- x:

  A `compiled_model` object.

- section:

  Character; which section to display. One of `"all"` (default),
  `"eval"` (log-posterior only), `"grad"` (gradient only), or `"data"`
  (data struct only).

- file:

  Optional file path. If provided, writes the full C++ source to this
  file.

- ...:

  Ignored.

## Value

Invisibly returns the code as a character string (when `file = NULL`) or
the file path (when `file` is given).

## Details

Implementation of `code.compiled_model`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`print.compiled_model`](https://robustecologies.github.io/lucifer/reference/print.compiled_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
spec <- model_spec("y ~ Normal(mu, sigma)\nmu ~ Normal(0, 100)\nsigma ~ HalfCauchy(25)")
Data <- spec$data_template(y = rnorm(30))
compiled <- compile(spec, Data)

# Print to console
code(compiled)

# Show only the gradient function
code(compiled, section = "grad")

# Save to file
code(compiled, file = "my_model.cpp")
} # }
```
