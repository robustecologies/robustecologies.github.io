# Plot the output of [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

This may be used to plot, or save plots of, the iterated history of the
parameters and variances, and if posterior samples were taken, density
plots of parameters and monitors in an object of class `vb`. The `type`
argument selects among several visualization styles.

## Usage

``` r
# S3 method for class 'vb'
plot(
  x,
  Data = NULL,
  PDF = FALSE,
  Parms = NULL,
  type = "diagnostics",
  ground_truth = NULL,
  col = NULL,
  ...
)
```

## Arguments

- x:

  This required argument is an object of class `vb`.

- Data:

  Optional. The list of data supplied to
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).
  Monitor variable names are recovered automatically from the fit
  object, so `Data` is no longer required for standard plot types.

- PDF:

  This logical argument indicates whether or not the user wants Lucifer
  to save the plots as a .pdf file.

- Parms:

  This argument accepts a vector of quoted strings to be matched for
  selecting parameters for plotting.

- type:

  Character. Plot type: `"diagnostics"` (default), `"elbo"`,
  `"posterior"`, `"intervals"`, `"pairs"`, or `"pathfinder"` (ELBO along
  L-BFGS trajectory).

- ground_truth:

  Optional named numeric vector of true parameter values.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

When `type = "diagnostics"` (the default), the plots are arranged in a
\\3 \times 3\\ matrix showing parameter history, variance history, and
posterior density for each parameter. The `"elbo"` type shows the ELBO
convergence trajectory. The `"posterior"` type shows only posterior
densities with optional ground truth overlays. The `"intervals"` type
delegates to
[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md).

## See also

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
### See the VariationalBayes function for an example.
} # }
```
