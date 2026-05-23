# Plot the output of [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)

This may be used to plot, or save plots of, the iterated history of the
parameters and, if posterior samples were taken, density plots of
parameters and monitors in an object of class `laplace`. Alternative
types provide posterior density, interval, and pairs plots.

## Usage

``` r
# S3 method for class 'laplace'
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

  This required argument is an object of class `laplace`.

- Data:

  Optional. The list of data supplied to
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md).
  Monitor variable names are recovered automatically from the fit
  object.

- PDF:

  This logical argument indicates whether or not the user wants Lucifer
  to save the plots as a .pdf file.

- Parms:

  This argument accepts a vector of quoted strings to be matched for
  selecting parameters for plotting. This argument defaults to `NULL`
  and selects every parameter for plotting. Each quoted string is
  matched to one or more parameter names with the `grep` function. For
  example, if the user specifies `Parms=c("eta", "tau")`, and if the
  parameter names are beta\[1\], beta\[2\], eta\[1\], eta\[2\], and tau,
  then all parameters will be selected, because the string `eta` is
  within `beta`. Since `grep` is used, string matching uses regular
  expressions, so beware of meta-characters, though these are
  acceptable: ".", "\[", and "\]".

- type:

  Character string selecting the type of plot. One of `"diagnostics"`
  (default), `"posterior"`, `"intervals"`, or `"pairs"`.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided,
  vertical dashed lines or diamond markers are drawn at the true values.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

The default `"diagnostics"` type arranges plots in a \\2 \times 2\\
matrix. The purpose of the iterated history plots is to show how the
value of each parameter and the deviance changed by iteration as the
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
attempted to maximize the logarithm of the unnormalized joint posterior
density. If the algorithm converged, and if `sir=TRUE` in
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
then density plots are produced for selected parameters and all
monitored variables.

## See also

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
### See the LaplaceApproximation function for an example.
} # }
```
