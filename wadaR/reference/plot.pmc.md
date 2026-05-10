# Plot samples from the output of PMC

This may be used to plot, or save plots of, samples in an object of
class `pmc`. The default type shows trace/density pairs for parameters,
deviance density, monitor densities, and convergence diagnostics.
Alternative types provide posterior density, convergence, mixture
probability, interval, and pairs plots.

## Usage

``` r
# S3 method for class 'pmc'
plot(
  x,
  BurnIn = 0,
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

  This required argument is an object of class `pmc`.

- BurnIn:

  This argument requires zero or a positive integer that indicates the
  number of iterations to discard as burn-in for the purposes of
  plotting.

- Data:

  Optional. The list of data supplied to
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).
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
  (default), `"posterior"`, `"convergence"`, `"mixture"`, `"intervals"`,
  or `"pairs"`.

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

The default type `"diagnostics"` arranges plots in a \\2 \times 2\\
matrix. Each row represents a parameter, the deviance, or a monitored
variable. For parameters, the left column displays trace plots and the
right column displays kernel density plots.

Following these plots are three plots for convergence. First, ESSN (red)
and perplexity (black) are plotted by iteration. The second plot shows
the distribution of the normalized importance weights by iteration. The
third plot appears only when multiple mixture components are used,
displaying the probabilities of each mixture component by iteration.

## See also

[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md) and
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### See the PMC function for an example.
} # }
```
