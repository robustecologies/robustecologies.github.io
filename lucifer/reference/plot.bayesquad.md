# Plot the output of [`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)

Produces diagnostic plots for a Bayesian quadrature fit. The default
type shows the GP surrogate, evaluation point locations, integral
posterior convergence, and parameter density estimates. Alternative
types provide posterior density, interval, and pairs plots.

## Usage

``` r
# S3 method for class 'bayesquad'
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

  An object of class `bayesquad`.

- Data:

  Optional. The data list supplied to `BayesianQuadrature`. Monitor
  variable names are recovered automatically from the fit object.

- PDF:

  Logical; save plots to a PDF file (default `FALSE`).

- Parms:

  Character vector of parameter name patterns to select for plotting
  (uses `grep` matching).

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

The default `"diagnostics"` type produces four types of plot:

For each selected parameter: (1) a history plot showing the parameter
value at each evaluation point, and (2) a density overlay showing the
reference measure, GP-weighted density, and SIR posterior density (if
available).

Additionally: (3) the integral posterior convergence, showing the
posterior mean and 95% credible band as a function of the number of
evaluation points, and (4) a deviance history.

## See also

[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## See the BayesianQuadrature function for an example.
} # }
```
