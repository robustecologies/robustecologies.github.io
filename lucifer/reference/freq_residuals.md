# Frequentist residual diagnostics from a lucifer fit

Re-evaluates the Model at the MLE, computes raw and standardized
residuals, and returns a diagnostic object with print and plot methods.
The plot method produces the four-panel display familiar from
`plot(lm(...))`: residuals vs fitted, normal Q-Q, scale-location, and
residuals vs leverage (Cook's distance).

## Usage

``` r
freq_residuals(object, Model, Data, n_calls = 50L)
```

## Arguments

- object:

  a fitted model object of class `laplace`, `iterquad`, `demonoid`, or
  `data_cloning`.

- Model:

  the Model function used for fitting.

- Data:

  the Data list used for fitting.

- n_calls:

  integer; number of Model evaluations to average for recovering the
  deterministic fitted values. Default 50.

## Value

An object of class `freq_residuals` with elements:

- fitted:

  numeric vector of fitted values (conditional mean)

- observed:

  numeric vector of observed responses

- residuals:

  numeric vector of raw residuals

- std.residuals:

  numeric vector of standardized residuals

- leverage:

  numeric vector of hat values

- cooks.distance:

  numeric vector of Cook's distances

- sigma:

  residual standard deviation estimate

- n:

  sample size

- p:

  number of parameters

## Details

For `laplace` and `iterquad` objects the mode is used as the MLE. For
`demonoid` objects the posterior mean of the stationary samples is used.
For `data_cloning` objects the stored MLE is used directly.

The Model is called once at the MLE to obtain `yhat` (the deterministic
fitted values). Because lucifer Models typically generate stochastic
`yhat` (via `rnorm`, `rbinom`, etc.), the function calls the Model
`n_calls` times and averages the `yhat` vectors to recover the
conditional mean. With the default `n_calls = 50` this is sufficient for
most models.

Raw residuals are \\e_i = y_i - \hat\mu_i\\. Standardized residuals
divide by \\\hat\sigma \sqrt{1 - h\_{ii}}\\ where \\h\_{ii}\\ is the
\$i\$-th diagonal element of the hat matrix \\H = X(X^TX)^{-1}X^T\\.
When no design matrix `X` is available in Data, leverage values are set
to \\p/n\\ and standardized residuals use the simpler formula \\e_i /
\hat\sigma\\.

## See also

[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md),
[`confint_compare`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- LaplaceApproximation(Model, Initial.Values, Data,
                             Iterations = 200, Method = "SPG",
                             CovEst = "Hessian", sir = FALSE)
rd <- freq_residuals(fit, Model, Data)
print(rd)
plot(rd)
} # }
```
