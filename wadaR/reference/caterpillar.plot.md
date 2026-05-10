# Caterpillar plot

A caterpillar plot is a horizontal plot of 3 quantiles of selected
distributions. This may be used to produce a caterpillar plot of
posterior samples (parameters and monitored variables) from an object
either of class `demonoid`, `iterquad`, `laplace`, `pmc`, `vb`, or a
matrix.

## Usage

``` r
caterpillar.plot(
  x,
  Parms = NULL,
  Title = NULL,
  ground_truth = NULL,
  col = NULL
)
```

## Arguments

- x:

  This required argument is an object of class `demonoid`, `iterquad`,
  `laplace`, `pmc`, `vb`, or a \\S \times J\\ matrix of \\S\\ samples
  and \\J\\ variables. For an object of class `demonoid`, the
  distributions of the stationary posterior summary (`Summary2`) will be
  attempted first, and if missing, then the parameters of all posterior
  samples (`Summary1`) will be plotted. For a multi-chain `demonoid`,
  each chain is plotted in a different color. For an object of class
  `laplace` or `vb`, the distributions in the posterior summary,
  `Summary`, are plotted according to the posterior draws, sampled with
  sampling importance resampling in the
  [`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md)
  function. When a generic matrix is supplied, unimodal 95% HPD
  intervals are estimated with the
  [`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)
  function.

- Parms:

  This argument accepts a vector of quoted strings to be matched for
  selecting parameters and monitored variables for plotting (though all
  parameters are selected when a generic matrix is supplied). This
  argument defaults to `NULL` and selects every parameter for plotting.
  Each quoted string is matched to one or more parameter names with the
  `grep` function. For example, if the user specifies
  `Parms=c("eta", "tau")`, and if the parameter names are beta\[1\],
  beta\[2\], eta\[1\], eta\[2\], and tau, then all parameters will be
  selected, because the string `eta` is within `beta`. Since `grep` is
  used, string matching uses regular expressions, so beware of
  meta-characters, though these are acceptable: ".", "\[", and "\]".

- Title:

  This argument accepts a title for the plot.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided,
  diamond markers are overlaid at the true values.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

## Value

See Details.

## Details

Caterpillar plots are popular plots in Bayesian inference for
summarizing the quantiles of posterior samples. A caterpillar plot is
similar to a horizontal boxplot, though without quartiles, making it
easier for the user to study more distributions in a single plot. The
following quantiles are plotted as a line for each parameter: 0.025 and
0.975, with the exception of a generic matrix, where unimodal 95% HPD
intervals are estimated (for more information, see
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)).
A vertical, gray line is included at zero. The median appears as a dot,
and the quantile line extends to the 2.5% and 97.5% quantiles. For
multi-chain `demonoid` objects, the color of the median and quantile
line differs by chain.

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md),
[`SIR`](https://robustecologies.github.io/lucifer/reference/SIR.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model and data
N <- 25
y <- rnorm(N, 2, 0.5)
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnormv(mu, 0, 1000, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  Monitor <- LP
  return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
    yhat = yhat, parm = parm))
}

# Fit model
fit <- lucifer(Model, Data, Initial.Values = c(0, 0),
  Iterations = 1000, Status = 200, Thinning = 1,
  Algorithm = "CHARM", Specs = NULL)

# Caterpillar plot of all parameters
caterpillar.plot(fit)

# Caterpillar plot of selected parameters
caterpillar.plot(fit, Parms = "mu")
} # }
```
