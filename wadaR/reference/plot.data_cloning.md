# Plot data cloning MLE diagnostics

Produces diagnostic visualizations for `data_cloning` objects. The
default `"convergence"` type shows posterior mean and standard error as
a function of clone size K for each parameter.

## Usage

``` r
# S3 method for class 'data_cloning'
plot(x, type = "convergence", Parms = NULL, conf_level = NULL, col = NULL, ...)
```

## Arguments

- x:

  An object of class `data_cloning`.

- type:

  Character string specifying the plot type. One of `"convergence"`
  (default), `"eigenvalue"`, `"scaled_variance"`, `"density"`,
  `"trace"`, `"pairs"`, `"profile"`, `"estimability"`,
  `"prior_sensitivity"`, or `"diagnostics"`.

- Parms:

  Character string or regex for subsetting parameters. Default `NULL`
  (all parameters).

- conf_level:

  Numeric. Confidence level for profile likelihood. If `NULL`, uses the
  value from the `data_cloning` object.

- col:

  Optional color vector override.

- ...:

  Additional arguments. For `type = "profile"`, accepts `n_grid`
  (integer, default 50) and `Model`/`Data` which are required for
  profile computation.

## Value

Invisibly returns the ggplot object(s).

## Details

The ten plot types are:

- `"convergence"`:

  Posterior mean +/- SE vs K (log scale). The primary diagnostic for
  assessing MLE convergence.

- `"eigenvalue"`:

  Standardized largest eigenvalue ratio vs K with 1/K reference line.
  Departure from 1/K indicates non-identifiability.

- `"scaled_variance"`:

  log(K \* Var_K) vs log(K) per parameter. Horizontal indicates
  estimability.

- `"density"`:

  Overlaid posterior densities for each K. Higher K values are more
  opaque.

- `"trace"`:

  Trace plots from the final K run. Delegates to `plot.demonoid`.

- `"pairs"`:

  Bivariate posterior scatter from the max-K run with MLE overlay.

- `"profile"`:

  Profile likelihood curves computed on demand. Requires `Model` and
  `Data` in `...`.

- `"estimability"`:

  Per-parameter posterior means from different priors at each K
  (requires n_priors \> 1).

- `"prior_sensitivity"`:

  Forest plot of MLE +/- SE for each prior at the highest K (requires
  n_priors \> 1).

- `"diagnostics"`:

  Multi-panel dashboard. Row 1: convergence and eigenvalue. Row 2:
  scaled variance and density. When `n_priors > 1`, a third row is added
  with prior sensitivity and estimability (3x2 grid).

## References

Lele, S.R., Dennis, B. & Lutscher, F. (2007). Data cloning: easy maximum
likelihood estimation for complex ecological models using Bayesian
Markov chain Monte Carlo methods. *Ecology Letters*, 10, 551-563.
[doi:10.1111/j.1461-0248.2007.01047.x](https://doi.org/10.1111/j.1461-0248.2007.01047.x)

Lele, S.R., Nadeem, K. & Schmuland, B. (2010). Estimability and
likelihood inference for generalized linear mixed models using data
cloning. *J. Amer. Statist. Assoc.*, 105(492), 1617-1625.
[doi:10.1198/jasa.2010.tm09757](https://doi.org/10.1198/jasa.2010.tm09757)

Campbell, D. & Lele, S. (2014). An ANOVA test for parameter estimability
using data cloning with application to statistical inference for dynamic
systems. *Comput. Statist. Data Anal.*, 70, 257-267.
[doi:10.1016/j.csda.2013.09.013](https://doi.org/10.1016/j.csda.2013.09.013)

Ponciano, J.M., Taper, M.L., Dennis, B. & Lele, S.R. (2009).
Hierarchical models in ecology: confidence intervals, hypothesis
testing, and model selection using data cloning. *Ecology*, 90, 356-362.
[doi:10.1890/07-1960.1](https://doi.org/10.1890/07-1960.1)

Ponciano, J.M., Burleigh, J.G., Braun, E.L. & Taper, M.L. (2012).
Assessing parameter identifiability in phylogenetic models using data
cloning. *Syst. Biol.*, 61(6), 955-972.
[doi:10.1093/sysbio/sys055](https://doi.org/10.1093/sysbio/sys055)

## See also

[`DataCloning`](https://robustecologies.github.io/lucifer/reference/DataCloning.md),
[`summary.data_cloning`](https://robustecologies.github.io/lucifer/reference/summary.data_cloning.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dc <- DataCloning(Model, Data, Initial.Values = 0,
                  n_clones = c(1, 2, 4, 8))
plot(dc)
plot(dc, type = "eigenvalue")
plot(dc, type = "density")
plot(dc, type = "diagnostics")
} # }
```
