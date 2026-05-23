# Mode estimation and modality testing

The mode is a measure of central tendency representing the value that
occurs most frequently in discrete data or has the highest density in
continuous data. `Mode` returns the single primary mode; `Modes` returns
all modes with their relative sizes. Both return objects of class
`mode_estimate` with associated print, summary, and plot methods. The
`is.*` functions test modality.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `bimodal`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `multimodal`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `trimodal`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `unimodal`.

See Details.

See Details.

## Usage

``` r
is.amodal(x, min.size = 0.1)

is.bimodal(x, min.size = 0.1)

is.multimodal(x, min.size = 0.1)

is.trimodal(x, min.size = 0.1)

is.unimodal(x, min.size = 0.1)

Mode(x, method = "auto")

Modes(x, min.size = 0.1, method = "kde")
```

## Arguments

- x:

  A numeric vector in which a mode (or modes) will be sought.

- min.size:

  The minimum relative size for a mode to be retained, where size is the
  proportion of density between adjacent antimodes. Default 0.1.

- method:

  Character string specifying the estimation method. `"auto"` selects
  `"kde"` for continuous data and table-based counting for discrete
  data. Options: `"kde"` (kernel density estimation), `"hsm"`
  (half-sample mode; Bickel & Fruhwirth 2006), `"venter"` (Venter's
  mode, shortest interval containing half the data), `"shorth"`
  (shortest half-range). Default `"auto"`.

## Value

`Mode` returns an object of class `mode_estimate` containing the primary
mode value. The object is backward-compatible: subscripting with `[1]`
and [`as.numeric()`](https://rdrr.io/r/base/numeric.html) return the
scalar mode value. `Modes` returns a `mode_estimate` list with
components `modes`, `mode.dens`, and `size` (accessible by name or
position). The `is.*` functions return logical scalars.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

`Mode` returns the single most prominent mode. `Modes` returns all
detected modes ordered by density, each with its relative size
(proportion of the density it captures). The `is.*` functions test
modality: `is.amodal` tests for no modes (constant input), `is.unimodal`
for exactly one mode, `is.bimodal` for two, `is.trimodal` for three, and
`is.multimodal` for more than one.

The half-sample mode (HSM) recursively selects the densest half of the
sample, converging on the mode location. It is robust to outliers and
does not require bandwidth selection. Venter's mode finds the shortest
interval containing 50% of the data and returns its midpoint. The
shortest half-range (shorth) is the minimum of the range of the shortest
half of the sorted sample.

Implementation of `is.bimodal`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.multimodal`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.trimodal`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.unimodal`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `Mode`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `Modes`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Bickel, D.R. and Fruhwirth, R. (2006). On a fast, robust estimator of
the mode: comparisons to other robust estimators with applications.
*Computational Statistics & Data Analysis*, 50(12), 3500-3530.
[doi:10.1016/j.csda.2005.07.014](https://doi.org/10.1016/j.csda.2005.07.014)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

[`as.double.mode_estimate`](https://robustecologies.github.io/lucifer/reference/as.double.mode_estimate.md),
[`plot.mode_estimate`](https://robustecologies.github.io/lucifer/reference/plot.mode_estimate.md),
[`print.mode_estimate`](https://robustecologies.github.io/lucifer/reference/print.mode_estimate.md),
[`summary.mode_estimate`](https://robustecologies.github.io/lucifer/reference/summary.mode_estimate.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Unimodal continuous data
x <- rnorm(1000, mean = 5)
m <- Mode(x)
print(m)
plot(m)

# Bimodal data
x <- c(rnorm(500, 0), rnorm(500, 6))
ms <- Modes(x)
print(ms)
summary(ms)
plot(ms)

# Compare estimation methods
plot(Modes(x), type = "comparison")

# Modality tests
is.unimodal(rnorm(1000))
is.bimodal(c(rnorm(500, 0), rnorm(500, 6)))
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.bimodal
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.multimodal
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.trimodal
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.unimodal
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Mode
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving Modes
} # }
```
