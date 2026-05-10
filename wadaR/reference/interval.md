# Constrain to interval

This function constrains the value(s) of a scalar, vector, matrix, or
array to a specified interval \\\[a,b\]\\. In Bayesian inference it is
often used to truncate a parameter to an interval, such as \\p(\theta)
\in \[a,b\]\\. The `interval` function is often used in conjunction with
[`dtrunc`](https://robustecologies.github.io/lucifer/reference/dist.Truncated.md)
to truncate the prior probability distribution associated with the
constrained parameter. This is unrelated to the probability interval
(see
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md)
and
[`LPL.interval`](https://robustecologies.github.io/lucifer/reference/LPL.interval.md)).

## Usage

``` r
interval(x, a = -Inf, b = Inf, reflect = TRUE)
```

## Arguments

- x:

  a scalar, vector, matrix, or array whose elements will be constrained
  to the interval \\\[a,b\]\\.

- a:

  the lower bound of the interval. Defaults to `-Inf`.

- b:

  the upper bound of the interval. Defaults to `Inf`.

- reflect:

  logical; when `TRUE` (the default), a value outside the interval is
  reflected back into the interval. When `FALSE`, a value outside the
  interval is assigned the nearest boundary.

## Value

A scalar, vector, matrix, or array matching the input `x`, with each
element constrained to the interval \\\[a,b\]\\.

## Details

The default method reflects an out-of-bounds proposal off of the
boundaries until the proposal is within the specified interval. This is
rare in the literature but works well in practice. The alternative
method (`reflect = FALSE`) sets the value equal to the violated
boundary, which is not generally recommended.

After constraining a parameter, the constrained value should be updated
back into the `parm` vector so the algorithm knows it has been
constrained. Alternatives include log-transforming positive-only
parameters or discarding out-of-bounds samples after the model update.

## See also

[`dtrunc`](https://robustecologies.github.io/lucifer/reference/dist.Truncated.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LPL.interval`](https://robustecologies.github.io/lucifer/reference/LPL.interval.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`p.interval`](https://robustecologies.github.io/lucifer/reference/p.interval.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
x <- 2
interval(x, 0, 1)
X <- matrix(runif(25, -2, 2), 5, 5)
interval(X, -1, 1)
} # }
```
