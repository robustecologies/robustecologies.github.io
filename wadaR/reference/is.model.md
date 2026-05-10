# Logical check of a model

This function provides a logical test of whether or not a `Model`
specification function meets minimum requirements to be considered as
such.

## Usage

``` r
is.model(Model, Initial.Values, Data)
```

## Arguments

- Model:

  A model specification function. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

- Initial.Values:

  A vector of initial values, or current parameter values. For more
  information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

- Data:

  A list of data. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

## Value

A logical value of `TRUE` when `Model` meets minimum criteria of a model
specification function, and `FALSE` otherwise.

## Details

This function tests for minimum criteria for `Model` to be considered a
model specification function. Specifically, it tests:

- `Model` must be a function

- `Model` must execute without errors

- `Model` must return a list

- `Model` must have five components in the list

- The first component must be named LP and have length 1

- The second component must be named Dev and have length 1

- The third component must be named Monitor

- The lengths of Monitor and mon.names must be equal

- The fourth component must be named yhat

- The fifth component must be named parm

- The lengths of parm and parm.names must be equal

This function is not extensive, and checks only for these minimum
criteria. Additional checks are conducted in
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

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

# Validate the model specification
is.model(Model, Initial.Values = c(0, 0), Data) # TRUE
} # }
```
