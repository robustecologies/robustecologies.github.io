# Logical check of constraints

This function provides a logical test of constraints for each initial
value or parameter for a model specification, given data.

## Usage

``` r
is.constrained(Model, Initial.Values, Data)
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

A logical vector, equal in length to the number of initial values. Each
element receives `TRUE` if the corresponding initial value changed due
to a constraint, or `FALSE` if it did not.

## Details

This function is useful for testing whether or not initial values
changed due to constraints when being passed through a `Model`
specification function. If any initial value changes, then the
constrained values that are output in the fifth component of the `Model`
specification are suitable as initial values, not the tested initial
values.

A parameter may be constrained and this function may not discover the
constraint, since the discovery depends on the initial values and
whether or not they change as they are passed through the model.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model with a constrained parameter (log-transformed sigma)
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

# Check if initial values are constrained
is.constrained(Model, Initial.Values = c(0, 0), Data)
} # }
```
