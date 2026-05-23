# MISS summary

This function summarizes posterior predictive distributions from an
object of class `miss`.

## Usage

``` r
# S3 method for class 'miss'
summary(object = NULL, ...)
```

## Arguments

- object:

  An object of class `miss` is required.

- ...:

  Additional arguments are unused.

## Value

This function returns a \\M \times 7\\ matrix, in which each row is the
posterior predictive distribution of one of \\M\\ missing values.
Columns are Mean, SD, MCSE, ESS, LB, Median, and UB.

## Details

This function summarizes the posterior predictive distributions from an
object of class `miss`.

## See also

[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model and data with missing values
N <- 25
y <- rnorm(N, 2, 0.5)
y[c(3, 10, 18)] <- NA
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE), na.rm = TRUE)
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

# Impute missing values
imp <- MISS(fit, Model, Data)

# Summarize imputations
summary(imp)
} # }
```
