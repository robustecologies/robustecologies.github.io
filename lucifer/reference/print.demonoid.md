# Print an object of class `demonoid` to the screen

This may be used to print the contents of an object of class `demonoid`
to the screen.

## Usage

``` r
# S3 method for class 'demonoid'
print(x, ...)
```

## Arguments

- x:

  An object of class `demonoid` is required.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. The side effect is the printed report.

## Details

For multi-chain runs (`Chains > 1`), `print.demonoid` additionally
displays the number of chains, CPUs used, per-chain acceptance rates,
and the maximum Gelman-Rubin PSRF.

## See also

[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md)
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

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

# Print results
print(fit)
} # }
```
