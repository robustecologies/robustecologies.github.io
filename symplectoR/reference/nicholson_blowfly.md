# Nicholson blowfly cage population series

Adult *Lucilia cuprina* (Australian sheep blowfly) counts from the
laboratory cage experiments of Nicholson (1957), in long tidy form. The
`main` cage is the 180-point canonical series; `rep1`, `rep2` and `rep3`
are three further cages run for different durations. The series is the
textbook example of delayed density dependence producing sustained,
aperiodic population cycles, and its estimation problem is the package's
reference for the global-then-local solver pipeline.

## Usage

``` r
nicholson_blowfly
```

## Format

A data frame with 677 rows and 3 variables:

- population:

  Cage identity, a factor with levels `main`, `rep1`, `rep2` and `rep3`.

- day:

  Day of the experiment, sampled every two days (integer).

- count:

  Adult-fly count on that day (numeric).

## Source

The cage experiments of Nicholson (1957), tabulated for analysis by
Brillinger, Guckenheimer, Guttorp and Oster (1980) and distributed with
the analyses of Wood (2010); the tabulation enters this package through
the RElab ecosystem data store.

## Details

The mechanistic model is the delayed-recruitment equation of Gurney,
Blythe and Nisbet (1980), \$\$\dot N(t) = P\\ N(t - \tau)\\ e^{-N(t -
\tau)/N_0} - \delta N(t),\$\$ in which recruitment depends on the adult
population one development time \\\tau\\ in the past and mortality is
instantaneous. On the two-day sampling grid its deterministic skeleton
reads \$\$N\_{t+1} = P\\ N\_{t-k}\\ e^{-N\_{t-k}/N_0} + N_t
e^{-\delta},\$\$ with \\k = 7\\ the delay in sampling steps.

The estimator matters here. The dynamics are aperiodic and the record is
noisy, so matching a simulated trajectory to the whole 360-day series
compares phases that decohere long before the end of the record; this is
precisely the failure that motivates the synthetic-likelihood treatment
of Wood (2010). Conditional one-step-ahead prediction avoids it: each
predicted value is conditioned on the observed state, so no phase error
accumulates. The resulting three-parameter surface (\\\log P\\, \\\log
N_0\\, \\\log \delta\\) is exactly at the dimension ceiling of the
quantum Hamiltonian descent solver, which makes this dataset the natural
demonstration of a global scan followed by symplectic refinement.

## References

Nicholson, A. J. (1957). The self-adjustment of populations to change.
*Cold Spring Harbor Symposia on Quantitative Biology*, 22, 153-173.
[doi:10.1101/sqb.1957.022.01.017](https://doi.org/10.1101/sqb.1957.022.01.017)

Gurney, W. S. C., Blythe, S. P., & Nisbet, R. M. (1980). Nicholson's
blowflies revisited. *Nature*, 287(5777), 17-21.
[doi:10.1038/287017a0](https://doi.org/10.1038/287017a0)

Brillinger, D. R., Guckenheimer, J., Guttorp, P., & Oster, G. (1980).
Empirical modelling of population time series data: the case of age and
density dependent vital rates. *Lectures on Mathematics in the Life
Sciences*, 13, 65-90.

Wood, S. N. (2010). Statistical inference for noisy nonlinear ecological
dynamic systems. *Nature*, 466(7310), 1102-1104.
[doi:10.1038/nature09319](https://doi.org/10.1038/nature09319)

## See also

[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[paramecium_didinium](https://robustecologies.github.io/symplectoR/reference/paramecium_didinium.md)

## Examples

``` r
data(nicholson_blowfly)
str(nicholson_blowfly)
#> 'data.frame':    677 obs. of  3 variables:
#>  $ population: Factor w/ 4 levels "main","rep1",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ day       : int  2 4 6 8 10 12 14 16 18 20 ...
#>  $ count     : num  948 942 911 858 801 676 504 397 248 146 ...
table(nicholson_blowfly$population)
#> 
#> main rep1 rep2 rep3 
#>  180  275  111  111 

if (FALSE) { # \dontrun{
## Conditional one-step-ahead estimation of the delayed-recruitment model
main <- subset(nicholson_blowfly, population == "main")
N <- main$count
k <- 7L
idx <- seq.int(k + 1L, length(N) - 1L)
predict_step <- function(theta) {
  log(exp(theta[1]) * N[idx - k] * exp(-N[idx - k] / exp(theta[2])) +
      N[idx] * exp(-exp(theta[3])) + 1)
}
obj <- sym_inverse(predict_step, data = log(N[idx + 1L] + 1),
                   obs_times = main$day[idx + 1L],
                   theta_bounds = list(
                     lo = c(logP = log(0.05), logN0 = log(50),    logdelta = log(0.005)),
                     hi = c(logP = log(50),   logN0 = log(20000), logdelta = log(3))))

## Global scan, then symplectic refinement
global  <- sym_optim(obj, method = "qhd", seed = 1,
                     control = sym_control("qhd", N_grid = 48, K = 1200))
refined <- sym_optim(obj, x0 = global$x_best, method = "slc_expo",
                     control = sym_control("slc_expo", C = 0.1, h = 1.5,
                                           max_iter = 3000, tol_grad = 1e-9))
exp(refined$x_best)

## Observations against one-step-ahead predictions, and the composite view
plot(refined, type = "observed")
plot(refined, type = "residuals")
plot(refined, type = "dashboard")
} # }
```
