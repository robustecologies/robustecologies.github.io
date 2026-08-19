# Paramecium and Didinium microcosm predator-prey experiment

Abundances of the prey *Paramecium aurelia* and its predator *Didinium
nasutum* in a laboratory microcosm at 0.5 Cerophyll concentration,
recorded roughly every twelve hours over thirty-five days. The series is
among the cleanest predator-prey records in ecology, because it comes
from a controlled culture rather than from a field census, and it is the
reference dataset of the package's trajectory-matching workflow.

## Usage

``` r
paramecium_didinium
```

## Format

A data frame with 71 rows and 3 variables:

- time:

  Observation time in days since the start of the experiment, nominally
  at half-day intervals (numeric).

- paramecium:

  *Paramecium aurelia* abundance, individuals per millilitre (numeric).

- didinium:

  *Didinium nasutum* abundance, individuals per millilitre (numeric).

## Source

The microcosm experiment of Veilleux (1979), digitised by Jost and
Ellner (2000); the digitised record enters this package through the
RElab ecosystem data store. The published record contains one time
inversion between the thirty-third and thirty-fourth observations, an
artefact of reading coordinates off the printed figure; the rows are
ordered by time here so that the series is directly usable for
trajectory matching. No value is altered, interpolated or dropped by
that reordering.

## Details

The canonical use is trajectory matching of a two-species dynamical
system through
[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md):
build a
[`janos::system_spec()`](https://robustecologies.github.io/janos/reference/system_spec.html)
for the candidate model, hand it to
[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md)
with `theta_bounds`, and minimise the resulting objective with any
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md)
method. Estimating the rate parameters on the logarithmic scale is
strongly recommended, because their natural scales span several orders
of magnitude and the resulting Hessian is severely ill conditioned
otherwise.

Two model families are worth contrasting on this record. The classical
Lotka-Volterra system, \\\dot P = r P - a P D\\ and \\\dot D = b a P D -
m D\\, has four parameters and a neutrally stable centre; the
Rosenzweig-MacArthur system adds prey self-limitation and a saturating
functional response and has six. The comparison of their attained losses
is a small model-selection exercise that the package's solvers perform
in a few seconds.

## References

Veilleux, B. G. (1979). An analysis of the predatory interaction between
*Paramecium* and *Didinium*. *The Journal of Animal Ecology*, 48(3),
787-803. [doi:10.2307/4195](https://doi.org/10.2307/4195)

Jost, C., & Ellner, S. P. (2000). Testing for predator dependence in
predator-prey dynamics: a non-parametric approach. *Proceedings of the
Royal Society B: Biological Sciences*, 267(1453), 1611-1620.
[doi:10.1098/rspb.2000.1186](https://doi.org/10.1098/rspb.2000.1186)

Rosenzweig, M. L., & MacArthur, R. H. (1963). Graphical representation
and stability conditions of predator-prey interactions. *The American
Naturalist*, 97(895), 209-223.
[doi:10.1086/282272](https://doi.org/10.1086/282272)

## See also

[`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md),
[`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md),
[nicholson_blowfly](https://robustecologies.github.io/symplectoR/reference/nicholson_blowfly.md)

## Examples

``` r
data(paramecium_didinium)
head(paramecium_didinium)
#>   time paramecium didinium
#> 1 0.00      15.65     5.76
#> 2 0.54      53.57     9.05
#> 3 1.02      73.34    17.26
#> 4 1.54      93.93    41.97
#> 5 2.02     115.40    55.97
#> 6 2.48      76.57    74.91

if (FALSE) { # \dontrun{
## Trajectory matching of a Lotka-Volterra system, estimated on the log scale
lv <- janos::system_spec(
  rhs = list(paramecium ~ exp(lr) * paramecium - exp(la) * paramecium * didinium,
             didinium   ~ exp(lb) * exp(la) * paramecium * didinium - exp(lm) * didinium),
  state_names = c("paramecium", "didinium"),
  parms = list(lr = 0, la = -4, lb = -0.7, lm = 0),
  init = c(paramecium = paramecium_didinium$paramecium[1],
           didinium   = paramecium_didinium$didinium[1])
)
obj <- sym_inverse(lv, data = paramecium_didinium,
                   loss = "custom",
                   loss_fn = function(pred, data) sum((log(pred + 1) - log(data + 1))^2),
                   theta_bounds = list(lo = c(lr = -3, la = -9, lb = -5, lm = -3),
                                       hi = c(lr =  3, la =  0, lb =  2, lm =  3)))
fit <- sym_optim(obj, x0 = c(lr = 0, la = -4, lb = -0.7, lm = 0), method = "slc_expo")
summary(fit)

## The fit remembers the record it was estimated from, so it plots against it
plot(fit, type = "observed")     # both species, observations and fitted trajectory
plot(fit, type = "residuals")
plot(fit, type = "dashboard")    # the two data views plus the solver diagnostics
} # }
```
