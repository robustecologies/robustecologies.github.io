# Is ssm model

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `ssm_model`.

These functions each provide a logical test of the class of an object.

Functions in Lucifer often assign a class to an output object. For
example, after updating a model with the
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
function, an output object is created of class `demonoid`. Likewise,
after passing a model to the
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
function, an output object is created, and it is of class `laplace`. The
class of these and other objects may be logically tested.

By assigning a class to an output object, the package is able to discern
which other functions are appropriate for it. For example, after
updating a model with
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
which creates an object of class `demonoid`, the user may desire to plot
its output. Since it is assigned a class, the user may use the generic
`plot` function, which internally selects the
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md)
function, which differs from
[`plot.laplace`](https://robustecologies.github.io/lucifer/reference/plot.laplace.md)
for objects of class `laplace`.

For more information on object classes, see the
[`class`](https://rdrr.io/r/base/class.html) function.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `blocks`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `bmk`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `demonoid`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `demonoid.ppc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `demonoid.val`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `hangartner`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `heidelberger`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `importance`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `iterquad`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `iterquad.ppc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `laplace`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `laplace.ppc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `miss`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `pmc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `pmc.ppc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `pmc.val`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `posteriorchecks`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `raftery`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `rejection`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `vb`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `vb.ppc`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `prescription`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `arena`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `consort`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `crucible`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `sde_model`.

S3 method: apply [`is()`](https://rdrr.io/r/methods/is.html) to objects
of class `sde_fit`.

## Usage

``` r
is.ssm_model(x)

is.bayesfactor(x)

is.blocks(x)

is.bmk(x)

is.demonoid(x)

is.demonoid.ppc(x)

is.demonoid.val(x)

is.hangartner(x)

is.heidelberger(x)

is.importance(x)

is.iterquad(x)

is.iterquad.ppc(x)

is.laplace(x)

is.laplace.ppc(x)

is.miss(x)

is.pmc(x)

is.pmc.ppc(x)

is.pmc.val(x)

is.posteriorchecks(x)

is.raftery(x)

is.rejection(x)

is.vb(x)

is.vb.ppc(x)

is.prescription(x)

is.arena(x)

is.consort(x)

is.crucible(x)

is.sde_model(x)

is.sde_fit(x)
```

## Arguments

- x:

  An object that will be subjected to a logical test of its class.

## Value

See Details.

Each function returns a logical value indicating whether or not the
supplied object is of the corresponding class.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

See Details.

## Details

Implementation of `is.ssm_model`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.bayesfactor`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.blocks`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.bmk`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `is.demonoid`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.demonoid.ppc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.demonoid.val`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.hangartner`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.heidelberger`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.importance`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.iterquad`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.iterquad.ppc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.laplace`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.laplace.ppc`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.miss`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.pmc`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `is.pmc.ppc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.pmc.val`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.posteriorchecks`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.raftery`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.rejection`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.vb`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `is.vb.ppc`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.prescription`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

Implementation of `is.arena`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.consort`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.crucible`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.sde_model`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

Implementation of `is.sde_fit`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`print.ssm_model`](https://robustecologies.github.io/lucifer/reference/print.ssm_model.md).

[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`Blocks`](https://robustecologies.github.io/lucifer/reference/Blocks.md),
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`class`](https://rdrr.io/r/base/class.html),
[`Hangartner.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Hangartner.Diagnostic.md),
[`Heidelberger.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Heidelberger.Diagnostic.md),
[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`Raftery.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Raftery.Diagnostic.md),
[`RejectionSampling`](https://robustecologies.github.io/lucifer/reference/RejectionSampling.md),
[`Validate`](https://robustecologies.github.io/lucifer/reference/Validate.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`plot.bmk`](https://robustecologies.github.io/lucifer/reference/plot.bmk.md).

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LOO.demonoid`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.demonoid`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.demonoid`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.demonoid`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`print.demonoid`](https://robustecologies.github.io/lucifer/reference/print.demonoid.md).

[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md),
[`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md).

[`Heidelberger.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Heidelberger.Diagnostic.md),
[`print.heidelberger`](https://robustecologies.github.io/lucifer/reference/print.heidelberger.md).

[`Importance`](https://robustecologies.github.io/lucifer/reference/Importance.md),
[`plot.importance`](https://robustecologies.github.io/lucifer/reference/plot.importance.md).

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`plot.iterquad`](https://robustecologies.github.io/lucifer/reference/plot.iterquad.md),
[`predict.iterquad`](https://robustecologies.github.io/lucifer/reference/predict.iterquad.md),
[`print.iterquad`](https://robustecologies.github.io/lucifer/reference/print.iterquad.md).

[`plot.iterquad.ppc`](https://robustecologies.github.io/lucifer/reference/plot.iterquad.ppc.md),
[`summary.iterquad.ppc`](https://robustecologies.github.io/lucifer/reference/summary.iterquad.ppc.md).

[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`plot.laplace`](https://robustecologies.github.io/lucifer/reference/plot.laplace.md),
[`predict.laplace`](https://robustecologies.github.io/lucifer/reference/predict.laplace.md),
[`print.laplace`](https://robustecologies.github.io/lucifer/reference/print.laplace.md).

[`plot.laplace.ppc`](https://robustecologies.github.io/lucifer/reference/plot.laplace.ppc.md),
[`summary.laplace.ppc`](https://robustecologies.github.io/lucifer/reference/summary.laplace.ppc.md).

[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md),
[`plot.miss`](https://robustecologies.github.io/lucifer/reference/plot.miss.md),
[`print.miss`](https://robustecologies.github.io/lucifer/reference/print.miss.md),
[`summary.miss`](https://robustecologies.github.io/lucifer/reference/summary.miss.md).

[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
[`LOO.pmc`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.pmc`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.pmc`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.pmc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.md),
[`predict.pmc`](https://robustecologies.github.io/lucifer/reference/predict.pmc.md),
[`print.pmc`](https://robustecologies.github.io/lucifer/reference/print.pmc.md).

[`plot.pmc.ppc`](https://robustecologies.github.io/lucifer/reference/plot.pmc.ppc.md),
[`summary.pmc.ppc`](https://robustecologies.github.io/lucifer/reference/summary.pmc.ppc.md).

[`Raftery.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Raftery.Diagnostic.md),
[`print.raftery`](https://robustecologies.github.io/lucifer/reference/print.raftery.md).

[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`LOO.vb`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.vb`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.vb`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.vb`](https://robustecologies.github.io/lucifer/reference/plot.vb.md),
[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md),
[`print.vb`](https://robustecologies.github.io/lucifer/reference/print.vb.md).

[`plot.vb.ppc`](https://robustecologies.github.io/lucifer/reference/plot.vb.ppc.md),
[`summary.vb.ppc`](https://robustecologies.github.io/lucifer/reference/summary.vb.ppc.md).

[`plot.prescription`](https://robustecologies.github.io/lucifer/reference/plot.prescription.md),
[`print.prescription`](https://robustecologies.github.io/lucifer/reference/print.prescription.md),
[`summary.prescription`](https://robustecologies.github.io/lucifer/reference/summary.prescription.md).

[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`plot.arena`](https://robustecologies.github.io/lucifer/reference/plot.arena.md),
[`print.arena`](https://robustecologies.github.io/lucifer/reference/print.arena.md),
[`summary.arena`](https://robustecologies.github.io/lucifer/reference/summary.arena.md).

[`plot.consort`](https://robustecologies.github.io/lucifer/reference/plot.consort.md),
[`print.consort`](https://robustecologies.github.io/lucifer/reference/print.consort.md),
[`summary.consort`](https://robustecologies.github.io/lucifer/reference/summary.consort.md).

[`Crucible`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
[`plot.crucible`](https://robustecologies.github.io/lucifer/reference/plot.crucible.md),
[`print.crucible`](https://robustecologies.github.io/lucifer/reference/print.crucible.md),
[`summary.crucible`](https://robustecologies.github.io/lucifer/reference/summary.crucible.md).

[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
[`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md),
[`plot.sde_model`](https://robustecologies.github.io/lucifer/reference/plot.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`simulate.sde_model`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md),
[`summary.sde_model`](https://robustecologies.github.io/lucifer/reference/summary.sde_model.md).

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.ssm_model
} # }
if (FALSE) { # \dontrun{
# Check class of a demonoid fit object
is.demonoid(fit)   # TRUE if fit is of class demonoid
is.laplace(fit)    # FALSE

# Check other classes
is.bayesfactor(bf) # TRUE if bf is of class bayesfactor
is.pmc(pmc_fit)    # TRUE if pmc_fit is of class pmc
is.vb(vb_fit)      # TRUE if vb_fit is of class vb
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.blocks
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.bmk
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.demonoid
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.demonoid.ppc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.demonoid.val
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.hangartner
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.heidelberger
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.importance
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.iterquad
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.iterquad.ppc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.laplace
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.laplace.ppc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.miss
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.pmc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.pmc.ppc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.pmc.val
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.posteriorchecks
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.raftery
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.rejection
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.vb
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.vb.ppc
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.prescription
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.arena
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.consort
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.crucible
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.sde_model
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving is.sde_fit
} # }
```
