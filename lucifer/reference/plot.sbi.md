# Plot method for SBI objects

Produces diagnostic plots for a simulation-based inference fit. The
default type shows marginal posterior density plots. Additional types
provide pairwise posterior plots, training curves, and calibration
diagnostics.

## Usage

``` r
# S3 method for class 'sbi'
plot(
  x,
  type = "posterior",
  ground_truth = NULL,
  true_values = NULL,
  sbc = NULL,
  coverage = NULL,
  ppc_data = NULL,
  tarp = NULL,
  c2st = NULL,
  col = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `sbi`.

- type:

  Character: `"posterior"` (default), `"pairs"`, `"training"`, `"sbc"`,
  `"coverage"`, `"ppc"`, `"tarp"`, or `"c2st"`.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided
  and `type = "posterior"`, vertical dashed lines are drawn at the true
  values.

- true_values:

  Alias for `ground_truth` (backward compatibility).

- sbc:

  An `sbi_sbc` object from
  [`SBC`](https://robustecologies.github.io/lucifer/reference/SBC.md)
  (required for `type = "sbc"`).

- coverage:

  An `sbi_coverage` object from
  [`expected_coverage`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md)
  (required for `type = "coverage"`).

- ppc_data:

  A list from
  [`sbi_ppc`](https://robustecologies.github.io/lucifer/reference/sbi_ppc.md)
  (for `type = "ppc"`).

- tarp:

  An `sbi_tarp` object from
  [`TARP`](https://robustecologies.github.io/lucifer/reference/TARP.md)
  (required for `type = "tarp"`).

- c2st:

  An `sbi_c2st` object from
  [`C2ST`](https://robustecologies.github.io/lucifer/reference/C2ST.md)
  (required for `type = "c2st"`).

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns `NULL`.

## Details

Produces diagnostic plots of a simulation-based inference fit produced
by [`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Cranmer, K., Brehmer, J., & Louppe, G. (2020). The frontier of
simulation-based inference. *PNAS*, 117(48), 30055-30062.
[doi:10.1073/pnas.1912789117](https://doi.org/10.1073/pnas.1912789117)

## See also

[`print.sbi`](https://robustecologies.github.io/lucifer/reference/print.sbi.md),
[`summary.sbi`](https://robustecologies.github.io/lucifer/reference/summary.sbi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.sbi
} # }
```
