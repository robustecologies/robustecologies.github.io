# Print method for SBI objects

Prints a concise summary of a simulation-based inference fit.

## Usage

``` r
# S3 method for class 'sbi'
print(x, ...)
```

## Arguments

- x:

  An object of class `sbi`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Details

Produces a concise one-screen console report of a simulation-based
inference fit produced by
[`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md).
Summary of the content is given below. The text is formatted for quick
visual triage in the terminal and intentionally elides large matrices
(covariance blocks, posterior samples) in favour of their structural
summaries (dimensions, diagonal, summary quantiles). For the full
numeric content, coerce the object with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or pass it to
[`summary`](https://rdrr.io/r/base/summary.html) or
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md).

## References

Cranmer, K., Brehmer, J., & Louppe, G. (2020). The frontier of
simulation-based inference. *PNAS*, 117(48), 30055-30062.
[doi:10.1073/pnas.1912789117](https://doi.org/10.1073/pnas.1912789117)

## See also

[`plot.sbi`](https://robustecologies.github.io/lucifer/reference/plot.sbi.md),
[`summary.sbi`](https://robustecologies.github.io/lucifer/reference/summary.sbi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving print.sbi
} # }
```
