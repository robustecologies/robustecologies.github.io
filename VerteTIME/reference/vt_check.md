# One-call smoke test of the installed VerteTIME package

Exercises every analytical family and every S3 lifecycle method against
the shipped
[vt_demo](https://robustecologies.github.io/VerteTIME/reference/vt_demo.md)
compilation and returns a tibble that records which calls succeeded. Use
the function as the single command that confirms an installed copy of
the package is healthy end-to-end, without depending on any
maintainer-private artefact.

## Usage

``` r
vt_check(verbose = TRUE, plots = TRUE)
```

## Arguments

- verbose:

  Logical, default `TRUE`. When `TRUE`, prints a one-screen summary of
  the call grid to the console as it runs.

- plots:

  Logical, default `TRUE`. When `TRUE`, every plot family is exercised
  (output saved into
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html)); when `FALSE`,
  plot calls are skipped (useful in headless CI without a graphics
  device).

## Value

A `tibble` with columns `family`, `call`, `status`
(`"ok"`/`"error"`/`"skip"`), `message`, `elapsed_s`. The result carries
an `"ok"` attribute that is `TRUE` when every executed row is `"ok"`.

## Details

The function loads `vt_demo` into a private environment, calls one
representative function per analytical family, and exercises `print`,
`summary` and `plot` on every resulting S3 object. Export functions are
exercised against [`tempfile()`](https://rdrr.io/r/base/tempfile.html)
paths. The returned tibble lists each call as a row, with the elapsed
time and the error message (when applicable).

The function is the canonical smoke test of the installed package: a
green tibble proves that `data(vt_demo)`, every metric, every S3 method,
every plot family, and the publication subsystem are all callable
without touching any maintainer-private folder.

## References

Wickham, H., & Bryan, J. (2023). *R Packages: Organize, Test, Document,
and Share Your Code* (2nd ed.). O'Reilly. <https://r-pkgs.org>.

## See also

[vt_demo](https://robustecologies.github.io/VerteTIME/reference/vt_demo.md),
[`vt_read()`](https://robustecologies.github.io/VerteTIME/reference/vt_read.md),
[vertetime](https://robustecologies.github.io/VerteTIME/reference/vertetime.md),
[`vt_validate()`](https://robustecologies.github.io/VerteTIME/reference/vt_validate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
v <- vt_check()
attr(v, "ok")
v[v$status != "ok", ]
} # }
```
