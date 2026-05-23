# Interpolate a palette to N colors (delegated to koloRo)

Returns a character vector of `n` hexadecimal colors interpolated from
the named palette. When koloRo is installed the call is forwarded to
[`koloRo::palette_ramp()`](https://rdrr.io/pkg/koloRo/man/palette_ramp.html);
otherwise only `"okabe_ito"` is supported via a local
[`colorRampPalette`](https://rdrr.io/r/grDevices/colorRamp.html)
fallback.

## Usage

``` r
palette_ramp(palette_name, n = 100, reverse = FALSE)
```

## Arguments

- palette_name:

  Character. Palette name from the koloRo catalogue.

- n:

  Integer. Number of colors to generate. Default 100.

- reverse:

  Logical. If `TRUE`, reverse the palette order before interpolation.
  Default `FALSE`.

## Value

A character vector of `n` hexadecimal color codes.

## References

Almaraz, P. (2026). *koloRo: a comprehensive color palette collection
for R*. <https://github.com/robustecologies/koloRo>.

## See also

[`palettes`](https://robustecologies.github.io/wadaR/reference/palettes.md)
for the underlying palette catalogue.

## Examples

``` r
# Built-in fallback (works without koloRo)
palette_ramp("okabe_ito", n = 5)
#> [1] "#E69F00" "#15A390" "#77AB7A" "#D26429" "#000000"
palette_ramp("okabe_ito", n = 5, reverse = TRUE)
#> [1] "#000000" "#D26429" "#78AB79" "#15A390" "#E69F00"

if (FALSE) { # \dontrun{
# Full catalogue (requires koloRo)
palette_ramp("viridis", n = 50)
} # }
```
