# Color palettes (delegated to koloRo)

Thin wrapper that exposes the full koloRo palette catalogue (282
palettes across scientific, colorblind-safe, alhambra, chameleons,
natural, cultural, artistic, seasonal and monochrome categories)
directly inside wadaR. When the koloRo package is installed the call is
forwarded verbatim to
[`koloRo::palettes()`](https://rdrr.io/pkg/koloRo/man/palettes.html);
when it is not, only the built-in fallback `"okabe_ito"` is available.

## Usage

``` r
palettes(category = "all", palette = NULL)
```

## Arguments

- category:

  Character. Palette category (forwarded to koloRo). Defaults to
  `"all"`.

- palette:

  Character. Specific palette name (forwarded to koloRo). If `NULL`,
  koloRo returns the full named list for the requested category.

## Value

A character vector of hexadecimal color codes when `palette` is set, or
the named list returned by
[`koloRo::palettes()`](https://rdrr.io/pkg/koloRo/man/palettes.html)
otherwise.

## Details

The koloRo dependency is declared in `Suggests`, so wadaR remains
installable on a system without koloRo and the bundled `wada_basins` /
`wada_straddle_method` / `wada_merging_method` plots continue to work
because they request only the colorblind-safe `"okabe_ito"` palette,
which is embedded as a fallback. Any other palette name triggers a
[`stop()`](https://rdrr.io/r/base/stop.html) with the install
instruction.

## References

Almaraz, P. (2026). *koloRo: a comprehensive color palette collection
for R*. <https://github.com/robustecologies/koloRo>.

Okabe, M., & Ito, K. (2008). Color universal design (CUD): how to make
figures and presentations that are friendly to colorblind people.
<https://jfly.uni-koeln.de/color/>.

## See also

[`palette_ramp`](https://robustecologies.github.io/wadaR/reference/palette_ramp.md)
for interpolating a palette to N colors.

## Examples

``` r
# Built-in fallback (works without koloRo)
palettes(palette = "okabe_ito")
#> [1] "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7"
#> [8] "#000000"

if (FALSE) { # \dontrun{
# Full catalogue (requires koloRo)
palettes(palette = "alhambra_nazari")
palettes(category = "scientific")
} # }
```
