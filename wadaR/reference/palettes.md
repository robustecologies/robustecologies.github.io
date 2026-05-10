# Comprehensive color palette collection

Returns a comprehensive collection of color palettes from various
sources including scientific visualization standards, extensive
colorblind-safe options, chameleon species coloration patterns, artistic
movements, natural phenomena, authentic Alhambra geometric patterns with
historical pigments, and contemporary design trends.

## Usage

``` r
palettes(category = "all", palette = NULL)
```

## Arguments

- category:

  Character string specifying the palette category. Options include:

  - `"all"` - Returns all available palettes (default)

  - `"scientific"` - Scientific visualization palettes (viridis, plasma,
    etc.)

  - `"colorblind"` - Extensive colorblind-safe palettes (40+ options)

  - `"alhambra"` - Authentic Alhambra palettes with historical pigments
    (50+)

  - `"chameleons"` - Chameleon species-based palettes (20+ species)

  - `"natural"` - Nature-inspired palettes (ocean, forest, sunset, etc.)

  - `"cultural"` - Cultural and heritage palettes (Japanese, Persian,
    etc.)

  - `"artistic"` - Art movement palettes (impressionist, bauhaus, etc.)

  - `"seasonal"` - Seasonal palettes (spring, summer, autumn, winter)

  - `"modern"` - Contemporary design palettes (neon, pastel, cyberpunk,
    etc.)

  - `"classic"` - Classic color theory palettes

  - `"monochrome"` - Single-hue variations

  - `"food"` - Culinary-inspired palettes

  - `"diverging"` - Diverging palettes for data with neutral midpoint

  - `"qualitative"` - Qualitative palettes for categorical data

  - `"cinema"` - Film-inspired palettes (Dune cinematography)

- palette:

  Character string specifying a specific palette name. If provided,
  returns only that palette. Use
  [`list_palettes()`](https://robustecologies.github.io/koloRo/reference/list_palettes.md)
  to see all available names.

## Value

A named list of color palettes, where each palette is a character vector
of hexadecimal color codes.

## References

Itten, J. (1961). *The art of color*. Reinhold Publishing. ISBN
978-0471289296

CIE (2004). Colorimetry (3rd ed.). CIE Publication 15:2004. Commission
Internationale de l'Eclairage. ISBN 978-3901906336

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

Okabe, M., & Ito, K. (2008). Color universal design (CUD): how to make
figures and presentations that are friendly to colorblind people.
*J\*Fly Data Depository for Drosophila Researchers*.
<https://jfly.uni-koeln.de/color/>

Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd
ed.). Springer-Verlag.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all colorblind-safe palettes
cb_palettes <- palettes(category = "colorblind")

# Get authentic Alhambra palettes
alhambra_palettes <- palettes(category = "alhambra")

# Get a specific palette
nazari <- palettes(palette = "alhambra_nazari")
okabe <- palettes(palette = "okabe_ito")

# Create a color ramp
colors <- colorRampPalette(palettes(palette = "viridis"))(100)
} # }
```
