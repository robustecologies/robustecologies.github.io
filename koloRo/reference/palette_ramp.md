# Create color ramp function from palette

Generates an interpolated sequence of colors from a named palette. This
is useful for creating smooth color gradients for continuous data
visualization or for generating more colors than the original palette
contains.

## Usage

``` r
palette_ramp(palette_name, n = 100, reverse = FALSE)
```

## Arguments

- palette_name:

  Character string specifying the palette name.

- n:

  Integer specifying the number of colors to generate (default 100).

- reverse:

  Logical. If TRUE, reverses the palette order (default FALSE).

## Value

A character vector of `n` hexadecimal color codes representing the
interpolated color ramp.

## Details

The function uses
[`colorRampPalette`](https://rdrr.io/r/grDevices/colorRamp.html) to
interpolate between the discrete colors of the original palette. For
perceptually uniform results with scientific data, prefer palettes from
the scientific category (e.g., viridis, batlow) which maintain uniform
perceptual spacing across interpolations.

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
# Create a 100-color gradient from viridis
gradient <- palette_ramp("viridis", n = 100)
length(gradient)

# Create a reversed gradient
gradient_rev <- palette_ramp("batlow", n = 50, reverse = TRUE)

# Visualize the ramp
plot_palette(palette_ramp("alhambra_nazari", n = 25))
} # }
```
