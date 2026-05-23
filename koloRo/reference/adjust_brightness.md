# Lighten or darken a color

Adjusts the brightness of one or more colors by blending toward white
(lightening) or black (darkening). This function is useful for creating
tints and shades of base colors for layered visualizations or 3D
effects.

## Usage

``` r
adjust_brightness(color, amount = 0.3)
```

## Arguments

- color:

  Hex color code or vector of hex color codes.

- amount:

  Numeric value specifying the adjustment amount, ranging from -1 to 1.
  Positive values lighten the color (blend toward white), negative
  values darken it (blend toward black). Default is 0.3.

## Value

A character vector of hex color codes with adjusted brightness, same
length as the input `color` vector.

## Details

The function converts input colors from hexadecimal representation to
their RGB channel components, then modifies each channel independently
to achieve the desired brightness shift. The HSV color model, introduced
by Smith (1978), separates chromatic information (hue and saturation)
from luminance (value), making brightness adjustment conceptually
independent of color identity; however, this implementation operates
directly in RGB space for efficiency. For lightening, each channel is
interpolated toward maximum intensity as
`new = old + amount * (1 - old)`, which preserves the relative channel
ratios and avoids hue shifts. For darkening, channels are scaled
proportionally as `new = old * (1 + amount)`, maintaining color balance
while reducing overall luminance. All resulting values are clamped to
the valid 0-1 range before conversion back to hexadecimal notation.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- "#1E4D8C"  # Alhambra blue

# Lighten the color
lighter <- adjust_brightness(base, 0.3)
lightest <- adjust_brightness(base, 0.6)

# Darken the color
darker <- adjust_brightness(base, -0.3)
darkest <- adjust_brightness(base, -0.6)

# Create a brightness scale
scale <- adjust_brightness(base, seq(-0.5, 0.5, by = 0.1))
plot_palette(scale)
} # }
```
